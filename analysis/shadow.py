"""Earth's shadow at the orbit's radius, seen from Palomar, on each plate's grid: how much of a plate's exposure each
location spends in shadow, and the shadow's edge on the sky."""
import warnings
from dataclasses import dataclass, replace
from functools import lru_cache

import numpy as np
from scipy.ndimage import gaussian_filter
import astropy.units as u
from astropy.coordinates import EarthLocation, get_sun
from erfa import ErfaWarning
from astropy.utils.exceptions import AstropyWarning
from astropy.utils import iers

from geometry import deproject, unit                                            # noqa: F401  (used here and re-exported)

iers.conf.auto_download = False     # the Sun's position needs no fresh Earth-rotation tables; downloading them stalls runs

warnings.simplefilter("ignore", ErfaWarning)
warnings.simplefilter("ignore", AstropyWarning)

PALOMAR = EarthLocation.from_geodetic(lon=-116.8650 * u.deg, lat=33.3563 * u.deg, height=1706 * u.m)
R_GEO = 42164.0
# geocentric half-angle of Earth's shadow at geostationary radius, asin(6371 / 42164): the 8.69 deg of Villarroel et al.
# (2025, PASP, Table 4)
SHADOW_DEG = 8.69


def distinct(t):
    """The distinct times in t, and where each element of t sits among them (the Sun is then found once per time)."""
    if t.isscalar:
        return t, None
    _, i, inv = np.unique(np.c_[t.jd1, t.jd2], axis=0, return_index=True, return_inverse=True)
    return t[i], inv.ravel()


def antisun_axis(t):
    """Unit vector of the anti-solar direction (GCRS: ICRS-aligned, J2000 axes)."""
    tu, inv = distinct(t)
    s = get_sun(tu)
    a = -unit(s.ra.deg, s.dec.deg)
    return a if inv is None else a[inv]


# ------------------------------------------------------------------ plate grid for the within-plate tests
# Each plate is gridded into 0.1 deg cells over the central 6 x 6 deg; the shadow share of a cell is taken from 3 x 3
# points inside it.
HALF, CELL = 3.0, 0.1
g = np.arange(-HALF + CELL / 2, HALF, CELL)
GX, GY = [a.ravel() for a in np.meshgrid(g, g)]
SUB = np.array([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]) * CELL / 3


def cell_counts(xi, eta, w=None):
    """Transients per grid cell from plate-frame offsets (deg); with w, their summed weights."""
    ix = np.floor((np.asarray(xi) + HALF) / CELL).astype(int); iy = np.floor((np.asarray(eta) + HALF) / CELL).astype(int)
    ok = (ix >= 0) & (ix < len(g)) & (iy >= 0) & (iy < len(g))
    return np.bincount(iy[ok] * len(g) + ix[ok], weights=None if w is None else np.asarray(w)[ok], minlength=GX.size)



SIGMAS = (0.5, 1, 1.5, 2, 3, 4, 6)                                           # Gaussian widths tried, in 0.1 deg cells


def smooth_map(hosts):
    """The hosts' average map (hosts: the transients' tables of the plates no shadow reaches), mean 1, smoothed by the
    width from SIGMAS that best predicts held-out hosts: ten folds by plate number, each fold's transients scored by the
    mean log density (by weight) of the map from the other nine. Returns the map, the width, and every width's score."""
    M = np.array([cell_counts(g.xi.values, g.eta.values, g.w.values) / g.w.sum() for g in hosts])
    n = int(round(M.shape[1] ** 0.5)); fold = np.array([int(g.plate.iloc[0][2:]) % 10 for g in hosts])
    cells = [np.clip(np.floor((g.eta.values + HALF) / CELL).astype(int), 0, n - 1) * n
             + np.clip(np.floor((g.xi.values + HALF) / CELL).astype(int), 0, n - 1) for g in hosts]
    sm = lambda m, s: (lambda x: x / x.mean())(gaussian_filter(m.reshape(n, n), s, mode="nearest").ravel())
    score = {}
    for s in SIGMAS:
        tot = W = 0.0
        for k in range(10):
            A = sm(M[fold != k].mean(0), s)
            for h in np.flatnonzero(fold == k):
                w = hosts[h].w.values; tot += w @ np.log(A[cells[h]]); W += w.sum()
        score[s] = tot / W
    best = max(score, key=score.get)
    return sm(M.mean(0), best), best, score


# ------------------------------------------------------------------ the shadow's edge on the sky
# topo: the shadow seen from Palomar, each line of sight followed out to the orbit's radius (the geometry the tests use);
# geo: seen from the centre of the Earth (for the page's comparison of the two).
def shadow_ring(ut, frame, step=1.5, r_orbit=R_GEO):
    """The shadow's edge on the sky at time ut, as (ra, dec) arrays (deg): the circle SHADOW_DEG from the shadow's axis on
    the sphere of radius r_orbit, seen from the centre of the Earth (geo) or from Palomar (topo)."""
    from astropy.time import Time
    t = Time(np.datetime64(np.asarray(ut, dtype="datetime64[ns]")))
    a = np.asarray(antisun_axis(t)).reshape(3)
    e1 = np.cross(a, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(a, e1)
    th = np.radians(np.arange(0, 360 + step / 2, step))[:, None]
    p = r_orbit * (np.cos(np.radians(SHADOW_DEG)) * a + np.sin(np.radians(SHADOW_DEG)) * (np.cos(th) * e1 + np.sin(th) * e2))
    if frame == "topo":
        p = p - PALOMAR.get_gcrs(t).cartesian.xyz.to(u.km).value.reshape(3)
    p = p / np.linalg.norm(p, axis=1, keepdims=True)
    return np.degrees(np.arctan2(p[:, 1], p[:, 0])) % 360, np.degrees(np.arcsin(p[:, 2]))


# ------------------------------------------------------------------ the shadow model
# A glint is a flash at one unknown instant, equally likely at any moment of the exposure. At that instant it is dark if
# the point where Palomar's line of sight meets the orbit is in Earth's shadow; a glint's brightness follows how much of
# the Sun's disc that point sees, so its chance of being recorded is taken to follow it too. A point on the plate is
# therefore in shadow for a share of the exposure: the average over the exposure of the share of the Sun's disc the Earth
# hides from the orbit there. Seen from Palomar the shadow slides across the plate during an exposure (Palomar moves
# about 1,100 km in 50 min); from the centre of the Earth it barely moves. Under the hypothesis that r is the rate of
# transients in shadow over the rate in sunlight, a point in shadow for a share s of the exposure has rate b (1 - (1 - r) s),
# b the rate with no effect; with s 0 or 1 everywhere this is the plain count in shadow against sunlight.
# Edge: "disc" (the hidden share of the Sun's disc, from umbra, 8.42 deg, to penumbra, 8.96 deg, at the geostationary
# radius, half hidden at 8.69); "penumbra" (dark wherever any of it is hidden), with earth > 1 enlarging the Earth for its
# atmosphere (Danjon's 1/85), gives the shadow at its widest, for deciding which plates it can reach.

R_EARTH, R_SUN = 6371.0, 695700.0
DANJON = 1 + 1 / 85
NT = 21                                                                          # instants per exposure
CHUNK = 50000                                                                    # points evaluated at once


@dataclass(frozen=True)
class Shadow:
    frame: str = "topo"
    edge: str = "disc"
    earth: float = 1.0
    r_orbit: float = R_GEO
    nt: int = NT

    def but(self, **k):
        return replace(self, **k)


MAIN = Shadow("topo")


def wide(frame, r_orbit=R_GEO):
    """The shadow at its widest (penumbra's outer edge, the Earth enlarged for its atmosphere), over the whole exposure:
    plates it misses are clear of the shadow."""
    return Shadow(frame, edge="penumbra", earth=DANJON, r_orbit=r_orbit)


def instants(p, S):
    """The times (datetime64) at which plate p's exposure is sampled: S.nt instants spread evenly over it."""
    mid = np.datetime64(np.asarray(p.ut, dtype="datetime64[ns]"))
    k = (np.arange(S.nt) + 0.5) / S.nt - 0.5
    return mid + (k * float(p.exp) * 60e9).astype("timedelta64[ns]")


@lru_cache(maxsize=4096)
def _at(t):
    """Palomar's and the Sun's positions (km, GCRS) and the anti-sun axis at time t (datetime64, as an int of ns)."""
    from astropy.time import Time
    tt = Time(np.datetime64(t, "ns"))
    return (PALOMAR.get_gcrs(tt).cartesian.xyz.to(u.km).value.reshape(3), get_sun(tt).cartesian.xyz.to(u.km).value.reshape(3),
            np.asarray(antisun_axis(tt)).reshape(3))


def orbit_point(d, t, S):
    """Where each line of sight d (unit vectors, n x 3) meets the orbit's sphere at time t (km, GCRS): from the centre of
    the Earth (geo) or from Palomar (topo)."""
    if S.frame == "geo":
        return S.r_orbit * d
    obs = _at(int(np.datetime64(t, "ns").astype(np.int64)))[0]
    a = d @ obs; rng = -a + np.sqrt(a ** 2 - (obs @ obs - S.r_orbit ** 2))
    return obs + rng[:, None] * d


def hidden(P, sun, earth=1.0):
    """Share of the Sun's disc hidden by the Earth, seen from points P (km, n x 3); sun: the Sun's position (km). The Sun
    is small (0.27 deg), so the overlap is the area shared by a small disc and the Earth's limb near it, the limb kept at
    its true angular distance and curvature."""
    dE = np.linalg.norm(P, axis=1); vS = sun[None, :] - P; dS = np.linalg.norm(vS, axis=1)
    rE = np.arcsin(np.clip(earth * R_EARTH / dE, -1, 1)); rS = np.arcsin(R_SUN / dS)
    th = np.arccos(np.clip(np.sum(-P * vS, 1) / (dE * dS), -1, 1))
    b = np.tan(rE); d = b + (th - rE); a = rS                                    # planar: limb's curvature, distance kept
    with np.errstate(invalid="ignore", divide="ignore"):
        lens = (a ** 2 * np.arccos(np.clip((d ** 2 + a ** 2 - b ** 2) / (2 * d * a), -1, 1))
                + b ** 2 * np.arccos(np.clip((d ** 2 + b ** 2 - a ** 2) / (2 * d * b), -1, 1))
                - 0.5 * np.sqrt(np.clip((-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b), 0, None)))
    return np.where(d >= a + b, 0.0, np.where(d <= b - a, 1.0, lens / (np.pi * a ** 2)))


def darkness(d, t, S):
    """How dark the orbit is along lines of sight d (unit vectors) at one time t, under shadow model S (0 sunlit, 1 dark)."""
    sun = _at(int(np.datetime64(t, "ns").astype(np.int64)))[1]
    P = orbit_point(d, t, S)
    h = hidden(P, sun, S.earth)
    return h if S.edge == "disc" else (h > 0).astype(float)


def dark_dirs(d, p, S=MAIN):
    """Share of plate p's exposure in shadow along lines of sight d (unit vectors, n x 3), under model S. p: a plate row
    (ut mid-exposure, exp in min). Taken in chunks of points, the instants summed as they come, so memory stays small."""
    ts = instants(p, S); out = np.empty(len(d))
    for a in range(0, len(d), CHUNK):
        acc = np.zeros(min(CHUNK, len(d) - a))
        for t in ts:
            acc += darkness(d[a:a + CHUNK], t, S)
        out[a:a + CHUNK] = acc / len(ts)
    return out


def dark(xi, eta, p, S=MAIN):
    """dark_dirs at plate-frame offsets (xi, eta) (deg) about p's centre (ra_c, dec_c)."""
    d = unit(*deproject(np.ravel(xi), np.ravel(eta), p.ra_c, p.dec_c)).reshape(-1, 3)
    return dark_dirs(d, p, S).reshape(np.shape(xi))


_SUB = {}
_NONE = np.zeros((GX.size, len(SUB)))
_NONE.flags.writeable = False


def sub_dark(p, S=MAIN):
    """Share in shadow at each cell's 3 x 3 sub-points (cells x 9): the points over which a plate's shadowed area is
    integrated. 0 everywhere, quickly, when the shadow cannot reach the plate during the exposure. Cached."""
    key = (p.name, str(p.ut), S)
    if key not in _SUB:                                                          # plates the shadow misses keep no array
        _SUB[key] = None if near_miss(p, S) else dark((GX[:, None] + SUB[:, 0]).ravel(), (GY[:, None] + SUB[:, 1]).ravel(), p, S).reshape(-1, len(SUB))
    return _NONE if _SUB[key] is None else _SUB[key]


def cell_dark(p, S=MAIN):
    """Each grid cell's share in shadow: its 3 x 3 sub-points averaged."""
    return sub_dark(p, S).mean(1)


def clearance(p, S):
    """How far (deg, about the shadow's axis on the orbit's sphere) plate p's central square stays outside the shadow's
    outer edge (penumbra, with the atmosphere) over its exposure; inf when it is far (near_miss), negative when inside."""
    if near_miss(p, S):
        return np.inf
    d = unit(*deproject(GX, GY, p.ra_c, p.dec_c)).reshape(-1, 3)
    outer = np.degrees(np.arcsin(min(1.0, DANJON * R_EARTH / S.r_orbit))) + np.degrees(np.arcsin(R_SUN / 1.47e8))
    best = np.inf
    for t in instants(p, S):
        P = orbit_point(d, t, S); ax = _at(int(np.datetime64(t, "ns").astype(np.int64)))[2]
        best = min(best, float(np.degrees(np.arccos(np.clip((P @ ax / np.linalg.norm(P, axis=1)).max(), -1, 1)))))
    return best - outer


def near_miss(p, S, margin=2.0):
    """True when plate p's central square cannot reach the shadow at any time of its exposure. Conservative: the square's
    lines of sight lie within HALF * sqrt(2) of its centre's, and meet the orbit within that angle times (r + R_earth) / r
    of where the centre's does, seen from the Earth's centre; the shadow's outer edge (penumbra, with the atmosphere) is
    taken at its widest, and the exposure's start, middle, and end are tested, with a margin for the motion between."""
    c = unit(p.ra_c, p.dec_c).reshape(1, 3)
    outer = np.degrees(np.arcsin(min(1.0, DANJON * R_EARTH / S.r_orbit))) + 0.3
    spread = (HALF * np.sqrt(2) + 0.1) * (S.r_orbit + R_EARTH) / S.r_orbit
    for t in instants(p, S.but(nt=3)):
        P = orbit_point(c, t, S)[0]; ax = _at(int(np.datetime64(t, "ns").astype(np.int64)))[2]
        if np.degrees(np.arccos(np.clip(P @ ax / np.linalg.norm(P), -1, 1))) < outer + spread + margin:
            return False
    return True

