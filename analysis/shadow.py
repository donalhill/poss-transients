import warnings
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, get_sun
from erfa import ErfaWarning
from astropy.utils.exceptions import AstropyWarning
from astropy.utils import iers

from geometry import deproject, unit                                            # noqa: F401  (used here and re-exported)

iers.conf.auto_download = False     # the Sun's position needs no fresh Earth-rotation tables; downloading them stalls runs

warnings.simplefilter("ignore", ErfaWarning)
warnings.simplefilter("ignore", AstropyWarning)

ROOT = Path(__file__).resolve().parents[1]
PALOMAR = EarthLocation.from_geodetic(lon=-116.8650 * u.deg, lat=33.3563 * u.deg, height=1706 * u.m)
R_GEO = 42164.0
# geocentric half-angle of Earth's shadow at geostationary radius, asin(6371 / 42164): the 8.69 deg of Villarroel et al.
# (2025, PASP, Table 4), whose shadow result the page tests. VASCO-ML's code for Bruehl et al. (2026) uses 8.50, between
# this and the umbra (8.42); shadow_pooled.py checks the page's test at that radius too
SHADOW_DEG = 8.69
SHADOW_DEG_ML = 8.50


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


def antisun_geocentric(ra, dec, t):
    """Angle (deg) between the line of sight and the geocentric anti-sun direction (no parallax)."""
    cosang = np.sum(unit(ra, dec) * antisun_axis(t), axis=-1)
    return np.degrees(np.arccos(np.clip(cosang, -1, 1)))


def sep_topocentric(ra, dec, t, r_orbit=R_GEO, return_alt=False):
    """Geocentric angle from the shadow axis of the point where the line of sight from Palomar
    meets the sphere of radius r_orbit (km). Comparable to SHADOW_DEG."""
    d = unit(ra, dec)
    tu, inv = distinct(t)
    obs = PALOMAR.get_gcrs(tu).cartesian.xyz.to(u.km).value.T  # (N, 3)
    obs = obs if inv is None else obs[inv]
    a = np.sum(obs * d, axis=-1)
    b = np.sum(obs**2, axis=-1) - r_orbit**2
    rng = -a + np.sqrt(a**2 - b)  # forward root; observer is inside the sphere
    p = obs + rng[:, None] * d
    cosang = np.sum(p * antisun_axis(t), axis=-1) / r_orbit
    sep = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
    if return_alt:
        return sep, np.degrees(np.arcsin(a / np.linalg.norm(obs, axis=-1)))
    return sep


# ------------------------------------------------------------------ plate grid for the within-plate tests
# Each plate is gridded into 0.1 deg cells over the central 6 x 6 deg, the plate size VASCO use (Villarroel et al. 2025);
# VASCO's random positions cover this square evenly on every plate. The shadow share of a cell is taken from 3 x 3
# points inside it.
HALF, CELL = 3.0, 0.1
g = np.arange(-HALF + CELL / 2, HALF, CELL)
GX, GY = [a.ravel() for a in np.meshgrid(g, g)]
SUB = np.array([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]) * CELL / 3


def cell_index(xi, eta):
    """The grid cell (index into GX, GY) of each plate-frame offset (deg), clipped to the grid."""
    ix = np.clip(np.floor((np.asarray(xi) + HALF) / CELL).astype(int), 0, len(g) - 1)
    iy = np.clip(np.floor((np.asarray(eta) + HALF) / CELL).astype(int), 0, len(g) - 1)
    return iy * len(g) + ix


def cell_counts(xi, eta, w=None):
    """Transients per grid cell from plate-frame offsets (deg); with w, their summed weights."""
    ix = np.floor((np.asarray(xi) + HALF) / CELL).astype(int); iy = np.floor((np.asarray(eta) + HALF) / CELL).astype(int)
    ok = (ix >= 0) & (ix < len(g)) & (iy >= 0) & (iy < len(g))
    return np.bincount(iy[ok] * len(g) + ix[ok], weights=None if w is None else np.asarray(w)[ok], minlength=GX.size)


# ------------------------------------------------------------------ the two frames for Earth's shadow
# geo: the shadow as VASCO place it, about the anti-solar direction seen from the centre of the Earth. topo: the shadow
# seen from Palomar, each line of sight followed out to the orbit's radius (the more exact geometry for an object there).
# Every shadow test runs in both, through these functions alone.
FRAMES = ("geo", "topo")


def shadow_sep(ra, dec, ut, frame, r_orbit=R_GEO):
    """Angle (deg) from the shadow's axis, comparable to the shadow's radius, of each line of sight (ra, dec) at time
    ut (one time, or one per line of sight), in the given frame."""
    from astropy.time import Time
    ra, dec = np.atleast_1d(ra), np.atleast_1d(dec)
    t = Time(np.broadcast_to(np.asarray(ut, dtype="datetime64[ns]"), ra.shape))
    return antisun_geocentric(ra, dec, t) if frame == "geo" else sep_topocentric(ra, dec, t, r_orbit)


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


def shadow_share(sub_xyz, ut, frame="geo"):
    """Share of each grid cell inside Earth's shadow at geostationary height, from its 3 x 3 sub-points (unit vectors)."""
    ra = np.degrees(np.arctan2(sub_xyz[:, 1], sub_xyz[:, 0])) % 360
    dec = np.degrees(np.arcsin(np.clip(sub_xyz[:, 2], -1, 1)))
    return (shadow_sep(ra, dec, ut, frame) < SHADOW_DEG).reshape(-1, len(SUB)).mean(1)


# ------------------------------------------------------------------ the shadow model
# A glint is a flash at one unknown instant, equally likely at any moment of the exposure. At that instant it is dark if
# the point where Palomar's line of sight meets the orbit is in Earth's shadow; a glint's brightness follows how much of
# the Sun's disc that point sees, so its chance of being recorded is taken to follow it too. A point on the plate is
# therefore in shadow for a share of the exposure: the average over the exposure of the share of the Sun's disc the Earth
# hides from the orbit there. Seen from Palomar the shadow slides across the plate during an exposure (Palomar moves
# about 1,100 km in 50 min); from the centre of the Earth it barely moves. Under the hypothesis that r is the rate of
# transients in shadow over the rate in sunlight, a point in shadow for a share s of the exposure has rate b (1 - (1 - r) s),
# b the rate with no effect; with s 0 or 1 everywhere this is the plain count in shadow against sunlight.
# Edge models: "disc" (default: the hidden share of the Sun's disc, from umbra, 8.42 deg, to penumbra, 8.96 deg, at the
# geostationary radius, half hidden at 8.69); "umbra" (dark only where the whole disc is hidden) and "penumbra" (dark
# wherever any of it is) bracket how glints' detection follows the light; "sharp" (a cone of the given radius about the
# shadow's axis, VASCO's construction). earth > 1 enlarges the Earth for its atmosphere (Danjon's 1/85).
from dataclasses import dataclass, replace                                       # noqa: E402

R_EARTH, R_SUN = 6371.0, 695700.0
DANJON = 1 + 1 / 85
NT = 21                                                                          # instants per exposure
CHUNK = 50000                                                                    # points evaluated at once


@dataclass(frozen=True)
class Shadow:
    frame: str = "topo"
    edge: str = "disc"
    radius: float = SHADOW_DEG                                                   # sharp edge only
    earth: float = 1.0
    instant: bool = False                                                        # one instant (mid-exposure, or as given)
    shift_min: float = 0.0                                                       # every plate's times moved by this
    r_orbit: float = R_GEO
    nt: int = NT

    def but(self, **k):
        return replace(self, **k)


MAIN = Shadow("topo")
GEO = Shadow("geo")


def wide(frame, r_orbit=R_GEO):
    """The shadow at its widest (penumbra's outer edge, the Earth enlarged for its atmosphere), over the whole exposure: no
    variant reaches beyond it, so plates it misses are clear of every one."""
    return Shadow(frame, edge="penumbra", earth=DANJON, r_orbit=r_orbit)


def variants(S):
    """The model's systematics: how glints' detection follows the light (umbra only, any dimming), the atmosphere, the
    plates' times (the asteroid trails put mid-exposure 0.5 min before ours, scatter 0.7 min), the exposure taken as one
    instant, a sharp edge at VASCO's 8.69 deg over the exposure, and the sampling of the exposure."""
    return {"umbra": S.but(edge="umbra"), "penumbra": S.but(edge="penumbra"), "atmosphere": S.but(earth=DANJON),
            "times_minus_1": S.but(shift_min=-1.0), "times_plus_1": S.but(shift_min=1.0),
            "times_minus_5": S.but(shift_min=-5.0), "times_plus_5": S.but(shift_min=5.0),
            "mid_exposure": S.but(instant=True), "sharp": S.but(edge="sharp"), "instants_41": S.but(nt=41)}


def instants(p, S):
    """The times (datetime64) at which plate p's exposure is sampled: NT instants spread evenly over it, or one."""
    mid = np.datetime64(np.asarray(p.ut, dtype="datetime64[ns]")) + np.timedelta64(int(round(S.shift_min * 60e9)), "ns")
    if S.instant:
        return np.array([mid])
    k = (np.arange(S.nt) + 0.5) / S.nt - 0.5
    return mid + (k * float(p.exp) * 60e9).astype("timedelta64[ns]")


from functools import lru_cache                                                  # noqa: E402


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
    _, sun, ax = _at(int(np.datetime64(t, "ns").astype(np.int64)))
    P = orbit_point(d, t, S)
    if S.edge == "sharp":
        return (P @ ax / np.linalg.norm(P, axis=1) > np.cos(np.radians(S.radius))).astype(float)
    h = hidden(P, sun, S.earth)
    return h if S.edge == "disc" else (h >= 1 - 1e-12).astype(float) if S.edge == "umbra" else (h > 0).astype(float)


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


def turn(xi, eta, d):
    """Plate-frame offsets turned by d quarter-turns about the plate centre (d % 4), and mirrored (d >= 4): the eight
    ways an outline can be laid on another plate."""
    xi, eta = np.asarray(xi), np.asarray(eta)
    for _ in range(d % 4):
        xi, eta = -eta, xi
    return (eta, xi) if d >= 4 else (xi, eta)


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


def cell_dark(p, S=MAIN, centres=False):
    """Each grid cell's share in shadow: its 3 x 3 sub-points averaged, or (centres, VASCO's construction) its centre
    alone."""
    if centres:
        return np.zeros(GX.size) if near_miss(p, S) else dark(GX, GY, p, S)
    return sub_dark(p, S).mean(1)


def clearance(p, S):
    """How far (deg, about the shadow's axis on the orbit's sphere) plate p's central square stays outside the shadow's
    outer edge (penumbra, with the atmosphere) over its exposure; inf when it is far (near_miss), negative when inside."""
    if near_miss(p, S):
        return np.inf
    d = unit(*deproject(GX, GY, p.ra_c, p.dec_c)).reshape(-1, 3)
    outer = max(np.degrees(np.arcsin(min(1.0, DANJON * R_EARTH / S.r_orbit))) + np.degrees(np.arcsin(R_SUN / 1.47e8)),
                S.radius if S.edge == "sharp" else 0.0)
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
    outer = max(np.degrees(np.arcsin(min(1.0, DANJON * R_EARTH / S.r_orbit))) + 0.3, S.radius if S.edge == "sharp" else 0.0)
    spread = (HALF * np.sqrt(2) + 0.1) * (S.r_orbit + R_EARTH) / S.r_orbit
    for t in instants(p, S if S.instant else S.but(nt=3)):
        P = orbit_point(c, t, S)[0]; ax = _at(int(np.datetime64(t, "ns").astype(np.int64)))[2]
        if np.degrees(np.arccos(np.clip(P @ ax / np.linalg.norm(P), -1, 1))) < outer + spread + margin:
            return False
    return True


def vasco_plate(plate, p, cat):
    """One plate of VASCO's catalogue on the grid: transient counts per cell and the sky position of each sub-point.
    p: the plate's row of plates_corrected.parquet (centre, corrected mid-exposure time)."""
    sra, sdec = deproject((GX[:, None] + SUB[:, 0]).ravel(), (GY[:, None] + SUB[:, 1]).ravel(), p.ra_c, p.dec_c)
    c = cat[cat.plate == plate]
    return dict(plate=plate, n=cell_counts(c.xi.values, c.eta.values), use=np.ones(GX.size, bool),
                sub=unit(sra, sdec), ut=p.ut)
