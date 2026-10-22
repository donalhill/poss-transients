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


def cell_counts(xi, eta, w=None):
    """Transients per grid cell from plate-frame offsets (deg); with w, their summed weights."""
    ix = np.floor((np.asarray(xi) + HALF) / CELL).astype(int); iy = np.floor((np.asarray(eta) + HALF) / CELL).astype(int)
    ok = (ix >= 0) & (ix < len(g)) & (iy >= 0) & (iy < len(g))
    return np.bincount(iy[ok] * len(g) + ix[ok], weights=None if w is None else np.asarray(w)[ok], minlength=GX.size)


def shadow_share(sub_xyz, ut):
    """Share of each cell inside Earth's shadow at geostationary height: the geocentric anti-sun direction at
    mid-exposure, as in VASCO's analysis."""
    from astropy.time import Time
    import pandas as pd
    ax = np.asarray(antisun_axis(Time(pd.Timestamp(ut).to_datetime64()))).reshape(3)
    ang = np.degrees(np.arccos(np.clip(sub_xyz @ ax, -1, 1)))
    return (ang < SHADOW_DEG).reshape(-1, len(SUB)).mean(1)


def vasco_plate(plate, p, cat):
    """One plate of VASCO's catalogue on the grid: transient counts per cell and the sky position of each sub-point.
    p: the plate's row of plates_corrected.parquet (centre, corrected mid-exposure time)."""
    sra, sdec = deproject((GX[:, None] + SUB[:, 0]).ravel(), (GY[:, None] + SUB[:, 1]).ravel(), p.ra_c, p.dec_c)
    c = cat[cat.plate == plate]
    return dict(plate=plate, n=cell_counts(c.xi.values, c.eta.values), use=np.ones(GX.size, bool),
                sub=unit(sra, sdec), ut=p.ut)
