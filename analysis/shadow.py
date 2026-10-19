import math
import warnings
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, get_sun
from erfa import ErfaWarning
from astropy.utils.exceptions import AstropyWarning
from astropy.utils import iers

iers.conf.auto_download = False     # the Sun's position needs no fresh Earth-rotation tables; downloading them stalls runs

warnings.simplefilter("ignore", ErfaWarning)
warnings.simplefilter("ignore", AstropyWarning)

ROOT = Path(__file__).resolve().parents[1]
PALOMAR = EarthLocation.from_geodetic(lon=-116.8650 * u.deg, lat=33.3563 * u.deg, height=1706 * u.m)
R_EARTH = 6371.0
R_GEO = 42164.0
SHADOW_DEG = 8.50  # geocentric half-angle used by VASCO / Doherty


def unit(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def sun_of_date(jd):
    """Doherty's low-precision Sun (RA/Dec referred to the equinox of date)."""
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360
    g = np.radians((357.528 + 0.9856003 * n) % 360)
    lam = np.radians(L + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    eps = np.radians(23.439 - 0.0000004 * n)
    ra = np.degrees(np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))) % 360
    dec = np.degrees(np.arcsin(np.sin(eps) * np.sin(lam)))
    return ra, dec


def antisun_axis(t, precess=True):
    """Unit vector of the anti-solar direction. precess=True: GCRS (ICRS-aligned, J2000 axes)."""
    if precess:
        s = get_sun(t)
        return -unit(s.ra.deg, s.dec.deg)
    ra, dec = sun_of_date(t.jd)
    return -unit(ra, dec)


def antisun_geocentric(ra, dec, t, precess=True):
    """Angle (deg) between the line of sight and the geocentric anti-sun direction (no parallax)."""
    cosang = np.sum(unit(ra, dec) * antisun_axis(t, precess), axis=-1)
    return np.degrees(np.arccos(np.clip(cosang, -1, 1)))


def sep_topocentric(ra, dec, t, r_orbit=R_GEO, return_alt=False):
    """Geocentric angle from the shadow axis of the point where the line of sight from Palomar
    meets the sphere of radius r_orbit (km). Comparable to SHADOW_DEG."""
    d = unit(ra, dec)
    obs = PALOMAR.get_gcrs(t).cartesian.xyz.to(u.km).value.T  # (N, 3)
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


def deproject(xi, eta, ra0, dec0):
    """Plate-frame offsets (deg) -> RA, Dec (deg), gnomonic."""
    xi, eta, ra0, dec0 = map(np.radians, (xi, eta, ra0, dec0))
    rho = np.hypot(xi, eta); c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def cell_counts(xi, eta):
    """Transients per grid cell from plate-frame offsets (deg)."""
    ix = np.floor((np.asarray(xi) + HALF) / CELL).astype(int); iy = np.floor((np.asarray(eta) + HALF) / CELL).astype(int)
    ok = (ix >= 0) & (ix < len(g)) & (iy >= 0) & (iy < len(g))
    return np.bincount(iy[ok] * len(g) + ix[ok], minlength=GX.size)


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
