import math
import warnings
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, get_sun
from erfa import ErfaWarning
from astropy.utils.exceptions import AstropyWarning

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
