"""Every property of a plate that could change how many transients it holds, used by nuclear.py: distance from the
Milky Way, dust, exposure time, airmass, the Moon (its brightness while up, and its distance from the plate), VASCO's
plate-quality flag, ecliptic latitude, and how deep the plate goes and how sharp its stars are (plate_depth.py). Plates
with no exposure in their header, or no depth measured, take the median."""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter("ignore", AstropyWarning)

ROOT = Path(__file__).resolve().parents[1]
TERMS = ["lebv", "lebv_plate", "exp", "airmass", "moon_up", "moon_sep", "quality_bad", "ecl_lat", "depth", "fwhm"]


def properties(cat, plates):
    """plates with n (its transients) and every property; a copy."""
    p = plates.copy()
    p["n"] = cat.groupby("plate").size().reindex(p.index).fillna(0)
    sky = SkyCoord(p.ra_c.values * u.deg, p.dec_c.values * u.deg)
    t = Time(p.ut.values)
    pal = EarthLocation.from_geodetic(-116.865 * u.deg, 33.3563 * u.deg, 1706 * u.m)
    moon = get_body("moon", t, pal)
    p["b"] = np.abs(sky.galactic.b.deg)
    dust = pd.read_csv(ROOT / "data/plate_dust.csv").set_index("plate").reindex(p.index)
    p["lebv"] = np.log10(dust.ebv)                                               # at the plate centre
    p["lebv_plate"] = np.log10(dust.ebv_plate)                                   # averaged over the plate's searched area
    p["exp"] = pd.read_parquet(ROOT / "data/xe_headers.parquet").drop_duplicates("plate").set_index("plate").exp.reindex(p.index)
    p["exp"] = p.exp.fillna(p.exp.median())
    p["airmass"] = 1 / np.sin(np.radians(np.clip(sky.transform_to(AltAz(obstime=t, location=pal)).alt.deg, 5, 90)))
    p["moon_up"] = (moon.transform_to(AltAz(obstime=t, location=pal)).alt.deg > 0) * (1 - np.cos(np.radians(get_sun(t).separation(moon).deg))) / 2
    p["moon_sep"] = sky.separation(moon).deg
    p["quality_bad"] = cat.assign(q=cat.plate_quality != "GOOD").groupby("plate").q.mean().reindex(p.index)
    p["ecl_lat"] = np.abs(sky.barycentrictrueecliptic.lat.deg)
    # the row of plates centred on the celestial equator: on the same sky, these plates hold about two-fifths as many
    # transients as their neighbours to the north, and over the whole plate about a quarter of what their pointing
    # otherwise predicts. The deficit is already in Solano et al.'s published list, so it comes from the detection pipeline
    p["eq_row"] = (np.abs(p.dec_c) < 3).astype(float)
    dp = pd.read_csv(ROOT / "data/plate_depth.csv").set_index("plate").reindex(p.index)
    p["depth"] = dp.completeness.fillna(dp.completeness.median())
    p["fwhm"] = dp.fwhm_arcsec.fillna(dp.fwhm_arcsec.median())
    return p


def standardised(p):
    """The model's inputs: n, b, and every other property per standard deviation."""
    return p[["n", "b", "eq_row"]].assign(**{k: (p[k] - p[k].mean()) / p[k].std() for k in TERMS})


POINTING = "bs(b, df=4) + lebv_plate + ecl_lat + eq_row"


def pointing(p, use=None, formula=POINTING):
    """The pointing curve: the transients each plate's pointing predicts, from a negative binomial fit of distance from the
    Milky Way (a smooth curve), dust averaged over the plate, distance from the ecliptic, and whether the plate is in the
    equatorial row. Fitted
    to the plates in use (every plate by default) and predicting every plate. Returns the expected counts, the design,
    and the dispersion."""
    import patsy
    import statsmodels.api as sm
    y = p.n.values.astype(float); u = np.ones(len(y), bool) if use is None else np.asarray(use)
    X = np.asarray(patsy.dmatrix(formula, standardised(p))); a = 0.5
    for _ in range(4):
        E = sm.GLM(y[u], X[u], family=sm.families.NegativeBinomial(alpha=a)).fit().predict(X)
        a = float(np.sum(((y[u] - E[u]) ** 2 - E[u]) / E[u] ** 2) / (u.sum() - X.shape[1]))
    return E, X, a
