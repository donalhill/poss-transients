"""Every property of a plate that could change how many transients it holds, and the model of a plate's expected count
from them, shared by nuclear.py and nuclear_fit.py.

The properties, chosen in advance and kept whichever way they move a result: distance from the Milky Way (a smooth
spline in the model), dust, exposure time, airmass, the Moon (its brightness while up, and its distance from the
plate), when in the survey the plate was taken (a smooth spline with a knot every half year: of knots every two years,
year, half year and quarter, the one with the lowest BIC; it cannot follow anything as brief as a few days), VASCO's
plate-quality flag, ecliptic latitude, and how deep the plate goes and how sharp its stars are (plate_depth.py), through
which the weather, the emulsion and the processing act. Plates with no exposure in their header, or no depth measured,
take the median. In the model each property other than the splines is per standard deviation."""
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
TERMS = ["lebv", "exp", "airmass", "moon_up", "moon_sep", "quality_bad", "ecl_lat", "depth", "fwhm"]
NAMES = {"b": "distance from the Milky Way", "lebv": "dust", "exp": "exposure time", "airmass": "height above the horizon",
         "moon_up": "moonlight", "moon_sep": "distance from the Moon", "t": "when in the survey it was taken", "quality_bad": "VASCO’s plate-quality flag",
         "ecl_lat": "ecliptic latitude", "depth": "how deep the plate goes", "fwhm": "how sharp its stars are"}
FORMULA = "n ~ bs(b, df=4) + bs(t, df=16) + " + " + ".join(TERMS)


def properties(cat, plates):
    """plates with n (its transients) and every property; a copy."""
    p = plates.copy()
    p["n"] = cat.groupby("plate").size().reindex(p.index).fillna(0)
    sky = SkyCoord(p.ra_c.values * u.deg, p.dec_c.values * u.deg)
    t = Time(p.ut.values)
    pal = EarthLocation.from_geodetic(-116.865 * u.deg, 33.3563 * u.deg, 1706 * u.m)
    moon = get_body("moon", t, pal)
    p["b"] = np.abs(sky.galactic.b.deg)
    p["lebv"] = np.log10(pd.read_csv(ROOT / "data/plate_dust.csv").set_index("plate").ebv.reindex(p.index))
    p["exp"] = pd.read_parquet(ROOT / "data/xe_headers.parquet").drop_duplicates("plate").set_index("plate").exp.reindex(p.index)
    p["exp"] = p.exp.fillna(p.exp.median())
    p["airmass"] = 1 / np.sin(np.radians(np.clip(sky.transform_to(AltAz(obstime=t, location=pal)).alt.deg, 5, 90)))
    p["moon_up"] = (moon.transform_to(AltAz(obstime=t, location=pal)).alt.deg > 0) * (1 - np.cos(np.radians(get_sun(t).separation(moon).deg))) / 2
    p["moon_sep"] = sky.separation(moon).deg
    p["t"] = (p.ut - pd.Timestamp("1949-01-01")).dt.days / 365.25
    p["quality_bad"] = cat.assign(q=cat.plate_quality != "GOOD").groupby("plate").q.mean().reindex(p.index)
    p["ecl_lat"] = np.abs(sky.barycentrictrueecliptic.lat.deg)
    dp = pd.read_csv(ROOT / "data/plate_depth.csv").set_index("plate").reindex(p.index)
    p["depth"] = dp.completeness.fillna(dp.completeness.median())
    p["fwhm"] = dp.fwhm_arcsec.fillna(dp.fwhm_arcsec.median())
    return p


def standardised(p):
    """The model's inputs: n, b, t, and every other property per standard deviation."""
    return p[["n", "b", "t"]].assign(**{k: (p[k] - p[k].mean()) / p[k].std() for k in TERMS})
