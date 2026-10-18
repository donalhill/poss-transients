"""Does the published 0.692% expected shadow fraction match the method used for the transients?

Applies several shadow definitions to both the plate-aware control sample (635 plates x 100)
and the 107,875 transient catalogue.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, GCRS, get_sun
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False
iers.conf.iers_degraded_accuracy = "warn"

sys.path.insert(0, str(Path(__file__).parent))
from shadow import ROOT, PALOMAR, R_GEO, SHADOW_DEG, antisun_geocentric, sep_topocentric

ctrl = pd.read_csv(ROOT / "VASCO-challenge/data/XE_final_random_sources.csv")
cat = pd.read_csv(ROOT / "VASCO-ML/data/SUPERVIKTIG_HELAVASCO_validated_v4.csv")
strict = set(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv")["obs_date"].astype(str).str[:10])

ctrl["ut"] = pd.to_datetime(ctrl["UTobservation"])
ctrl["date"] = ctrl["ut"].dt.strftime("%Y-%m-%d")
print(f"control: {len(ctrl)} positions, {ctrl.Name.nunique()} plates, UT hours {ctrl.ut.dt.hour.min()}-{ctrl.ut.dt.hour.max()}")

def frac(mask):
    return f"{mask.mean()*100:.3f}% ({mask.sum()})"

# 1. Published transient method: geocentric anti-sun, 08:00 UT on obs_date, no precession
for name, df, dcol in [("control", ctrl, "date"), ("catalogue", cat, "obs_date")]:
    t08 = Time(pd.to_datetime(df[dcol].astype(str).str[:10]) + pd.Timedelta(hours=8))
    for prec in (False, True):
        sep = antisun_geocentric(df["RA" if name == "control" else "ra"].values,
                                 df["Dec" if name == "control" else "dec"].values, t08, precess=prec)
        m = sep < SHADOW_DEG
        print(f"{name:9s} geocentric 08UT {'J2000' if prec else 'of-date'}: {frac(m)}")
        if name == "catalogue" and not prec:
            ss = cat["obs_date"].astype(str).str[:10].isin(strict)
            print(f"           strict dates only (N={ss.sum()}): {frac(m[ss.values])}")

# 2. Control sample with actual plate UT and proper topocentric GEO shadow
t = Time(ctrl["ut"].values)
sep_geo = antisun_geocentric(ctrl.RA.values, ctrl.Dec.values, t, precess=True)
sep_topo = sep_topocentric(ctrl.RA.values, ctrl.Dec.values, t)
print(f"control geocentric, plate UT, J2000: {frac(sep_geo < SHADOW_DEG)}")
print(f"control topocentric GEO, plate UT:   {frac(sep_topo < SHADOW_DEG)}")
