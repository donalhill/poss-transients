"""Which POSS-I red plates, from Palomar's own log, were crossed by Earth's shadow at geostationary height?

For every plate in the log (VizieR VI/25, data/palomar_log_vi25.csv): mid-exposure time in UT (local start time is
Pacific Standard Time, exposure from the log), the share of the plate's central 6 x 6 deg in shadow over its exposure
under the page's shadow model (shadow.MAIN, seen from Palomar), and whether the plate is in VASCO's catalogue of 107,875
transients. Plates crossed by the
shadow but missing from the catalogue are the ones that would enlarge the shadow test.

python analysis/shadow_plates_log.py   -> data/shadow_plates_log.csv"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import shadow as SH   # noqa: E402

log = pd.read_csv(ROOT / "data/palomar_log_vi25.csv").dropna(subset=["red_start_local", "ra_j2000", "dec_j2000"])
c = SkyCoord(log.ra_j2000.values, log.dec_j2000.values, unit=(u.hourangle, u.deg))
log["ra"], log["dec"] = c.ra.deg, c.dec.deg
hh = log.red_start_local.str.split(":").apply(lambda v: int(v[0]) + int(v[1]) / 60)
day = pd.to_datetime(log.evening_date) + pd.to_timedelta((hh < 12).astype(int), unit="D")   # after midnight
log["ut"] = (day + pd.to_timedelta(hh, unit="h") + pd.Timedelta(hours=8)
             + pd.to_timedelta(log.red_exposure_min.fillna(45) / 2, unit="m"))


log["exp"] = log.red_exposure_min.fillna(45.0)
frac = [float(SH.sub_dark(pd.Series({"ra_c": p.ra, "dec_c": p.dec, "ut": p.ut, "exp": p.exp}, name=i), SH.MAIN).mean())
        for i, p in log.iterrows()]
log["shadow_fraction"] = frac

pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
inc = pl[pl.index.isin(set(cat.plate.unique()) - {""})]


def in_catalogue(p):
    d = np.hypot((inc.ra_c - p.ra) * np.cos(np.radians(p.dec)), inc.dec_c - p.dec)
    return bool(((d < 1.0) & ((inc.ut - p.ut).abs() < pd.Timedelta(days=1.5))).any())


sh = log[log.shadow_fraction > 0].copy()
sh["in_catalogue"] = [in_catalogue(p) for _, p in sh.iterrows()]
cols = ["plate_number", "evening_date", "ut", "ra", "dec", "shadow_fraction", "in_catalogue"]
sh[cols].sort_values(["in_catalogue", "shadow_fraction"], ascending=[True, False]).round(3).to_csv(ROOT / "data/shadow_plates_log.csv", index=False)
for flag, name in ((True, "in VASCO's catalogue"), (False, "not in the catalogue")):
    s = sh[sh.in_catalogue == flag]
    print(f"{name}: {len(s)} plates crossed by the shadow, {s.shadow_fraction.sum():.1f} plates' worth of shadowed sky, "
          f"declinations {s.dec.min():+.0f} to {s.dec.max():+.0f}")
