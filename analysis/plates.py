"""Assign each catalogue transient to a POSS-I XE plate (same night, nearest centre, inside footprint).

Writes data/transients_with_plates.parquet and data/plates.parquet.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from geometry import tangent, unit

ROOT = Path(__file__).resolve().parents[1]

OUT = ROOT / "data"
OUT.mkdir(exist_ok=True)

ctrl = pd.read_csv(ROOT / "data/external/XE_final_random_sources.csv")
cat = pd.read_csv(ROOT / "data/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv",
                  usecols=["ra", "dec", "obs_date", "plate_quality", "snr", "near_defect"])

v = unit(ctrl.RA.values, ctrl.Dec.values)
ctrl[["x", "y", "z"]] = v
plates = ctrl.groupby("Name").agg(x=("x", "mean"), y=("y", "mean"), z=("z", "mean"), ut=("UTobservation", "first"))
norm = np.sqrt(plates.x**2 + plates.y**2 + plates.z**2)
plates["ra_c"] = np.degrees(np.arctan2(plates.y, plates.x)) % 360
plates["dec_c"] = np.degrees(np.arcsin(plates.z / norm))
plates["ut"] = pd.to_datetime(plates.ut)
plates["ut_date"] = plates.ut.dt.strftime("%Y-%m-%d")

# footprint half-width from the control scatter in a tangent plane, per plate
c = ctrl.join(plates[["ra_c", "dec_c"]], on="Name")
xi, eta = tangent(c.RA, c.Dec, c.ra_c, c.dec_c)
print(f"control footprint |xi| max {np.abs(xi).max():.2f}, |eta| max {np.abs(eta).max():.2f} deg "
      f"(99th pct {np.percentile(np.abs(xi), 99):.2f}, {np.percentile(np.abs(eta), 99):.2f})")
HALF = float(np.ceil(max(np.percentile(np.abs(xi), 99.9), np.percentile(np.abs(eta), 99.9)) * 10) / 10)

def assign(date_shift):
    cat_dates = (pd.to_datetime(cat.obs_date) + pd.Timedelta(days=date_shift)).dt.strftime("%Y-%m-%d")
    best = np.full(len(cat), "", dtype=object)
    bxi = np.full(len(cat), np.nan)
    beta = np.full(len(cat), np.nan)
    for d, idx in cat.groupby(cat_dates).groups.items():
        pl = plates[plates.ut_date == d]
        if pl.empty:
            continue
        sub = cat.loc[idx]
        dmin = np.full(len(sub), np.inf)
        for name, p in pl.iterrows():
            x, e = tangent(sub.ra.values, sub.dec.values, p.ra_c, p.dec_c)
            r = np.maximum(np.abs(x), np.abs(e))
            better = (r < dmin) & (r <= HALF + 0.2)
            dmin = np.where(better, r, dmin)
            best[np.asarray(idx)[better]] = name
            bxi[np.asarray(idx)[better]] = x[better]
            beta[np.asarray(idx)[better]] = e[better]
    return best, bxi, beta

for shift in (0, -1, 1):
    b, _, _ = assign(shift)
    print(f"obs_date {shift:+d} day -> assigned {np.mean(b != '')*100:.1f}%")

best, bxi, beta = assign(0)
cat["plate"] = best
cat["xi"], cat["eta"] = bxi, beta
cat.to_parquet(OUT / "transients_with_plates.parquet")
plates["half"] = HALF
plates.to_parquet(OUT / "plates.parquet")
un = cat[cat.plate == ""]
print(f"half-width used {HALF} deg; unassigned {len(un)} ({len(un)/len(cat)*100:.1f}%) on {un.obs_date.nunique()} dates")
print(f"plates with transients: {cat[cat.plate != ''].plate.nunique()} of {len(plates)}")
n = cat[cat.plate != ""].groupby("plate").size()
print(f"transients per plate: mean {n.mean():.0f}, median {n.median():.0f}, max {n.max()}, var/mean {n.var()/n.mean():.0f}")
