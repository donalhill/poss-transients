"""Pan-STARRS DR2 and AllWISE sources around a plate centre, from the catalogues' home archives:
AllWISE from IRSA's TAP service (boxes), Pan-STARRS from MAST's catalogue API (cones of 0.5 deg, its limit; objects detected at least 3 times, since single detections are mostly spurious).
python analysis/catalog_field.py XE320 RA DEC HALFWIDTH_DEG [--wise] [--cone RADIUS_DEG]   -> data/scans/{ps1,wise}_<plate>.parquet"""
import io
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyvo
import requests

ROOT = Path(__file__).resolve().parents[1]
plate, ra0, dec0, hw = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])


def retry(fn, what):
    for attempt in range(6):
        try:
            return fn()
        except Exception as e:
            print(f"{what}: {type(e).__name__}, retrying", flush=True); time.sleep(15)
    raise SystemExit(f"{what} failed")


# ---- AllWISE from IRSA (only with --wise; skipped if already fetched)
irsa = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
w = hw / np.cos(np.radians(abs(dec0) + hw))
q = (f"SELECT ra, dec FROM allwise_p3as_psd WHERE dec BETWEEN {dec0 - hw} AND {dec0 + hw} "
     f"AND ra BETWEEN {ra0 - w} AND {ra0 + w}")
if "--wise" in sys.argv and not (ROOT / f"data/scans/wise_{plate}.parquet").exists():   # not part of the selection; diagnostic only
    wise = retry(lambda: irsa.run_async(q, maxrec=5_000_000).to_table().to_pandas(), "AllWISE")
    wise.to_parquet(ROOT / f"data/scans/wise_{plate}.parquet")
    print(f"wise: {len(wise):,} sources", flush=True)

# ---- Pan-STARRS DR2 (mean object table) from MAST, cones of 0.5 deg on a grid
API = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
CONE = float(sys.argv[sys.argv.index("--cone") + 1]) if "--cone" in sys.argv else 0.5   # radius; smaller in dense fields
step = 1.3 * CONE
centres = []
for d in np.arange(dec0 - hw, dec0 + hw + step, step):
    for r in np.arange(ra0 - w, ra0 + w + step / np.cos(np.radians(d)), step / np.cos(np.radians(d))):
        if d > -30.5:                                     # Pan-STARRS covers dec > -30 only
            centres.append((r % 360, d))


def cone(c):
    cache = ROOT / f"data/scans/ps1_parts/{plate}_{c[0]:.3f}_{c[1]:.3f}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    params = {"ra": c[0], "dec": c[1], "radius": CONE, "columns": "[raMean,decMean]", "nDetections.gte": 3, "pagesize": 500000}
    t = retry(lambda: pd.read_csv(io.StringIO(requests.get(API, params=params, timeout=180).text)), f"PS1 cone {c}")
    t = t.rename(columns={"raMean": "ra", "decMean": "dec"}).dropna()
    if len(t) >= 500000:                                  # the API's page limit: the cone is incomplete
        raise SystemExit(f"PS1 cone {c} hit the 500,000-row limit; use a smaller --cone")
    cache.parent.mkdir(exist_ok=True); t.to_parquet(cache)
    return t


with ThreadPoolExecutor(8) as ex:
    ps1 = pd.concat(list(ex.map(cone, centres))).drop_duplicates()
ps1.to_parquet(ROOT / f"data/scans/ps1_{plate}.parquet")
print(f"ps1: {len(ps1):,} sources from {len(centres)} cones", flush=True)
