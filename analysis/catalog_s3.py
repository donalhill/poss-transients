"""Gaia DR3 and Pan-STARRS DR2 around a plate centre, read straight from STScI's public copies on AWS (HATS Parquet,
s3://stpubdata, no account needed): only the sky partitions that overlap the plate, and only the columns we use.
Same outputs as gaia_field.py and catalog_field.py (Pan-STARRS: objects detected at least 3 times).
python analysis/catalog_s3.py XE320 RA DEC HALFWIDTH_DEG   -> data/scans/{gaia,ps1}_<plate>.parquet"""
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from pyarrow import fs

ROOT = Path(__file__).resolve().parents[1]
plate, ra0, dec0, hw = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
S3 = fs.S3FileSystem(anonymous=True, region="us-east-1")
CATS = {"gaia": ("stpubdata/gaia/gaia_dr3/public/hats/gaia", ["source_id", "ra", "dec", "phot_g_mean_mag", "pmra", "pmdec"],
                 "ra", "dec", None),
        "ps1": ("stpubdata/panstarrs/ps1/public/hats/otmo/otmo", ["raMean", "decMean", "nDetections"],
                "raMean", "decMean", pc.field("nDetections") >= 3)}
w = hw / np.cos(np.radians(min(89.0, abs(dec0) + hw)))
d1, d2 = dec0 - hw, dec0 + hw
corners = [(ra0 - w, d1), (ra0 + w, d1), (ra0 + w, d2), (ra0 - w, d2)]
vec = np.array([hp.ang2vec(r % 360, d, lonlat=True) for r, d in corners])


def partitions(root):
    with S3.open_input_file(f"{root}/partition_info.csv") as f:
        info = pd.read_csv(f)
    keep = []
    for order, grp in info.groupby("Norder"):
        hit = set(hp.query_polygon(2 ** order, vec, inclusive=True, nest=True).tolist())
        keep += [(order, p) for p in grp.Npix if p in hit]
    return [f"{root}/dataset/Norder={o}/Dir={p // 10000 * 10000}/Npix={p}.parquet" for o, p in keep]


def read(path, cols, rac, decc, extra):
    for attempt in range(5):
        try:
            flt = (pc.field(decc) >= d1) & (pc.field(decc) <= d2)
            flt = flt & extra if extra is not None else flt
            t = ds.dataset(path, filesystem=S3, format="parquet").to_table(columns=cols, filter=flt).to_pandas()
            dra = (t[rac] - ra0 + 180) % 360 - 180
            return t[np.abs(dra) <= w]
        except Exception as e:
            print(f"{path.split('/')[-1]}: {type(e).__name__}, retrying", flush=True); time.sleep(10)
    raise SystemExit(f"{path} failed")


for key, (root, cols, rac, decc, extra) in CATS.items():
    out = ROOT / f"data/scans/{key}_{plate}.parquet"
    if out.exists():
        continue
    t0 = time.time()
    paths = partitions(root)
    with ThreadPoolExecutor(8) as ex:
        t = pd.concat(list(ex.map(lambda p: read(p, cols, rac, decc, extra), paths)), ignore_index=True)
    if key == "ps1":
        t = t.rename(columns={"raMean": "ra", "decMean": "dec"})[["ra", "dec"]]
    t.to_parquet(out)
    print(f"{key}: {len(t):,} sources from {len(paths)} partitions in {time.time() - t0:.0f} s", flush=True)
