"""Gaia DR3 sources around a plate centre, fetched in boxes small enough to stay under the archive's 50,000-row cap.
python analysis/gaia_field.py XE320 RA DEC HALFWIDTH_DEG [N_BOXES]"""
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyvo
from astroquery.gaia import Gaia

ROOT = Path(__file__).resolve().parents[1]
plate, ra0, dec0, hw = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
N = int(sys.argv[5]) if len(sys.argv) > 5 else 6        # N x N boxes in (dec, ra/cos dec); more for dense fields


IRSA = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")


def run(q):
    """IRSA's copy of Gaia DR3 (faster, and up when ESA's archive was down), falling back to ESA's archive."""
    try:
        return IRSA.run_sync(q.replace("gaiadr3.gaia_source", "gaia_dr3_source"), maxrec=200000).to_table().to_pandas()
    except Exception:
        return Gaia.launch_job_async(q).get_results().to_pandas()


def box(ij):
    i, j = ij
    d1, d2 = dec0 - hw + 2 * hw * i / N, dec0 - hw + 2 * hw * (i + 1) / N
    w = hw / np.cos(np.radians(max(abs(d1), abs(d2))))
    r1, r2 = ra0 - w + 2 * w * j / N, ra0 - w + 2 * w * (j + 1) / N
    q = (f"SELECT source_id, ra, dec, phot_g_mean_mag, pmra, pmdec FROM gaiadr3.gaia_source "
         f"WHERE dec BETWEEN {d1} AND {d2} AND ra BETWEEN {r1} AND {r2}")
    cache = ROOT / f"data/scans/gaia_parts/{plate}_{i}_{j}.parquet" if N == 6 else ROOT / f"data/scans/gaia_parts/{plate}_N{N}_{i}_{j}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    for attempt in range(5):                          # the archive drops connections now and then
        try:
            t = run(q); break
        except Exception as e:
            print(f"box {i},{j}: {type(e).__name__}, retrying", flush=True); time.sleep(10)
    else:
        raise SystemExit(f"box {i},{j} failed")
    if len(t) >= 50000:                               # dense field: split the box into four and fetch those
        print(f"box {i},{j} hit the row cap; splitting", flush=True)
        return split(d1, d2, r1, r2, cache.stem)
    cache.parent.mkdir(exist_ok=True); t.to_parquet(cache)
    return t


def split(d1, d2, r1, r2, tag):
    out = []
    for a, (e1, e2) in enumerate([(d1, (d1 + d2) / 2), ((d1 + d2) / 2, d2)]):
        for b, (s1, s2) in enumerate([(r1, (r1 + r2) / 2), ((r1 + r2) / 2, r2)]):
            cache = ROOT / f"data/scans/gaia_parts/{tag}_{a}{b}.parquet"
            if cache.exists():
                out.append(pd.read_parquet(cache)); continue
            q = (f"SELECT source_id, ra, dec, phot_g_mean_mag, pmra, pmdec FROM gaiadr3.gaia_source "
                 f"WHERE dec BETWEEN {e1} AND {e2} AND ra BETWEEN {s1} AND {s2}")
            for attempt in range(5):
                try:
                    t = run(q); break
                except Exception as e:
                    print(f"{tag}_{a}{b}: {type(e).__name__}, retrying", flush=True); time.sleep(10)
            else:
                raise SystemExit(f"{tag}_{a}{b} failed")
            if len(t) >= 50000:
                out.append(split(e1, e2, s1, s2, f"{tag}_{a}{b}")); continue
            t.to_parquet(cache); out.append(t)
    return pd.concat(out)


with ThreadPoolExecutor(6) as ex:
    parts = list(ex.map(box, [(i, j) for i in range(N) for j in range(N)]))
g = pd.concat(parts).drop_duplicates("source_id")
g.to_parquet(ROOT / f"data/scans/gaia_{plate}.parquet")
print(len(g), "Gaia sources; densest box", max(len(p) for p in parts))
