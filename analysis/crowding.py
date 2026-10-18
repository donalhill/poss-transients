"""How crowded is a plate for the "no Gaia or Pan-STARRS source within 5 arcsec" requirement? Share of random positions
within 0.25 deg of the plate centre that have a catalogued source within 5". Uses whichever of Gaia and Pan-STARRS
have been fetched; Pan-STARRS alone already decides the dense fields.
python analysis/crowding.py XE701 RA DEC   -> prints the share, and appends it to data/crowding.csv"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
plate, ra0, dec0 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
rng = np.random.default_rng(7)


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


n = 3000
pr = rng.uniform(-0.25, 0.25, n) / np.cos(np.radians(dec0)) + ra0
pd_ = rng.uniform(-0.25, 0.25, n) + dec0
hit = np.zeros(n, bool)
used = []
for key in ("ps1", "gaia"):
    f = ROOT / f"data/scans/{key}_{plate}.parquet"
    if not f.exists():
        continue
    c = pd.read_parquet(f, columns=["ra", "dec"])
    c = c[(np.abs(c.dec - dec0) < 0.3) & (np.abs((c.ra - ra0 + 180) % 360 - 180) * np.cos(np.radians(dec0)) < 0.3)]
    d, _ = cKDTree(xyz(c.ra.values, c.dec.values)).query(xyz(pr, pd_))
    hit |= np.degrees(2 * np.arcsin(d / 2)) * 3600 < 5
    used.append(key)
out = ROOT / "data/crowding.csv"
with open(out, "a") as fh:
    if out.stat().st_size == 0:
        fh.write("plate,catalogues,veto\n")
    fh.write(f"{plate},{'+'.join(used)},{hit.mean():.3f}\n")
print(f"{plate}: {hit.mean():.1%} of random positions have a {' or '.join(used)} source within 5\"")
