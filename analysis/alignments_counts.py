"""How many r-point alignments occur in the 107k catalogue, and how many does a clustering-preserving null produce?
r-point alignment: a pair (i, j) with |ij| <= LMAX and >= r-2 other transients between them within a strip of full width W.
Null: each transient displaced uniformly within a 60" disk (keeps structure > 1', destroys arcsec collinearity)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
LMAX = 900.0
rng = np.random.default_rng(11)
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
cat = cat[cat.plate != ""]
cat["x"], cat["y"] = cat.xi * 3600, cat.eta * 3600


def count(df, widths=(2.0, 5.0, 10.0, 15.0)):
    out = {(w, r): 0 for w in widths for r in (3, 4, 5)}
    for _, g in df.groupby("plate"):
        p = np.c_[g.x.values, g.y.values]
        tree = cKDTree(p)
        nb = tree.query_ball_point(p, LMAX)
        for i, j in tree.query_pairs(LMAX, output_type="ndarray"):
            d = p[j] - p[i]; L = np.hypot(*d)
            if L < 5:
                continue
            u = d / L
            ks = [k for k in nb[i] if k != i and k != j]
            if not ks:
                continue
            v = p[ks] - p[i]
            t = v @ u / L
            h = np.abs(v[:, 0] * u[1] - v[:, 1] * u[0])
            mid = (t > 0) & (t < 1)
            for w in widths:
                m = int((mid & (h < w / 2)).sum())
                for r in (3, 4, 5):
                    if m >= r - 2:
                        out[(w, r)] += 1
    return out


data = count(cat)
nulls = []
for s in range(3):
    jit = cat.copy()
    rad = 60 * np.sqrt(rng.uniform(0, 1, len(jit))); ang = rng.uniform(0, 2 * np.pi, len(jit))
    jit["x"] += rad * np.cos(ang); jit["y"] += rad * np.sin(ang)
    nulls.append(count(jit))
print(f"pairs up to {LMAX/60:.0f}' long; counts are (i,j) endpoint pairs with >= r-2 points in the strip")
print(" width   r   data    null (mean of 3 jitters)   data/null")
for (w, r), v in data.items():
    nm = np.mean([n[(w, r)] for n in nulls])
    print(f" {w:4.0f}\"  {r}  {v:7d}   {nm:9.0f}                 {v/nm if nm else float('nan'):.2f}")
