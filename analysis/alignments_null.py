"""Null distribution for the alignment counts of analysis/alignments_counts.py, from many shuffled copies of the catalogue.

Each shuffled copy moves every transient to a random point within JIT arcseconds of where it was (default 60"). That
keeps the bunching of transients on scales above an arcminute and breaks any alignment within the 2-15" strips being
counted. It also smooths the bunching on 20-60" scales, which adds a few per cent of chance alignments: the real
catalogue sits below every 60" copy. A 20" shuffle (python analysis/alignments_null.py 20) keeps that bunching and
lands on the real counts, but leaves part of any genuine 10-15" alignment intact, so the 60" copies are the reference.
Counts are also split by the length of the pair that defines the strip (shorter or longer than 60"), to see where any
difference between the real catalogue and the shuffled copies comes from.

Closure check: 100 artificial alignments (4-5 points, 1" scatter, 2-10' long) are added to the real catalogue; the count
on the affected plates should rise, and after shuffling the excess should be gone.
Writes data/alignment_null.csv (or alignment_null_<JIT>.csv) (one row per catalogue: real, shuffle k) and prints the closure check."""
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
LMAX, NSHUF = 900.0, 48
JIT = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
KINDS = [(3, 2.0), (4, 5.0), (5, 10.0), (5, 15.0)]   # (points in the strip, full width in arcsec)

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
cat = cat[cat.plate != ""][["plate", "xi", "eta"]].assign(x=lambda d: d.xi * 3600, y=lambda d: d.eta * 3600)


def count(df):
    out = {f"{r}_{w:g}_{band}": 0 for r, w in KINDS for band in ("short", "long")}
    for _, g in df.groupby("plate"):
        p = np.c_[g.x.values, g.y.values]
        tree = cKDTree(p)
        nb = tree.query_ball_point(p, LMAX)
        for i, j in tree.query_pairs(LMAX, output_type="ndarray"):
            d = p[j] - p[i]; L = np.hypot(*d)
            if L < 5:
                continue
            ks = [k for k in nb[i] if k != i and k != j]
            if not ks:
                continue
            u = d / L; v = p[ks] - p[i]
            t = v @ u / L
            h = np.abs(v[:, 0] * u[1] - v[:, 1] * u[0])
            mid = (t > 0) & (t < 1)
            band = "short" if L < 60 else "long"
            for r, w in KINDS:
                if int((mid & (h < w / 2)).sum()) >= r - 2:
                    out[f"{r}_{w:g}_{band}"] += 1
    return out


def jitter(df, seed):
    rng = np.random.default_rng(seed)
    rad = JIT * np.sqrt(rng.uniform(0, 1, len(df))); ang = rng.uniform(0, 2 * np.pi, len(df))
    return df.assign(x=df.x.values + rad * np.cos(ang), y=df.y.values + rad * np.sin(ang))


def shuffled(seed):
    return {"catalogue": f"shuffle {seed}", **count(jitter(cat, 1000 + seed))}


def injected():
    rng, rows = np.random.default_rng(7), []
    plates = cat.plate.unique()
    for _ in range(100):
        g = cat[cat.plate == rng.choice(plates)]
        r = rng.integers(4, 6); L = rng.uniform(120, 600); th = rng.uniform(0, np.pi)
        x0, y0 = g.x.iloc[int(rng.integers(len(g)))], g.y.iloc[int(rng.integers(len(g)))]
        s = np.sort(rng.uniform(0, L, r)); perp = rng.normal(0, 1.0, r)
        rows.append(pd.DataFrame({"plate": g.plate.iloc[0], "x": x0 + s * np.cos(th) - perp * np.sin(th),
                                  "y": y0 + s * np.sin(th) + perp * np.cos(th)}))
    art = pd.concat(rows)
    base = cat[cat.plate.isin(art.plate.unique())]
    both = pd.concat([base[["plate", "x", "y"]], art])
    tot = lambda c: {f"{r}_{w:g}": c[f"{r}_{w:g}_short"] + c[f"{r}_{w:g}_long"] for r, w in KINDS}
    before = {k: tot(count(df)) for k, df in (("real", base), ("with lines", both))}
    after = {k: tot(count(jitter(df, 99))) for k, df in (("real", base), ("with lines", both))}
    print(f"closure check on the {base.plate.nunique()} plates that received artificial alignments:")
    for r, w in KINDS:
        k = f"{r}_{w:g}"
        print(f"  {r} within {w:g}\": added lines raise the count by {before['with lines'][k] - before['real'][k]:+d}; "
              f"after shuffling, by {after['with lines'][k] - after['real'][k]:+d}")


if __name__ == "__main__":
    with Pool(8) as pool:
        res = [{"catalogue": "real", **count(cat)}] + pool.map(shuffled, range(NSHUF))
    df = pd.DataFrame(res)
    df.to_csv(ROOT / ("data/alignment_null.csv" if JIT == 60 else f"data/alignment_null_{JIT:g}.csv"), index=False)
    real, null = df.iloc[0], df.iloc[1:]
    for r, w in KINDS:
        for band in ("short", "long", None):
            cols = [f"{r}_{w:g}_{b}" for b in (("short", "long") if band is None else (band,))]
            rv, nv = real[cols].sum(), null[cols].sum(axis=1)
            print(f"{r} within {w:>2g}\" {band or 'all':5s}: real {rv:6d}  shuffled {nv.mean():8.0f} +/- {nv.std():5.0f}  "
                  f"(range {nv.min()}-{nv.max()})")
    injected()
