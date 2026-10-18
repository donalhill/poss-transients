"""Collinearity of POSS-I transients, self-calibrated against clustering.

For every pair (i, j) on the same plate with separation <= DMAX, every third transient k lying between them along the line
(0 < t < 1) contributes its perpendicular offset h from the line. If transients are merely clustered (no preference for
lines) dN/dh is flat for h much smaller than the clustering scale; genuine collinearity appears as an excess at h ~ 1-3".
Baseline: flat fit over 15-90". Orientation of the excess lines measured in the plate tangent frame (xi = east, eta = north).

Also: how often VASCO's Pearson |r(RA, Dec)| > 0.99 criterion passes a perfect line, as a function of its orientation."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
DMAX = 600.0   # arcsec, max length of the aligned segment (their boxes: a few to ~30')
HMAX = 90.0    # arcsec, perpendicular offsets recorded
rng = np.random.default_rng(1)

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
cat = cat[cat.plate != ""]
cat["x"], cat["y"] = cat.xi * 3600, cat.eta * 3600  # arcsec, x = east, y = north


def triples(x, y):
    """Return arrays (h, angle_deg, i, j) for all (i, j, k) with |ij| <= DMAX, k between i and j, h <= HMAX."""
    p = np.c_[x, y]
    tree = cKDTree(p)
    pairs = tree.query_pairs(DMAX, output_type="ndarray")
    if len(pairs) == 0:
        return np.empty(0), np.empty(0)
    nb = tree.query_ball_point(p, DMAX)
    H, A = [], []
    for i, j in pairs:
        d = p[j] - p[i]; L = np.hypot(*d)
        if L < 5:  # near-duplicates
            continue
        u = d / L
        ks = np.array([k for k in nb[i] if k != i and k != j])
        if ks.size == 0:
            continue
        v = p[ks] - p[i]
        t = v @ u / L
        h = np.abs(v[:, 0] * u[1] - v[:, 1] * u[0])
        m = (t > 0) & (t < 1) & (h <= HMAX)
        if m.any():
            H.append(h[m]); A.append(np.full(m.sum(), np.degrees(np.arctan2(u[1], u[0])) % 180))
    if not H:
        return np.empty(0), np.empty(0)
    return np.concatenate(H), np.concatenate(A)


hs, angs = [], []
for plate, g in cat.groupby("plate"):
    h, a = triples(g.x.values, g.y.values)
    hs.append(h); angs.append(a)
h = np.concatenate(hs); a = np.concatenate(angs)
np.savez_compressed(ROOT / "data/alignment_triples.npz", h=h, a=a)

edges = np.arange(0, HMAX + 1, 1.0)
n, _ = np.histogram(h, edges)
base = n[15:90].mean()
print(f"(i,j,k) configurations with h <= {HMAX:.0f}\": {len(h):,}; flat baseline {base:.1f} per 1\" bin (15-90\")")
for lo in range(0, 10):
    print(f"  h {lo}-{lo+1}\": {n[lo]:6d}  ratio to baseline {n[lo]/base:5.2f}")
exc = n[:3].sum() - 3 * base
print(f"excess in h < 3\": {exc:.0f} over expected {3*base:.0f} ({exc/np.sqrt(3*base):.1f} sigma Poisson; NB configurations share points)")

# orientation of aligned (h < 3") vs control (20 < h < 90") configurations, folded to 0-90 deg from E-W
def fold(x):
    x = x % 180
    return np.where(x > 90, 180 - x, x)
al, co = fold(a[h < 3]), fold(a[(h > 20) & (h < 90)])
bins = np.arange(0, 91, 15)
na, _ = np.histogram(al, bins); nc, _ = np.histogram(co, bins)
print("\norientation (deg from east-west): aligned share / control share, ratio")
for b0, b1, x1, x2 in zip(bins[:-1], bins[1:], na / na.sum(), nc / nc.sum()):
    print(f"  {b0:2d}-{b1:2d}: {x1:.3f} / {x2:.3f}  -> {x1/x2:.2f}")

# VASCO criterion: Pearson |r(RA, Dec)| > 0.99 for perfect 3-5 point lines with 1" scatter, length 5', vs orientation
print("\nPearson |r| > 0.99 pass rate for a 4-point line (5' long, 1\" scatter) vs orientation from east-west:")
for th in (0, 2, 5, 10, 20, 45, 70, 85, 88, 90):
    passes = 0
    for _ in range(2000):
        s = np.sort(rng.uniform(0, 300, 4))
        x = s * np.cos(np.radians(th)) + rng.normal(0, 1, 4)
        y = s * np.sin(np.radians(th)) + rng.normal(0, 1, 4)
        passes += abs(np.corrcoef(x, y)[0, 1]) > 0.99
    print(f"  {th:2d} deg: {passes/2000:.2f}")
