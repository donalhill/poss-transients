"""VASCO's own alignment score, applied to every alignment in the catalogue and in shuffled copies.

Villarroel et al. (2025, PASP) score each shortlisted alignment by the expected number of chance r-point alignments in a
box around it (Edmunds & George 1985). Their Table 2 is reproduced by
    mu_r = pi * N * n * (2 n p)^(r-2) * d^r / (r * (r-2)!),   n = N / A,
with N transients in a box of area A, strip half-width p and length d (candidates 1, 4* and 5*: 0.0436, 0.0033 and
0.00012 against their 0.044, 0.003 and 0.0001). Their best candidate (5*) has mu = 0.0001, the "1 in 10,000".

Here every r-point alignment (r = 3..6, up to 15' long) is found and scored the same way: p is the largest perpendicular
offset of its members from the line joining the end points (at least 1", their smallest strip), d is its length, and the
box is a square centred on it, max(10', 1.5 d) on a side, as in their table. Only groups passing VASCO's selection
(|Pearson r| > 0.99 between the members' two coordinates) are kept. Counts are also given with the plate edges (beyond
2.5 deg from the centre) left out, where artefacts pile up. The question: how many alignments in the
real catalogue score as improbable as VASCO's candidates, and how many do catalogues with no real lines produce?

python analysis/alignments_vasco_metric.py           -> data/alignments_vasco_scores.npz (histograms of log10 mu per
catalogue, track size and zone; the real catalogue first, then 48 copies with every transient moved by up to 20")."""
from math import factorial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
LMAX = 900.0                        # arcsec, longest alignment considered
RS = (3, 4, 5, 6)
NSHUF = 200
INTERIOR = 2.5 * 3600                # arcsec from the plate centre; outside this lies the edge zone where artefacts pile up

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
cat = cat[cat.plate != ""][["plate", "xi", "eta"]].assign(x=lambda d: d.xi * 3600, y=lambda d: d.eta * 3600)


def mu(r, N, A, p, d):
    n = N / A
    return np.pi * N * n * (2 * n * p) ** (r - 2) * d ** r / (r * factorial(r - 2))


BINS = np.arange(-8.0, 0.0001, 0.1)       # log10(mu); tracks with mu > 1 are unremarkable even by VASCO's measure


def score(df):
    """Histogram of log10(mu) for every track passing VASCO's selection, by size r and zone (whole plate / inner).
    Each track is counted once, from its unique pair of end points."""
    H = {(r, z): np.zeros(len(BINS) - 1, int) for r in RS for z in ("all", "inner")}
    for _, g in df.groupby("plate"):
        P = np.c_[g.x.values, g.y.values]
        tree = cKDTree(P)
        nb = tree.query_ball_point(P, LMAX)
        for i, j in tree.query_pairs(LMAX, output_type="ndarray"):
            dv = P[j] - P[i]; L = np.hypot(*dv)
            if L < 5:
                continue
            ks = np.array([k for k in nb[i] if k != i and k != j])
            if ks.size == 0:
                continue
            u = dv / L; v = P[ks] - P[i]
            t = v @ u / L; h = np.abs(v[:, 0] * u[1] - v[:, 1] * u[0])
            mid = (t > 0) & (t < 1) & (h < 30)
            if not mid.any():
                continue
            order = np.argsort(h[mid]); hk, kk = h[mid][order], ks[mid][order]
            side = max(600.0, 1.5 * L)                                  # arcsec
            N = len(tree.query_ball_point((P[i] + P[j]) / 2, side / 2, p=np.inf))
            for r in RS:
                if len(hk) < r - 2:
                    break
                m = mu(r, N, (side / 60.0) ** 2, max(hk[r - 3], 1.0) / 60.0, L / 60.0)
                if m > 1:
                    continue
                mem = np.r_[i, j, kk[:r - 2]]
                if abs(np.corrcoef(P[mem, 0], P[mem, 1])[0, 1]) <= 0.99:   # VASCO's selection: |Pearson r| > 0.99
                    continue
                b_ = min(max(int((np.log10(m) - BINS[0]) / 0.1), 0), len(BINS) - 2)
                H[(r, "all")][b_] += 1
                if np.max(np.abs(P[mem])) <= INTERIOR:
                    H[(r, "inner")][b_] += 1
    return H


def jitter(rad, seed):
    rng = np.random.default_rng(seed)
    rr = rad * np.sqrt(rng.uniform(0, 1, len(cat))); a = rng.uniform(0, 2 * np.pi, len(cat))
    return cat.assign(x=cat.x.values + rr * np.cos(a), y=cat.y.values + rr * np.sin(a))


def run(job):
    label, rad, seed = job
    return label, score(cat if rad == 0 else jitter(rad, seed))


if __name__ == "__main__":
    jobs = [("real", 0, 0)] + [(f"shuffle 20\" {s}", 20, 2000 + s) for s in range(NSHUF)]
    with Pool(8) as pool:
        res = pool.map(run, jobs)
    labels = [l for l, _ in res]
    out = {f"r{r}_{z}": np.array([H[(r, z)] for _, H in res]) for r in RS for z in ("all", "inner")}
    np.savez(ROOT / "data/alignments_vasco_scores.npz", bins=BINS, labels=np.array(labels), **out)
    for z in ("all", "inner"):
        print(f"--- {z} ---")
        for r in RS:
            A = out[f"r{r}_{z}"]; cum = np.cumsum(A, axis=1)          # tracks with log10(mu) below each bin edge
            for lvl in (-4, -3, -2, -1, 0):
                k = int(round((lvl - BINS[0]) / 0.1)) - 1
                real, null = cum[0, k], cum[1:, k]
                print(f"r={r} mu<=1e{lvl}: real {real:6d}  null mean {null.mean():8.1f}  share>=real {np.mean(null >= real):.2f}")
