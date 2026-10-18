"""Alignments as a track-reconstruction problem.

Hits: the transients on a plate (positions in arcsec on the plate's tangent plane). Track candidates: for every pair of
hits up to 15' apart (the end points), the r - 2 hits lying between them that sit closest to the line joining them, for
r = 3..6, kept if every hit lies within 15" of that line and no two hits are closer than 10" along it (two-hit resolution). No direction cut: VASCO's |Pearson r| > 0.99 selection is not
applied. Each candidate is fitted with a straight line (total least squares) and scored by
    chi2 = sum(perpendicular residual^2) / sigma^2,   ndof = r - 2,
with sigma = 1" for the plates' positional precision. In log10(chi2/ndof) a different sigma only shifts the axis, so the
comparisons below do not depend on it.

Background (fake tracks): the same reconstruction on copies of the catalogue in which every transient moves by up to
20", which keeps the transients bunched where they are but breaks any genuine line. Signal templates: 100 artificial
glint tracks (4-6 hits, 2-10' long, with 1", 0.5" or 2" scatter) added to the real catalogue; the extra candidates they
produce.

python analysis/alignments_tracks.py   -> data/alignment_tracks.npz (histograms of log10(chi2/ndof) per catalogue, size
                                          and zone: the real catalogue, the catalogue with tracks added, 200 copies)."""
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
LMAX, HMAX = 900.0, 15.0             # arcsec: longest track, widest residual allowed
MINSEP = 10.0                        # arcsec: two-hit resolution along the track (the pipeline does not separate closer
                                     # sources, so shuffled copies would otherwise gain near-coincident hits the real one lacks)
INTERIOR = 2.5 * 3600                # arcsec from the plate centre; outside lies the edge band where defects pile up
EDGES = np.round(np.arange(-4.0, 3.0001, 0.1), 3)   # bins of log10(chi2 / ndof) with sigma = 1"
NSHUF = int(sys.argv[1]) if len(sys.argv) > 1 else 200
RS = (3, 4, 5, 6)

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
cat = cat[cat.plate != ""][["plate", "xi", "eta"]].assign(x=lambda d: d.xi * 3600, y=lambda d: d.eta * 3600)


@njit(cache=True)
def plate_tracks(P, pairs, nb_ptr, nb_idx, H, lmin, step, interior, minsep):
    """Fill H[r-3, zone, bin] with the candidates' log10(chi2/ndof). zone 0 = whole plate, 1 = inner part only."""
    nbin = H.shape[2]
    hk = np.empty(4); kk = np.empty(4, np.int64); tk = np.empty(4); sv = np.empty(6)
    for q in range(pairs.shape[0]):
        i, j = pairs[q, 0], pairs[q, 1]
        dx, dy = P[j, 0] - P[i, 0], P[j, 1] - P[i, 1]
        L = np.sqrt(dx * dx + dy * dy)
        if L < minsep:
            continue
        ux, uy = dx / L, dy / L
        nk = 0
        for a in range(nb_ptr[i], nb_ptr[i + 1]):
            k = nb_idx[a]
            if k == i or k == j:
                continue
            vx, vy = P[k, 0] - P[i, 0], P[k, 1] - P[i, 1]
            t = (vx * ux + vy * uy) / L
            if t <= 0.0 or t >= 1.0:
                continue
            h = abs(vx * uy - vy * ux)
            if h >= 15.0:
                continue
            # keep the four smallest offsets, sorted
            if nk < 4:
                pos = nk; nk += 1
            elif h < hk[3]:
                pos = 3
            else:
                continue
            while pos > 0 and hk[pos - 1] > h:
                hk[pos] = hk[pos - 1]; kk[pos] = kk[pos - 1]; tk[pos] = tk[pos - 1]; pos -= 1
            hk[pos] = h; kk[pos] = k; tk[pos] = t * L
        for r in range(3, 7):
            m = r - 2
            if nk < m:
                break
            # two-hit resolution: every pair of hits at least minsep apart along the track
            sv[0] = 0.0; sv[1] = L
            for a in range(m):
                sv[2 + a] = tk[a]
            srt = np.sort(sv[:r])
            ok = True
            for a in range(r - 1):
                if srt[a + 1] - srt[a] < minsep:
                    ok = False
            if not ok:
                continue
            # members: i, j and the m closest; straight-line fit by total least squares
            sx = P[i, 0] + P[j, 0]; sy = P[i, 1] + P[j, 1]
            for a in range(m):
                sx += P[kk[a], 0]; sy += P[kk[a], 1]
            cx, cy = sx / r, sy / r
            sxx = syy = sxy = 0.0
            inner = max(abs(P[i, 0]), abs(P[i, 1]), abs(P[j, 0]), abs(P[j, 1])) <= interior
            for a in range(-2, m):
                b = i if a == -2 else (j if a == -1 else kk[a])
                ex, ey = P[b, 0] - cx, P[b, 1] - cy
                sxx += ex * ex; syy += ey * ey; sxy += ex * ey
                if a >= 0 and max(abs(P[b, 0]), abs(P[b, 1])) > interior:
                    inner = False
            tr, det = sxx + syy, sxx * syy - sxy * sxy
            lam = 0.5 * (tr - np.sqrt(max(tr * tr - 4.0 * det, 0.0)))    # sum of squared perpendicular residuals
            v = np.log10(max(lam, 1e-12) / m)                            # sigma = 1", ndof = r - 2
            b_ = int((v - lmin) / step)
            if b_ < 0:
                b_ = 0
            if b_ >= nbin:
                continue
            H[r - 3, 0, b_] += 1
            if inner:
                H[r - 3, 1, b_] += 1


def histograms(df):
    H = np.zeros((len(RS), 2, len(EDGES) - 1), np.int64)
    for _, g in df.groupby("plate"):
        P = np.ascontiguousarray(np.c_[g.x.values, g.y.values])
        if len(P) < 3:
            continue
        tree = cKDTree(P)
        pairs = tree.query_pairs(LMAX, output_type="ndarray")
        if len(pairs) == 0:
            continue
        nb = tree.query_ball_point(P, LMAX)
        ptr = np.zeros(len(nb) + 1, np.int64); ptr[1:] = np.cumsum([len(n) for n in nb])
        idx = np.fromiter((k for n in nb for k in n), np.int64, ptr[-1])
        plate_tracks(P, pairs.astype(np.int64), ptr, idx, H, EDGES[0], 0.1, INTERIOR, MINSEP)
    return H


def jitter(seed, rad=20.0):
    rng = np.random.default_rng(seed)
    rr = rad * np.sqrt(rng.uniform(0, 1, len(cat))); a = rng.uniform(0, 2 * np.pi, len(cat))
    return cat.assign(x=cat.x.values + rr * np.cos(a), y=cat.y.values + rr * np.sin(a))


def with_tracks(n=100, seed=7, scatter=1.0):
    rng, rows, plates = np.random.default_rng(seed), [], cat.plate.unique()
    for _ in range(n):
        g = cat[cat.plate == rng.choice(plates)]
        r = int(rng.integers(4, 7)); L = rng.uniform(120, 600); th = rng.uniform(0, np.pi)
        x0, y0 = g.x.iloc[int(rng.integers(len(g)))], g.y.iloc[int(rng.integers(len(g)))]
        s = np.sort(rng.uniform(0, L, r)); perp = rng.normal(0, scatter, r); along = rng.normal(0, scatter, r)
        rows.append(pd.DataFrame({"plate": g.plate.iloc[0], "x": x0 + (s + along) * np.cos(th) - perp * np.sin(th),
                                  "y": y0 + (s + along) * np.sin(th) + perp * np.cos(th)}))
    return pd.concat([cat[["plate", "x", "y"]]] + rows)


def run(job):
    kind, seed = job
    df = cat if kind == "real" else with_tracks(scatter=float(kind[6:])) if kind.startswith("signal") else jitter(seed)
    return histograms(df)


if __name__ == "__main__":
    histograms(cat[cat.plate == cat.plate.iloc[0]])                  # compile once before forking
    jobs = [("real", 0), ("signal1.0", 0), ("signal0.5", 0), ("signal2.0", 0)] + [("shuffle", 3000 + s) for s in range(NSHUF)]
    with Pool(8) as pool:
        H = np.array(pool.map(run, jobs))                            # (catalogues, size, zone, bin)
    np.savez(ROOT / "data/alignment_tracks.npz", edges=EDGES, real=H[0], signal=H[1], signal_05=H[2], signal_2=H[3], null=H[4:])
    for z, name in ((0, "whole plates"), (1, "plate edges left out")):
        real, null = H[0, :, z].sum(0), H[4:, :, z].sum(1)
        print(f"{name}: real {real.sum():,} candidates, background {null.sum(1).mean():,.0f} "
              f"(+/- {null.sum(1).std():,.0f}); signal adds {(H[1, :, z] - H[0, :, z]).sum():,}")
