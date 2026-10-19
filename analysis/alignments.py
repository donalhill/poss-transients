"""Alignments in VASCO's catalogue: the straight-line track fit and VASCO's own score, in one pass over the same track
candidates, the same copies of the catalogue and the same injected tracks.

Hits: the transients on a plate (positions in arcsec on the plate's tangent plane). Track candidates: for every pair of
hits up to 15' apart (the end points), the r - 2 hits lying between them that sit closest to the line joining them, for
r = 3..6, kept if every hit lies within 15" of that line and no two hits are closer than 8" along it (two-hit
resolution: the pipeline does not separate closer sources, so the copies would otherwise gain pairs the real one lacks).

Track fit: each candidate is fitted with a straight line (total least squares) and scored by
    chi2 = sum(perpendicular residual^2) / sigma^2,   ndof = r - 2,   sigma = 1" (the plates' positional precision).
No direction cut.

VASCO's score (Villarroel et al. 2025, PASP; Edmunds & George 1985): the expected number of chance r-point alignments in
a box around the candidate,
    mu_r = pi * N * n * (2 n p)^(r-2) * d^r / (r * (r-2)!),   n = N / A,
with N transients in a square box max(10', 1.5 d) on a side centred on the candidate, d its length and p the largest
offset of its members from the line joining the end points (at least 1", their narrowest strip; under 10", their widest).
This reproduces their Table 2 (candidates 1, 4* and 5*: 0.0436, 0.0033 and 0.00012 against their 0.044, 0.003 and
0.0001). Only candidates passing their selection, |Pearson r| > 0.99 between the members' two coordinates, and with
mu <= 1 are scored.

Track score: VASCO's expected number of chance alignments, built from the same box, but in every direction (no
correlation screen) and with the strip the widest residual from the fitted line rather than from the line through the
end points (at least 1", the plates' precision), and only for candidates straight to that precision: chi2 probability
at least 1% (with 1" per transient). Like the track fit, it is calibrated against the copies.

Isolation: the number of other transients within 2' of the segment between the end points (the 2D counterpart of the
occupancy of a cone around a reconstructed track); the track fit's distributions are also kept in bins of isolation.

Multiple candidates: one line of transients yields many overlapping candidates (every shorter line inside it). Each test
arbitrates its own candidates: best score first (chi2/ndof, or mu), and any candidate sharing a transient with one already
kept is dropped, within each size for the histograms; for VASCO's counts, across sizes (mu is comparable between sizes).
The real catalogue, every copy and every injected catalogue go through the same code.

Null: copies with every real line broken but the clustering kept: transients within 20" of one another (friends of
friends) move together as a rigid clump, rotated about its centre and shifted by up to 20", and lone transients are
shifted by up to 20" (alignments_vasco_null.py compares this with moving every transient independently). Signal: 100
artificial glint tracks (4-6 hits at least 8" apart, 2-10' long, with 1", 0.5" or 2" scatter, starting anywhere on a
plate) added to the real catalogue.

python analysis/alignments.py [ncopies] [independent]
    -> data/alignments.npz    histograms per catalogue and size, after arbitration: log10(chi2/ndof) (track fit) and
                              log10(mu) (VASCO's score), for the real catalogue, the catalogue with tracks added, and
                              the copies; and log10(mu) against log10(chi2/ndof) before arbitration
       data/alignments.json   VASCO's score: track candidates (all; arbitrated; arbitrated and straight, chi2/ndof <= 2)
                              scoring 1e-4 and 1e-3 or better, real catalogue against the copies; the real catalogue's
                              arbitrated candidates scoring 1e-3 or better
    with "independent", copies with every transient moved independently, to data/alignments_independent.npz."""
import json
import sys
from math import factorial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit
from scipy.spatial import cKDTree
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
LMAX = 900.0                         # arcsec: longest track candidate
HMAX, HMAX_V = 15.0, 10.0            # arcsec: widest offset from the line, track fit and VASCO's score
MINSEP = 8.0                         # arcsec: two-hit resolution along the track; the released catalogue has no two
                                     # transients on a plate closer than this, bar two pairs
SHIFT, LINK = 20.0, 20.0             # arcsec: largest move of a transient or clump in the copies; clump linking length
EDGES = np.round(np.arange(-4.0, 3.0001, 0.1), 3)   # log10(chi2/ndof) bins
BINS = np.round(np.arange(-8.0, 0.0001, 0.1), 3)    # log10(mu) bins; candidates with mu > 1 are unremarkable even to VASCO
RISO = 120.0                         # arcsec: isolation, other transients within this distance of the track
ISO = np.array([0, 1, 2, 3, 5, 10])  # isolation bins: 0, 1, 2, 3-4, 5-9, 10 or more other transients
CHI2MAX = np.r_[0.0, chi2.isf(0.01, np.arange(1, 5))]   # track score: chi2 probability at least 1%, by ndof (1-4)
KEEP, STRAIGHT = 1e-3, 2.0           # VASCO's counts: scores this low or better; chi2/ndof at or below which a candidate is straight
FACT = np.array([float(factorial(max(r - 2, 0))) for r in range(7)])   # (r-2)! by r
RS = (3, 4, 5, 6)
# injected catalogues, "signal:scatter:seed:direction": 1" scatter (ten seeds) in any direction, 0.5" and 2", and at fixed
# angles from east-west (three seeds each)
ANGLES = (0, 1, 2, 3, 5, 10, 20, 45, 70, 80, 85, 87, 88, 89, 90)
SIGNALS = ([f"signal:1.0:{7 + k}:any" for k in range(10)] + ["signal:0.5:7:any", "signal:2.0:7:any"]
           + [f"signal:1.0:{7 + k}:{a}" for a in ANGLES for k in range(3)])
NSIG = len(SIGNALS)
NSHUF = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 1000
INDEPENDENT = "independent" in sys.argv

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
cat = cat[cat.plate != ""][["plate", "xi", "eta"]].assign(x=lambda d: d.xi * 3600, y=lambda d: d.eta * 3600)


@njit(cache=True)
def plate_kernel(P, xs_order, xs, pairs, nb_ptr, nb_idx, minsep, hmax, hmax_v, riso, chi2max, out_i, out_f, out_mem):
    """Every track candidate on a plate: out_i (size, within VASCO's strip, passes VASCO's selection, isolation: other
    transients within riso of the segment between the end points), out_f (log10 chi2/ndof, log10 mu or +inf, strip p
    in arcsec, length in arcsec, log10 of the track score), out_mem (members, -1 padded).
    Returns the number written."""
    nout = 0
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
            if h >= hmax:
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
        N = -1                                                         # transients in VASCO's box, counted when first needed
        C = -1                                                         # transients within riso of the segment, likewise
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
            # isolation: transients within riso of the segment i-j (the members lie inside it)
            if C < 0:
                a0 = np.searchsorted(xs, min(P[i, 0], P[j, 0]) - riso, side="left")
                a1 = np.searchsorted(xs, max(P[i, 0], P[j, 0]) + riso, side="right")
                C = 0
                for a in range(a0, a1):
                    k = xs_order[a]
                    vx, vy = P[k, 0] - P[i, 0], P[k, 1] - P[i, 1]
                    t = (vx * ux + vy * uy) / L
                    if t < 0.0:
                        d = np.sqrt(vx * vx + vy * vy)
                    elif t > 1.0:
                        d = np.sqrt((P[k, 0] - P[j, 0]) ** 2 + (P[k, 1] - P[j, 1]) ** 2)
                    else:
                        d = abs(vx * uy - vy * ux)
                    if d <= riso:
                        C += 1
            # straight-line fit by total least squares
            sx = P[i, 0] + P[j, 0]; sy = P[i, 1] + P[j, 1]
            for a in range(m):
                sx += P[kk[a], 0]; sy += P[kk[a], 1]
            cx, cy = sx / r, sy / r
            sxx = syy = sxy = 0.0
            for a in range(-2, m):
                b = i if a == -2 else (j if a == -1 else kk[a])
                ex, ey = P[b, 0] - cx, P[b, 1] - cy
                sxx += ex * ex; syy += ey * ey; sxy += ex * ey
            tr, det = sxx + syy, sxx * syy - sxy * sxy
            lam = 0.5 * (tr - np.sqrt(max(tr * tr - 4.0 * det, 0.0)))    # sum of squared perpendicular residuals
            # transients in the box around the candidate (VASCO's box), for both scores
            if N < 0:
                side = max(600.0, 1.5 * L); hs = side / 2
                mx, my = (P[i, 0] + P[j, 0]) / 2, (P[i, 1] + P[j, 1]) / 2
                a0 = np.searchsorted(xs, mx - hs, side="left"); a1 = np.searchsorted(xs, mx + hs, side="right")
                N = 0
                for a in range(a0, a1):
                    if abs(P[xs_order[a], 1] - my) <= hs:
                        N += 1
                n = N / (side / 60.0) ** 2
            # VASCO's score, for candidates within their strip: widest offset from the line through the end points
            strip = hk[m - 1] < hmax_v
            lmu = np.inf
            if strip:
                mu = np.pi * N * n * (2.0 * n * max(hk[m - 1], 1.0) / 60.0) ** m * (L / 60.0) ** r / (r * FACT[r])
                lmu = np.log10(mu) if mu > 0 else -99.0
            # the track score: the same expected number, with the strip the widest residual from the fitted line
            phi = 0.5 * np.arctan2(2.0 * sxy, sxx - syy); nx_, ny_ = -np.sin(phi), np.cos(phi)
            w = 0.0
            for a in range(-2, m):
                b = i if a == -2 else (j if a == -1 else kk[a])
                w = max(w, abs((P[b, 0] - cx) * nx_ + (P[b, 1] - cy) * ny_))
            lam_ = np.pi * N * n * (2.0 * n * max(w, 1.0) / 60.0) ** m * (L / 60.0) ** r / (r * FACT[r])
            llam = np.log10(lam_) if lam_ > 0 else -99.0
            if lam > chi2max[m]:                                         # not straight to the plates' precision: no track
                llam = np.inf
            out_i[nout, 0] = r; out_i[nout, 1] = strip; out_i[nout, 3] = C - r
            out_i[nout, 2] = not (sxx > 0.0 and syy > 0.0 and abs(sxy / np.sqrt(sxx * syy)) <= 0.99)
            out_f[nout, 0] = np.log10(max(lam, 1e-12) / m); out_f[nout, 1] = lmu; out_f[nout, 4] = llam
            out_f[nout, 2] = max(hk[m - 1], 1.0); out_f[nout, 3] = L
            out_mem[nout, 0] = i; out_mem[nout, 1] = j
            for a in range(4):
                out_mem[nout, 2 + a] = kk[a] if a < m else -1
            nout += 1
    return nout


@njit(cache=True)
def arbitrate(order, mem, r, npts):
    """Best first (order), keep a candidate only if none of its transients is already used by a kept one."""
    used = np.zeros(npts, np.bool_); keep = np.zeros(len(r), np.bool_)
    for q in order:
        ok = True
        for a in range(r[q]):
            if used[mem[q, a]]:
                ok = False; break
        if ok:
            keep[q] = True
            for a in range(r[q]):
                used[mem[q, a]] = True
    return keep


def binned(v, edges):
    return np.clip(np.floor((v - edges[0]) / 0.1).astype(int), 0, len(edges) - 2)


def score(df, detail=False):
    """Both tests on one catalogue: Ht[r-3, chi2 bin] and Hv[r-3, mu bin] after arbitration within each size, and
    Hi[r-3, isolation bin, chi2 bin] for the track fit's kept candidates;
    H2[r-3, screen, mu bin, chi2 bin] before arbitration (screen 0 = within VASCO's strip, 1 = also passing their
    selection); and VASCO's candidates scoring KEEP or better, before arbitration."""
    Ht = np.zeros((len(RS), len(EDGES) - 1), np.int64); Hv = np.zeros((len(RS), len(BINS) - 1), np.int64)
    Hi = np.zeros((len(RS), len(ISO), len(EDGES) - 1), np.int64); Hl = np.zeros((len(RS), len(BINS) - 1), np.int64)
    both = []                                                          # (log10 mu, log10 track score, size) where both score
    keep_t = []                                                        # track-score candidates scoring KEEP or better
    H2 = np.zeros((len(RS), 2, len(BINS) - 1, len(EDGES) - 1), np.int64); keep = []
    for plate, g in df.groupby("plate"):
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
        order = np.argsort(P[:, 0], kind="stable"); n = len(pairs) * 4
        out_i, out_f, out_mem = np.zeros((n, 4), np.int64), np.zeros((n, 5)), np.zeros((n, 6), np.int64)
        k = plate_kernel(P, order, P[order, 0].copy(), pairs.astype(np.int64), ptr, idx, MINSEP, HMAX, HMAX_V, RISO, CHI2MAX, out_i, out_f, out_mem)
        r, strip, pear, iso = out_i[:k, 0], out_i[:k, 1].astype(bool), out_i[:k, 2].astype(bool), out_i[:k, 3]
        lchi, lmu, llam, mem = out_f[:k, 0], out_f[:k, 1], out_f[:k, 4], out_mem[:k]
        scored = strip & (lmu <= 0)                                    # within VASCO's strip, mu <= 1
        vasco = scored & pear
        for s in RS:
            sel = np.flatnonzero(r == s)
            kept = arbitrate(sel[np.argsort(lchi[sel], kind="stable")], mem, r, len(P))
            inr = kept & (lchi < EDGES[-1])
            np.add.at(Ht[s - 3], binned(lchi[inr], EDGES), 1)
            np.add.at(Hi[s - 3], (np.searchsorted(ISO, iso[inr], side="right") - 1, binned(lchi[inr], EDGES)), 1)
            sel = np.flatnonzero((r == s) & vasco)
            kept = arbitrate(sel[np.argsort(lmu[sel], kind="stable")], mem, r, len(P))
            np.add.at(Hv[s - 3], binned(lmu[kept], BINS), 1)
            sel = np.flatnonzero((r == s) & (llam <= 0))
            kept = arbitrate(sel[np.argsort(llam[sel], kind="stable")], mem, r, len(P))
            np.add.at(Hl[s - 3], binned(llam[kept], BINS), 1)
        if detail:
            both.append(np.c_[lmu, llam, r, vasco, np.isfinite(llam) & (llam <= 0)])
        for scr, msk in ((0, scored), (1, vasco)):
            np.add.at(H2, (r[msk] - 3, scr, binned(lmu[msk], BINS), binned(lchi[msk], EDGES)), 1)
        for c in np.flatnonzero(np.isfinite(llam) & (llam <= np.log10(KEEP))):
            m_ = mem[c][mem[c] >= 0]
            keep_t.append({"plate": plate, "ids": g.index.values[m_].tolist(), "r": int(r[c]), "mu": float(10 ** llam[c]), "iso": int(iso[c]),
                           "chi2ndof": float(10 ** lchi[c]), "d_arcmin": float(out_f[c, 3] / 60), "vasco": bool(vasco[c]),
                           "x": P[m_, 0].tolist(), "y": P[m_, 1].tolist()})
        for c in np.flatnonzero(vasco & (lmu <= np.log10(KEEP))):
            m_ = mem[c][mem[c] >= 0]
            keep.append({"plate": plate, "ids": g.index.values[m_].tolist(), "r": int(r[c]), "mu": float(10 ** lmu[c]), "iso": int(iso[c]),
                         "chi2ndof": float(10 ** lchi[c]), "p_arcsec": float(out_f[c, 2]), "d_arcmin": float(out_f[c, 3] / 60),
                         "x": P[m_, 0].tolist(), "y": P[m_, 1].tolist()})
    return (Ht, Hv, H2, Hi, Hl, keep, keep_t) + ((np.concatenate(both),) if detail else ())


def compare(d):
    """VASCO's score against the track score, candidate by candidate, in the real catalogue (before arbitration): the
    joint histogram where both score, their rank correlation, how many of each one's 50 strongest the other has among
    its 50 strongest, and the candidates only one of them scores (track score: not straight to 1"; VASCO: their
    correlation screen, or outside their 10" strip)."""
    lmu, llam, r, vasco, track = d[:, 0], d[:, 1], d[:, 2], d[:, 3].astype(bool), d[:, 4].astype(bool)
    both = vasco & track
    from scipy.stats import spearmanr
    i_v, i_t = set(np.argsort(np.where(both, lmu, np.inf))[:50]), set(np.argsort(np.where(both, llam, np.inf))[:50])
    return {"ml_hist": np.histogram2d(lmu[both], llam[both], [BINS, BINS])[0].astype(np.int64),
            "ml_rho": spearmanr(lmu[both], llam[both])[0], "ml_rho_46": spearmanr(lmu[both & (r >= 4)], llam[both & (r >= 4)])[0],
            "ml_top50": len(i_v & i_t), "ml_counts": np.array([vasco.sum(), track.sum(), both.sum()]),
            "ml_vasco_strong": np.array([(vasco & (lmu <= -3)).sum(), (vasco & (lmu <= -3) & ~track).sum()])}


def arbitrated(keep, lvl):
    """VASCO's candidates of any size scoring lvl or better, arbitrated across sizes: best score first, and any candidate
    sharing a transient with one already kept is dropped."""
    used, out = set(), []
    for a in sorted((a for a in keep if a["mu"] <= lvl), key=lambda a: a["mu"]):
        if used.isdisjoint(a["ids"]):
            out.append(a); used.update(a["ids"])
    return out


LEVELS_T = (-3.0, -4.0, -4.5, -5.0)   # log10 track score: levels for the counts of the strongest tracks


def tally_t(keep_t):
    """Track-score candidates arbitrated across sizes (best score first, each transient once), counted at each level."""
    A_ = arbitrated(keep_t, 10 ** LEVELS_T[0])
    return [int(sum(np.log10(a["mu"]) <= lv for a in A_)) for lv in LEVELS_T]


TALLY = [("all candidates", "1e-4"), ("all candidates", "1e-3"), ("candidates", "1e-4"), ("candidates", "1e-3"),
         ("straight candidates", "1e-4"), ("straight candidates", "1e-3")]


def tally(keep):
    A4, A3 = arbitrated(keep, 1e-4), arbitrated(keep, 1e-3)
    return ([sum(a["mu"] <= lvl for a in keep) for lvl in (1e-4, 1e-3)] + [len(A4), len(A3)]
            + [sum(a["chi2ndof"] <= STRAIGHT for a in A_) for A_ in (A4, A3)])


@njit(cache=True)
def clump_labels(n, pairs):
    """Friends-of-friends groups from linked pairs (union-find); labels 0..G-1."""
    parent = np.arange(n)
    for q in range(pairs.shape[0]):
        i, j = pairs[q, 0], pairs[q, 1]
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        while parent[j] != j:
            parent[j] = parent[parent[j]]; j = parent[j]
        if i != j:
            parent[i] = j
    lab = np.empty(n, np.int64); remap = -np.ones(n, np.int64); g = 0
    for i in range(n):
        k = i
        while parent[k] != k:
            k = parent[k]
        if remap[k] < 0:
            remap[k] = g; g += 1
        lab[i] = remap[k]
    return lab, g


def clumps(df, link=LINK):
    off = pd.factorize(df.plate)[0] * 1e7                             # keeps plates apart in one tree
    return clump_labels(len(df), cKDTree(np.c_[df.x.values + off, df.y.values]).query_pairs(link, output_type="ndarray").astype(np.int64))


def clump_copy(df, seed, link=LINK, shift=SHIFT):
    """A copy with every real line broken but the clustering kept: transients within link of one another (friends of
    friends, never across plates) form a clump, each clump rotated about its centre by a random angle and moved by up
    to shift as a rigid unit; lone transients are just moved."""
    rng = np.random.default_rng(seed)
    x, y = df.x.values, df.y.values
    lab, G = clumps(df, link)
    cnt = np.bincount(lab, minlength=G); cx = np.bincount(lab, x, G) / cnt; cy = np.bincount(lab, y, G) / cnt
    th = rng.uniform(0, 2 * np.pi, G); rr = shift * np.sqrt(rng.uniform(0, 1, G)); a = rng.uniform(0, 2 * np.pi, G)
    c, s = np.cos(th)[lab], np.sin(th)[lab]; dx, dy = x - cx[lab], y - cy[lab]
    return df.assign(x=cx[lab] + c * dx - s * dy + (rr * np.cos(a))[lab], y=cy[lab] + s * dx + c * dy + (rr * np.sin(a))[lab])


def independent_copy(df, seed, shift=SHIFT):
    """A copy with every transient moved independently by up to shift."""
    rng = np.random.default_rng(seed)
    rr = shift * np.sqrt(rng.uniform(0, 1, len(df))); a = rng.uniform(0, 2 * np.pi, len(df))
    return df.assign(x=df.x.values + rr * np.cos(a), y=df.y.values + rr * np.sin(a))


def with_tracks(n=100, seed=7, scatter=1.0, angle=None):
    """The real catalogue with n artificial glint tracks added (4-6 hits, 2-10' long, given scatter; any direction, or at
    a given angle in degrees from east-west, 0 being the track a geostationary object would leave). Tracks start anywhere
    on a plate's searched area (within 2.9 deg of the centre), as glints from orbit would, and their hits are at least
    MINSEP apart along the track, as the plates can resolve."""
    rng, rows, plates, nid = np.random.default_rng(seed), [], cat.plate.unique(), 10_000_000
    for _ in range(n):
        plate = rng.choice(plates)
        r = int(rng.integers(4, 7)); L = rng.uniform(120, 600); th = rng.uniform(0, np.pi) if angle is None else np.radians(angle) * rng.choice([-1, 1])
        x0, y0 = rng.uniform(-2.9 * 3600 + L, 2.9 * 3600 - L, 2)
        s = np.sort(rng.uniform(0, L, r))
        while np.min(np.diff(s)) < MINSEP:
            s = np.sort(rng.uniform(0, L, r))
        perp = rng.normal(0, scatter, r); along = rng.normal(0, scatter, r)
        rows.append(pd.DataFrame({"plate": plate, "x": x0 + (s + along) * np.cos(th) - perp * np.sin(th),
                                  "y": y0 + (s + along) * np.sin(th) + perp * np.cos(th)}, index=range(nid, nid + r))); nid += r
    return pd.concat([cat[["plate", "x", "y"]]] + rows)


def run(job):
    kind, seed = job
    if kind == "real":
        df = cat
    elif kind.startswith("signal"):
        sc, sd, d_ = kind.split(":")[1:]
        df = with_tracks(seed=int(sd), scatter=float(sc), angle=None if d_ == "any" else float(d_))
    else:
        df = independent_copy(cat, seed) if INDEPENDENT else clump_copy(cat, seed)
    Ht, Hv, H2, Hi, Hl, keep, keep_t, *det = score(df, detail=kind == "real")
    return (Ht, Hv, H2, tally(keep), (keep if kind == "real" else None), Hi, Hl, (det[0] if det else None),
            tally_t(keep_t), (keep_t if kind == "real" else None))


if __name__ == "__main__":
    score(cat[cat.plate == cat.plate.iloc[0]])                       # compile once before forking
    jobs = [("real", 0)] + [(k_, 0) for k_ in SIGNALS] + [("copy", 3000 + s) for s in range(NSHUF)]
    with Pool(8) as pool:
        R = pool.map(run, jobs)
    Ht = np.array([r_[0] for r_ in R]); Hv = np.array([r_[1] for r_ in R]); T = np.array([r_[3] for r_ in R[NSIG + 1:]])
    Hi = np.array([r_[5] for r_ in R]); Hl = np.array([r_[6] for r_ in R]); c0 = NSIG + 1
    Tt = np.array([r_[8] for r_ in R[c0:]]); tt_real, keep_t = R[0][8], R[0][9]
    out = ROOT / ("data/alignments_independent.npz" if INDEPENDENT else "data/alignments.npz")
    np.savez(out, edges=EDGES, bins=BINS, iso=ISO, signals=np.array(SIGNALS),
             real_t=Ht[0], signal_t=Ht[1:c0], null_t=Ht[c0:], real_v=Hv[0], signal_v=Hv[1:c0], null_v=Hv[c0:],
             real_l=Hl[0], signal_l=Hl[1:c0], null_l=Hl[c0:], real_i=Hi[0], signal_i=Hi[1:c0], null_i=Hi[c0:],
             real_h2=R[0][2], null_h2_mean=np.mean([r_[2] for r_ in R[c0:]], 0), **compare(R[0][7]))
    print(f"{out.name}: track fit, real {Ht[0].sum():,} candidates after arbitration, copies {Ht[c0:].sum((1, 2)).mean():,.0f}; "
          f"VASCO's score, real {Hv[0].sum():,}, copies {Hv[c0:].sum((1, 2)).mean():,.0f}; track score, real {Hl[0].sum():,}, "
          f"copies {Hl[c0:].sum((1, 2)).mean():,.0f}")
    t_real, keep = R[0][3], R[0][4]; counts = {}
    for c, (kind, lvl) in enumerate(TALLY):
        counts[f"{kind.replace(' ', '_')}_{lvl}"] = {"real": t_real[c], "copies_mean": round(float(T[:, c].mean()), 2),
                                                     "share_copies_ge_real": round(float(np.mean(T[:, c] >= t_real[c])), 3)}
        print(f"VASCO's score, {kind} scoring {lvl} or better: real {t_real[c]}, copies {T[:, c].mean():.2f}; share of copies with as many {np.mean(T[:, c] >= t_real[c]):.3f}")
    track = {"levels": list(LEVELS_T), "real": tt_real, "copies_mean": Tt.mean(0).round(3).tolist(),
             "share_copies_ge_real": [round(float(np.mean(Tt[:, j] >= tt_real[j])), 3) for j in range(len(LEVELS_T))],
             "candidates": [{"plate": a["plate"], "r": a["r"], "score": a["mu"], "chi2ndof": round(a["chi2ndof"], 2), "d_arcmin": round(a["d_arcmin"], 2),
                             "iso": a["iso"], "vasco": a["vasco"], "x": [round(v, 2) for v in a["x"]], "y": [round(v, 2) for v in a["y"]]}
                            for a in arbitrated(keep_t, 10 ** LEVELS_T[0])]}
    for j, lv in enumerate(LEVELS_T):
        print(f"track score, strongest tracks (across sizes) at 10^{lv} or better: real {tt_real[j]}, copies {Tt[:, j].mean():.2f}; share of copies with as many {np.mean(Tt[:, j] >= tt_real[j]):.3f}")
    if not INDEPENDENT:
        rows = [{"plate": a["plate"], "r": a["r"], "mu": a["mu"], "p_arcsec": round(a["p_arcsec"], 1), "d_arcmin": round(a["d_arcmin"], 2),
                 "chi2ndof": round(a["chi2ndof"], 2), "xi": round(float(np.mean(a["x"]) / 3600), 4), "eta": round(float(np.mean(a["y"]) / 3600), 4),
                 "x": [round(v, 2) for v in a["x"]], "y": [round(v, 2) for v in a["y"]]} for a in arbitrated(keep, KEEP)]
        (ROOT / "data/alignments.json").write_text(json.dumps({"counts": counts, "candidates": rows, "track": track, "ncopies": NSHUF,
                                                               "shift_arcsec": SHIFT, "link_arcsec": LINK, "hmax_vasco_arcsec": HMAX_V}))
