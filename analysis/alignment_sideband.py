"""The sideband of bent candidates: a check of the copies (the background-only model) on the real catalogue where a glint
cannot appear. A glint's track is straight to the plates' precision, so candidates failing the straightness test (chi2
probability below 1%, with 1" per transient) hold only chance alignments, built from the same pairs, clumps, and local
density as the chance straight candidates. Straight candidates are dropped before anything is counted.
For each version of the copies (LINK = BOX, SHIFT) and two bands of straightness (chi2 probability 1e-3 to 1e-2, just
failing; 1e-6 to 1e-3), the bent candidates of each size, arbitrated within the band (best track score first, each
transient used once), in bins of log10 track score (computed as for a straight candidate), in the real catalogue and in
NCOPY copies; and the counts at 1e-3, 1e-2, and 1e-1 or better.
POSS_SAMPLE=DIR python analysis/alignment_sideband.py -> DIR/alignment_sideband.json"""
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).parent))
sys.argv = sys.argv[:1]
import alignments as AL                                                      # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
cat = AL.cat[["plate", "x", "y"]]
CONFIGS = ((20.0, 20.0), (20.0, 12.0), (25.0, 20.0), (25.0, 12.0), (30.0, 20.0))   # (LINK = BOX, SHIFT)
BANDS = ((1e-3, 1e-2), (1e-6, 1e-3))                                     # chi2 probability
NCOPY = int(os.environ.get("NCOPY", "100"))
LCUT = (-3.0, -2.0, -1.0)


def bent(df):
    """Hist[band, size, log10 track score bin] of bent candidates, arbitrated within each band and size."""
    H = np.zeros((len(BANDS), len(AL.RS), len(AL.BINS) - 1), np.int64)
    for plate, g in df.groupby("plate"):
        P = np.ascontiguousarray(np.c_[g.x.values, g.y.values])
        if len(P) < 3:
            continue
        tree = cKDTree(P); pairs = tree.query_pairs(AL.LMAX, output_type="ndarray")
        if len(pairs) == 0:
            continue
        src, dst = np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]]
        o = np.lexsort((dst, src)); idx = dst[o].astype(np.int64)
        ptr = np.zeros(len(P) + 1, np.int64); ptr[1:] = np.cumsum(np.bincount(src, minlength=len(P)))
        order = np.argsort(P[:, 0], kind="stable"); n = len(pairs) * 4
        lab = np.r_[AL.clump_labels(len(P), tree.query_pairs(AL.LINK, output_type="ndarray").astype(np.int64))[0], -1]
        buf = [a_ for _ in range(2) for a_ in (np.zeros((n, 4), np.int64), np.zeros((n, 6)), np.zeros((n, 6), np.int64))]
        _, kb = AL.plate_kernel(P, order, P[order, 0].copy(), pairs.astype(np.int64), ptr, idx, AL.MINSEP, AL.BOX, AL.HMAX, AL.HMAX_V, AL.RISO, AL.CHI2MAX, *buf)
        oi, of, om = buf[3][:kb], buf[4][:kb], buf[5][:kb]
        L = lab[om]; ok = ((L[:, :, None] == L[:, None, :]) & (L[:, None, :] >= 0)).sum(2).max(1) < 3
        r, ndof = oi[:, 0], oi[:, 0] - 2
        p = chi2.sf(10 ** of[:, 0] * ndof, ndof)                          # chi2 probability, 1" per transient
        raw = of[:, 5]
        straight = np.isfinite(of[:, 4])                                  # the signal region: never counted
        for b, (lo, hi) in enumerate(BANDS):
            inb = ok & ~straight & (p >= lo) & (p < hi) & np.isfinite(raw)
            for s in AL.RS:
                sel = np.flatnonzero(inb & (r == s))
                kept = AL.arbitrate(sel[np.argsort(raw[sel], kind="stable")], om, r, len(P))
                np.add.at(H[b, s - 3], AL.binned(raw[kept], AL.BINS), 1)
    return H


def job(j):
    (link, shift), seed = j; AL.LINK, AL.BOX = link, link
    df = cat if seed < 0 else AL.clump_copy(cat, seed, link=link, shift=shift)
    return (link, shift), seed, bent(df)


def at(H, lv):                                                            # counts at lv or better, all sizes
    return int(H[:, :np.searchsorted(AL.BINS, lv)].sum())


if __name__ == "__main__":
    AL.score(cat[cat.plate == cat.plate.iloc[0]])
    with Pool(8) as pool:
        R = pool.map(job, [(c, -1) for c in CONFIGS] + [(c, 3000 + k) for c in CONFIGS for k in range(NCOPY)])
    out = {"bands_chi2_prob": BANDS, "copies": NCOPY, "bins_log10": AL.BINS.tolist(), "configs": {}}
    for c in CONFIGS:
        real = next(H for c_, s, H in R if c_ == c and s < 0); C = np.array([H for c_, s, H in R if c_ == c and s >= 0])
        res = {}
        for b in range(len(BANDS)):
            rr = np.array([at(real[b], lv) for lv in LCUT]); cc = np.array([[at(Ck[b], lv) for lv in LCUT] for Ck in C])
            res[f"band{b}"] = {"real": rr.tolist(), "copies_mean": cc.mean(0).round(1).tolist(), "copies_sd": cc.std(0).round(1).tolist(),
                               "ratio": (rr / cc.mean(0)).round(3).tolist(),
                               "share_copies_ge_real": ((cc >= rr).mean(0)).round(3).tolist(),
                               "real_hist": real[b].tolist(), "copies_hist_mean": C[:, b].mean(0).round(2).tolist()}
        out["configs"][f"link{int(c[0])}_shift{int(c[1])}"] = res
    out["levels_log10"] = LCUT
    (D / "alignment_sideband.json").write_text(json.dumps(out, indent=1))
    for k, v in out["configs"].items():
        print(k, {b: {q: v[b][q] for q in ("real", "copies_mean", "ratio", "share_copies_ge_real")} for b in v})
