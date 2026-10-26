"""Is there a glint signal in the track score distributions? A fit to the four distributions (3-6 transients) at once,
for glint tracks of each number of flashes in turn.

Background: each size's track score histogram (alignments.py, after arbitration within each size), averaged over the
1,000 random catalogues. Signal: for each number of flashes (3 to 6), a library of single artificial glint tracks
(alignments.with_tracks: that many flashes over 2-10', 1" scatter, any direction), each added to the real catalogue on
its own and its plate scored with and without it. A track's change to all four histograms is kept together: one glint
gives candidates of several sizes at once (a 6-flash glint also gives 5-, 4- and 3-transient candidates, after the same
arbitration as the search). The library's average t is the signal's shape. The real histograms are taken as Poisson
counts with mean b + S t, S >= 0 the number of glint tracks, fitted by maximum likelihood. The four histograms share
transients across sizes and are not strictly Poisson, so every p-value comes from the random catalogues, never from
asymptotic formulae:
- The background alone: the Poisson deviance of the real histograms from b, against the same for each random catalogue
  measured from the mean of the others.
- Evidence for a signal: q0 = 2 (ln L(S hat) - ln L(0)), against its values in the random catalogues.
Every flash of a library track is in the catalogue, so S counts glint tracks whose every flash the catalogue holds; with
no measure of how often the catalogue's selection keeps a real flash, no limit on the number of glint tracks is set.
POSS_FIT_FLASHES=3 runs one number of flashes (-> alignment_fit_part3.json; alignment_fit_merge.py joins the parts).
POSS_SAMPLE=DIR python analysis/alignment_fit.py -> DIR/alignment_fit.json"""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

import alignments as AL

FLASHES = (3, 4, 5, 6)
NLIB = 1500                                                            # single glint tracks in each library

OUT = AL.OUT
Z = np.load(OUT / "alignments.npz")
R, N = Z["real_l"].astype(float), Z["null_l"].astype(float)
K = len(N); tot = N.sum(0); b = tot / K
loo = lambda k: (tot - N[k]) / (K - 1)                                  # background for random catalogue k: the others


def library(flashes, scatter=1.0):
    """Each single glint track's change to the four track score histograms: its plate scored with and without it. And
    what a random catalogue keeps of it: a copy of the plate made the page's way, scored with and without the track's
    own transients moved as the copies move them (a transient by up to 20", transients within 20" of one another as one
    piece). Each transient moves by far less than a track is long, so the moved transients can still line up weakly."""
    base, L, A = {}, [], []
    for i in range(NLIB):
        df = AL.with_tracks(n=1, seed=100_000 * flashes + i, scatter=scatter, hits=(flashes, flashes)); q = df.plate.iloc[-1]
        if q not in base:
            base[q] = AL.score(AL.cat[AL.cat.plate == q][["plate", "x", "y"]])[4]
        L.append(AL.score(df[df.plate == q])[4] - base[q])
        cp = AL.clump_copy(AL.cat[AL.cat.plate == q][["plate", "x", "y"]], 7_000_000 + 10 * flashes + i)
        tr = AL.clump_copy(df.iloc[-flashes:], 8_000_000 + 10 * flashes + i)
        A.append(AL.score(pd.concat([cp, tr]))[4] - AL.score(cp)[4])
    return np.array(L, float), np.array(A, float)


def nll(n, b_, t_, s):
    mu = np.maximum(b_ + s * t_, 1e-12)
    return float(np.sum(mu - n * np.log(mu)))


def s_hat_many(Nm, Bm, t, tol=1e-10):
    """The maximum-likelihood number of tracks, at least 0, for many histograms at once (rows of Nm with backgrounds Bm,
    all bins flattened): the Poisson likelihood of b + S t is concave in S, so the zero of its derivative is found by
    safeguarded Newton (a Newton step whenever it stays inside the bracket [lo, hi] on [0, 500], else bisection),
    converging quadratically to far below the toys' precision."""
    t = np.ravel(t)
    def g_h(s_):                                                       # the score and its (negative) derivative
        mu = np.maximum(Bm + s_[:, None] * t[None, :], 1e-12)
        return (t[None, :] * (Nm / mu - 1)).sum(1), -(t[None, :] ** 2 * Nm / mu ** 2).sum(1)
    lo, hi = np.zeros(len(Nm)), np.full(len(Nm), 500.0)
    g0 = g_h(lo)[0]; zero = g0 <= 0; s_ = np.full(len(Nm), 1.0)
    for _ in range(60):
        g, h = g_h(s_)
        lo = np.where(g > 0, s_, lo); hi = np.where(g > 0, hi, s_)
        nt = s_ - g / np.where(h < 0, h, -1e-12)
        s_new = np.where((nt > lo) & (nt < hi), nt, (lo + hi) / 2)
        if np.all(np.abs(s_new - s_) < tol):
            s_ = s_new; break
        s_ = s_new
    return np.where(zero, 0.0, s_)


def fit(n, b_, t_):
    """S hat (at least 0) and q0."""
    s_hat = float(s_hat_many(np.ravel(n)[None, :], np.ravel(b_)[None, :], t_)[0])
    return s_hat, max(2 * (nll(n, b_, t_, 0) - nll(n, b_, t_, s_hat)), 0.0)


def deviance(n, b_):
    u = b_ > 0; n, b_ = n[u], b_[u]
    return float(2 * np.sum(b_ - n + np.where(n > 0, n * np.log(np.where(n > 0, n, 1) / b_), 0)))


def library_cached(f):
    """library(f), cached beside the sample (rebuilt when the search's output is newer than the cache, or the code has
    changed)."""
    import hashlib
    c = OUT / f"alignment_library_{f}.npz"
    key = hashlib.sha256(b"".join(Path(p).read_bytes() for p in (__file__, AL.__file__))).hexdigest()
    if c.exists() and c.stat().st_mtime > (OUT / "alignments.npz").stat().st_mtime:
        z = np.load(c)
        if "key" in z and str(z["key"]) == key:
            return z["LIB"], z["ABS"]
    LIB, ABS = library(f)
    np.savez_compressed(c, LIB=LIB, ABS=ABS, key=key); return LIB, ABS


def run_size(f):
    """One number of flashes: the fit and its significance against the copies."""
    LIB, ABS = library_cached(f); t = np.tensordot(np.full(len(LIB), 1 / len(LIB)), LIB, 1)   # every track counted equally
    s_hat, q0 = fit(R, b, t); q_null = np.array([fit(N[k], loo(k), t)[1] for k in range(K)])
    return {"candidates": round(float(t.sum()), 2), "kept_by_copies": round(float(ABS.mean(0).sum()), 2), "s_hat": round(s_hat, 2), "q0": round(q0, 3),
            "p": round(float((np.sum(q_null >= q0) + 1) / (K + 1)), 3), "tracks": int(len(LIB))}


if __name__ == "__main__":
    only = os.environ.get("POSS_FIT_FLASHES"); sizes = [int(only)] if only else list(FLASHES)
    out = {"copies": K, "flashes": {}}
    if not only or only == "3":
        d_null = np.array([deviance(N[k], loo(k)) for k in range(K)])
        out["background"] = {"deviance": round(deviance(R, b), 1), "copies_median": round(float(np.median(d_null)), 1),
                             "p": round(float((np.sum(d_null >= deviance(R, b)) + 1) / (K + 1)), 3),
                             "p_by_size": [round(float((np.sum(np.array([deviance(N[k][s_], loo(k)[s_]) for k in range(K)]) >= deviance(R[s_], b[s_])) + 1) / (K + 1)), 3)
                                           for s_ in range(4)]}
        if "real_axis" in Z:                                     # angles of the strongest candidates from the scan's axes
            na = Z["null_axis"]; out["axis_angles"] = {"bins_deg": list(range(0, 46, 5)), "real": Z["real_axis"].tolist(),
                                                       "copies_mean": np.round(na.mean(0), 2).tolist(),
                                                       "p_first_bin": round(float((np.sum(na[:, 0] >= Z["real_axis"][0]) + 1) / (len(na) + 1)), 3)}
    for f in sizes:
        out["flashes"][str(f)] = run_size(f)
        print(f, out["flashes"][str(f)], flush=True)
    (OUT / (f"alignment_fit_part{only}.json" if only else "alignment_fit.json")).write_text(json.dumps(out, indent=1))
