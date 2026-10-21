"""Is there a glint signal in the track score distributions? A fit to the four distributions (3-6 transients) at once,
for glint tracks of each number of flashes in turn.

Background: each size's track score histogram (alignments.py, after arbitration within each size), averaged over the
1,000 random catalogues. Signal: for each number of flashes (3 to 6), a library of single artificial glint tracks
(alignments.with_tracks: that many flashes over 2-10', 1" scatter, any direction), each added to the real catalogue on
its own and its plate scored with and without it. A track's change to all four histograms is kept together: one glint
gives candidates of several sizes at once (a 6-flash glint also gives 5-, 4- and 3-transient candidates, after the same
arbitration as the search). The library's average t is the signal's shape, and carries the search's efficiency. The real
histograms are taken as Poisson counts with mean b + S t, S >= 0 the number of glint tracks, fitted by maximum
likelihood. The four histograms share transients across sizes and are not strictly Poisson, so every p-value and limit
comes from toys built from the random catalogues, never from asymptotic formulae.
- The background alone: the Poisson deviance of the real histograms from b, against the same for each random catalogue
  measured from the mean of the others.
- Evidence for a signal: q0 = 2 (ln L(S hat) - ln L(0)), against its values in the random catalogues.
- The limit, for each number of flashes: toys built like the data at each true number of tracks S, a random catalogue's
  histograms (the background, with its real fluctuations and correlations) plus S whole tracks from that library, fitted
  the same way. Random catalogues built from a catalogue that held real glint tracks would take in a little of them
  (alignment_null_checks.py, section 5), little of the strongest but much of the weakest, since each transient moves far
  less than a track is long: each toy's background carries what copies keep of its own tracks, track by track (library).
  CL(s+b) is the share of toys with S tracks that fit as little signal as the real catalogue.
  Bayesian (shown on the page): with every whole number of tracks equally likely beforehand, the posterior on S is
  CL(s+b) normalised, and the 95% credible limit is where its running sum reaches 95%. (CL(s+b) is the probability of
  fitting at most the real catalogue's S hat; with S hat = 0, as here, that is the probability of the fit itself, so the
  likelihood.) CLs (a check): CL(s+b) / CL(b), with S ruled out at 95% confidence when it falls below 0.05. A check on the
  plates' precision: the limits again with every glint's positions scattered by 2" rather than 1".

python analysis/alignment_fit.py -> data/alignment_fit.json"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import alignments as AL

ROOT = Path(__file__).resolve().parents[1]
FLASHES = (3, 4, 5, 6)
SMAX = 80                                                              # most glint tracks tried in the toys
NTOY = 20000
NLIB = 1500                                                            # single glint tracks in each library

Z = np.load(ROOT / "data/alignments.npz")
R, N = Z["real_l"].astype(float), Z["null_l"].astype(float)
K = len(N); tot = N.sum(0); b = tot / K
loo = lambda k: (tot - N[k]) / (K - 1)                                  # background for random catalogue k: the others
# a check on the library's model of what random catalogues keep: the 100 tracks of 4 to 6 transients added to the real
# catalogue and the mean of its 1,000 copies (alignment_null_checks.py, section 5), against the library's prediction
NCK = json.loads((ROOT / "data/alignment_null_checks.json").read_text())["absorption"]
KEPT = np.array(NCK["hist_injected_copies_mean"]) - N.mean(0)                         # what the copies took in, by bin
BANDS = [-8.0, -5.0, -4.0, -3.0, -2.0, 0.0]
band = lambda h: [round(float(h[:, (AL.BINS[:-1] >= lo) & (AL.BINS[:-1] < hi)].sum()), 2) for lo, hi in zip(BANDS[:-1], BANDS[1:])]


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


def fit(n, b_, t_):
    """S hat (at least 0) and q0."""
    s_hat = minimize_scalar(lambda s: nll(n, b_, t_, s), bounds=(0, 500), method="bounded", options={"xatol": 1e-3}).x
    s_hat = 0.0 if nll(n, b_, t_, 0) <= nll(n, b_, t_, s_hat) else float(s_hat)
    return s_hat, max(2 * (nll(n, b_, t_, 0) - nll(n, b_, t_, s_hat)), 0.0)


def limits(LIB, ABS, t, s_hat, rng):
    """CL(s+b) at each true number of tracks, until it is negligible; the posterior; the Bayesian and CLs limits."""
    clsb = []
    for st in range(SMAX + 1):
        hits = 0
        for k in rng.integers(0, K, NTOY):
            j = rng.integers(0, NLIB, st); sig = LIB[j].sum(0)
            hits += fit(np.clip(N[k] + sig, 0, None), loo(k) + ABS[j].sum(0), t)[0] <= s_hat + 1e-9
        clsb.append(hits / NTOY)
        if st > 0 and clsb[-1] < 1e-3:
            break
    clsb = np.array(clsb); post = clsb / clsb.sum(); cum = np.cumsum(post); cls = clsb / clsb[0]
    return clsb, post, int(np.flatnonzero(cum >= 0.95)[0]), int(np.flatnonzero(cls < 0.05)[0]) - 1


def deviance(n, b_):
    u = b_ > 0; n, b_ = n[u], b_[u]
    return float(2 * np.sum(b_ - n + np.where(n > 0, n * np.log(np.where(n > 0, n, 1) / b_), 0)))


if __name__ == "__main__":
    rng = np.random.default_rng(5)
    out = {"copies": K, "flashes": {}, "kept_check": {"bands": BANDS, "measured": band(KEPT), "library": None}}
    kept_lib = {}
    # the background alone, all four sizes and each
    d_null = np.array([deviance(N[k], loo(k)) for k in range(K)])
    out["background"] = {"deviance": round(deviance(R, b), 1), "copies_median": round(float(np.median(d_null)), 1),
                         "p": round(float((np.sum(d_null >= deviance(R, b)) + 1) / (K + 1)), 3),
                         "p_by_size": [round(float((np.sum(np.array([deviance(N[k][s], loo(k)[s]) for k in range(K)]) >= deviance(R[s], b[s])) + 1) / (K + 1)), 3)
                                       for s in range(4)]}
    for f in FLASHES:
        LIB, ABS = library(f); t = LIB.mean(0); kept_lib[f] = ABS.mean(0)
        # the fit, and how often the random catalogues give as much evidence for a signal
        s_hat, q0 = fit(R, b, t)
        q_null = np.array([fit(N[k], loo(k), t)[1] for k in range(K)])
        clsb, post, s95, s95_cls = limits(LIB, ABS, t, s_hat, rng)
        # the same with the glints' positions scattered by 2" rather than 1", as a check on the plates' precision
        LIB2, ABS2 = library(f, 2.0); s_hat2 = fit(R, b, LIB2.mean(0))[0]
        s95_2 = limits(LIB2, ABS2, LIB2.mean(0), s_hat2, rng)[2]
        out["flashes"][str(f)] = {"candidates": round(float(t.sum()), 2), "kept_by_copies": round(float(ABS.mean(0).sum()), 2), "s_hat": round(s_hat, 2), "q0": round(q0, 3),
                                  "p": round(float((np.sum(q_null >= q0) + 1) / (K + 1)), 3), "clsb": np.round(clsb, 4).tolist(),
                                  "p_post": np.round(post, 4).tolist(), "s95": s95, "s95_cls": s95_cls, "s95_scatter2": s95_2,
                                  "post_at_s95": round(float(np.cumsum(post)[s95]), 4)}
        print(f, {k_: v_ for k_, v_ in out["flashes"][str(f)].items() if k_ not in ("clsb", "p_post")}, flush=True)
    out["kept_check"]["library"] = band(100 * np.mean([kept_lib[f] for f in (4, 5, 6)], 0))     # 100 tracks, 4 to 6 transients
    print("kept by copies, by band:", out["kept_check"], flush=True)
    (ROOT / "data/alignment_fit.json").write_text(json.dumps(out, indent=1))
    print("background", out["background"])
