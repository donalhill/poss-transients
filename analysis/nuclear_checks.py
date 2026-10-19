"""Two checks on the nuclear-test excess, run after the frozen analysis (NUCLEAR_DESIGN.md, addendum).

1. Blue plates. Specks (sources with no Gaia star within 5", plate_blue.py) counted on the scans of the red plate and
   of the blue plate of the same field, taken the same night. The within-run comparison of nuclear_fit.py is run on
   each. Something in the emulsion, or settling on the plates, would add specks to both; transients seen only in red
   light, as VASCO's are, would not add them to the blue plates.
2. Prediction out of sample. Each observing night held out in turn: the within-run model is fitted to every other
   plate, the held-out night's run level comes from its run's other nights, and its plates are predicted with and
   without the test-date feature. The gain in held-out log-likelihood is compared with the gain from the same feature
   given to random nights within runs.
3. Tracks. Each plate's number of straight track candidates (alignments.py's track score at 1 in 10, and at 1 in 100,
   or better; each transient in at most one track, best first), with the plate's number of transients as a further
   covariate: do transients near tests line up more often than their numbers imply?
4. The richest nights (a diagnostic found after the frozen results): each night's transients over what the model
   expects of them; how many nights are as rich as the richest near tests, and how many of those fall near tests; and
   r without the two richest nights near tests; and each informative observing run on its own.

python analysis/nuclear_checks.py -> data/nuclear_checks.json"""
import json

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import hypergeom

from nuclear_model import DAY0, GAP, P, ROOT, Trials, day, full_fit, pk, r3, reassign, runs, window, y_raw

NREP, NRAND = 2000, 200
run0 = runs(GAP)
S0 = window()
out = {}

# 1. specks on the red and blue plates of the same fields
B = pd.read_csv(ROOT / "data/plate_blue.csv").set_index("plate").reindex(P.index)
keep = (B.patches == 2).values & (B.blue_date == B.red_date).values
out["blue"] = {"plates": int(keep.sum()), "plates_near": int((keep & (S0 == 1)).sum())}
for colour in ("red", "blue"):
    y = B[f"{colour}_specks"].fillna(0).values.astype(float)
    r, se, a, dof, mu = full_fit(y, run0, S0, keep=keep)
    o = np.ones(len(y)); o[keep] = mu
    A = Trials(o, run0, a, dof, keep=keep).r(reassign(S0, run0, NREP, 5), y)
    out["blue"][colour] = {"r": r3(r), "count_lo": r3(r * np.exp(-se)), "count_hi": r3(r * np.exp(se)), "p_runs": r3(pk(A, r)),
                           "specks_per_plate": round(float(y[keep].mean()), 1)}
# the same plates' VASCO transients, for comparison on exactly this subset
r, se, a, dof, mu = full_fit(y_raw, run0, S0, keep=keep)
o = np.ones(len(y_raw)); o[keep] = mu
out["blue"]["catalogue"] = {"r": r3(r), "p_runs": r3(pk(Trials(o, run0, a, dof, keep=keep).r(reassign(S0, run0, NREP, 6), y_raw), r))}


# 2. prediction out of sample: each night held out in turn
def nb_ll(y, mu, a):
    """Negative binomial log-likelihood of counts y with means mu and scatter a."""
    k = 1 / a
    return gammaln(y + k) - gammaln(k) - gammaln(y + 1) + k * np.log(k / (k + mu)) + y * np.log(mu / (k + mu))


_, _, alpha, dof, mu0 = full_fit(y_raw, run0, np.zeros(len(y_raw)))      # the plate properties, with no test feature
features = {"VASCO's window": S0, "the day before": window((-1,)), "the test's date": window((0,)), "the day after": window((1,))}
out["prediction"] = []
for name, S in features.items():
    F = np.c_[S, reassign(S, run0, NRAND, 9)]                          # the real feature, then features on random nights
    gain = np.zeros(F.shape[1]); n_nights = 0
    for d in np.unique(day):
        test = day == d; train = ~test
        if not (train & (run0 == run0[test][0])).any():                  # the night's run has no other night
            continue
        n_nights += 1
        g, lam, a = Trials(mu0, run0, alpha, dof, keep=train).r(F, y_raw, levels=True)
        u = np.unique(run0[train]); j = np.searchsorted(u, run0[test][0])
        base = mu0[test][:, None] * np.exp(lam[j])[None, :]
        with_f = base * np.where(F[test] == 1, g[None, :], 1.0)
        # without the feature: the run level refitted with no feature at all
        _, lam0, a0 = Trials(mu0, run0, alpha, dof, keep=train).r(np.zeros((len(y_raw), 1)), y_raw, levels=True)
        without = mu0[test] * np.exp(lam0[j, 0])
        gain += (nb_ll(y_raw[test][:, None], with_f, a[None, :]) - nb_ll(y_raw[test], without, a0[0])[:, None]).sum(0)
    out["prediction"].append({"feature": name, "plates": int(S.sum()), "nights_held_out": n_nights, "gain": r3(gain[0]),
                              "random_median": r3(np.median(gain[1:])), "p": r3(pk(gain[1:], gain[0]))})

# 3. tracks per plate
import alignments as AL                                                       # noqa: E402
from scipy.spatial import cKDTree                                             # noqa: E402


def tracks(level):
    """Each plate's straight track candidates scoring 10**level or better, each transient in at most one, best first."""
    n = pd.Series(0, index=P.index)
    for plate, g in AL.cat.groupby("plate"):
        Q = np.ascontiguousarray(np.c_[g.x.values, g.y.values])
        if len(Q) < 3:
            continue
        tree = cKDTree(Q); pairs = tree.query_pairs(AL.LMAX, output_type="ndarray")
        if len(pairs) == 0:
            continue
        nb = tree.query_ball_point(Q, AL.LMAX)
        ptr = np.zeros(len(nb) + 1, np.int64); ptr[1:] = np.cumsum([len(x) for x in nb])
        idx = np.fromiter((k for x in nb for k in x), np.int64, ptr[-1])
        order = np.argsort(Q[:, 0], kind="stable"); m = len(pairs) * 4
        oi, of, om = np.zeros((m, 4), np.int64), np.zeros((m, 5)), np.zeros((m, 6), np.int64)
        k = AL.plate_kernel(Q, order, Q[order, 0].copy(), pairs.astype(np.int64), ptr, idx, AL.MINSEP, AL.HMAX, AL.HMAX_V, AL.RISO, AL.CHI2MAX, oi, of, om)
        llam = of[:k, 4]; sel = np.flatnonzero(np.isfinite(llam) & (llam <= level))
        n[plate] = int(AL.arbitrate(sel[np.argsort(llam[sel], kind="stable")], om[:k], oi[:k, 0], len(Q)).sum())
    return n.values.astype(float)


out["tracks"] = []
logn = np.log(np.maximum(y_raw, 1))[:, None]
for level, lab in ((-1, "1 in 10"), (-2, "1 in 100")):
    y = tracks(level)
    r, se, a, dof, mu = full_fit(y, run0, S0, extra=logn)
    A = Trials(mu, run0, a, dof).r(reassign(S0, run0, NREP, 13), y)
    out["tracks"].append({"score": lab, "tracks": int(y.sum()), "tracks_near": int(y[S0 == 1].sum()), "r": r3(r),
                          "count_lo": r3(r * np.exp(-se)), "count_hi": r3(r * np.exp(se)), "p_runs": r3(pk(A, r))})

# 4. the richest nights
r0, _, a0, dof0, mu0w = full_fit(y_raw, run0, S0)
N = pd.DataFrame({"day": day, "n": y_raw, "e": mu0w, "S": S0}).groupby("day").agg(n=("n", "sum"), e=("e", "sum"), S=("S", "max"))
N["ratio"] = N.n / N.e
top2 = N[N.S == 1].sort_values("ratio", ascending=False).head(2)
cut = round(float(top2.ratio.min()), 1) - 0.3                         # nights about as rich as the second-richest near tests
rich = N[N.ratio >= cut]
keep2 = ~np.isin(day, top2.index.values)
rk, sek, ak, dofk, muk = full_fit(y_raw, run0, S0, keep=keep2)
o = np.ones(len(y_raw)); o[keep2] = muk
Ak = Trials(o, run0, ak, dofk, keep=keep2).r(reassign(S0 * keep2, run0, NREP, 17), y_raw)
date = lambda d: str((DAY0 + pd.Timedelta(days=int(d))).date())
out["rich_nights"] = {"top_near": [{"date": date(d), "ratio": r3(r_.ratio), "n": int(r_.n), "plates": list(P.index[day == d])} for d, r_ in top2.iterrows()],
                      "threshold": cut, "nights": int(len(N)), "nights_near": int(N.S.sum()), "rich": int(len(rich)), "rich_near": int(rich.S.sum()),
                      "rich_elsewhere": [date(d) for d in rich.index[rich.S == 0]],
                      "p_chance": r3(hypergeom.sf(int(rich.S.sum()) - 1, len(N), int(N.S.sum()), len(rich))),
                      "r_without": r3(rk), "p_without": r3(pk(Ak, rk))}
# each informative run on its own: plates near tests against the run's other plates, each against its expected count
# from the plate properties alone, nothing shared between runs; with 68% counting ranges
import patsy                                                                  # noqa: E402
import statsmodels.api as sm                                                  # noqa: E402
from scipy.stats import gamma as gamma_dist                                   # noqa: E402
import plate_properties as PP                                                 # noqa: E402
Xp = np.asarray(patsy.dmatrix(PP.FORMULA.split("~")[1].replace(" + bs(t, df=16)", ""), PP.standardised(P)))
e_props = sm.GLM(y_raw, Xp, family=sm.families.NegativeBinomial(alpha=a0)).fit().fittedvalues
out["runs"] = []
for rn in np.unique(run0[S0 == 1]):
    k = run0 == rn; w, o = k & (S0 == 1), k & (S0 == 0)
    if not o.any():
        continue
    nw, ew, no, eo = y_raw[w].sum(), e_props[w].sum(), y_raw[o].sum(), e_props[o].sum()
    c = (nw / ew) / (no / eo)
    # counting range on the contrast, from the share of the run's transients on the plates near tests
    f = nw / (nw + no); lo_f, hi_f = gamma_dist.ppf(0.16, nw) / (gamma_dist.ppf(0.16, nw) + no), gamma_dist.ppf(0.84, nw + 1) / (gamma_dist.ppf(0.84, nw + 1) + no)
    to_c = lambda x: (x / (1 - x)) * (eo / ew)
    out["runs"].append({"from": date(day[k].min()), "to": date(day[k].max()), "near": int(w.sum()), "others": int(o.sum()),
                        "n_near": int(nw), "n_others": int(no), "near_ratio": r3(nw / ew), "others_ratio": r3(no / eo),
                        "contrast": r3(c), "lo": r3(to_c(lo_f)), "hi": r3(to_c(hi_f))})
cs = np.array([x["contrast"] for x in out["runs"]])
out["runs_summary"] = {"runs": len(cs), "above": int((cs > 1).sum()), "within_10pct": int((np.abs(cs - 1) < 0.1).sum()), "median": r3(np.median(cs))}

(ROOT / "data/nuclear_checks.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
