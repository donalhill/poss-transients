"""Two checks on the nuclear-test excess, run after the frozen analysis (NUCLEAR_DESIGN.md, addendum).

1. Blue plates. Specks (sources with no Gaia star within 5", plate_blue.py) counted on the scans of the red plate and
   of the blue plate of the same field, taken the same night. The within-run comparison of nuclear_fit.py is run on
   each. Something in the emulsion, or settling on the plates, would add specks to both; transients seen only in red
   light, as VASCO's are, would not add them to the blue plates.
2. Prediction out of sample. Each observing night held out in turn: the within-run model is fitted to every other
   plate, the held-out night's run level comes from its run's other nights, and its plates are predicted with and
   without the test-date feature. The gain in held-out log-likelihood is compared with the gain from the same feature
   given to random nights within runs.

python analysis/nuclear_checks.py -> data/nuclear_checks.json"""
import json

import numpy as np
import pandas as pd
from scipy.special import gammaln

from nuclear_model import GAP, P, ROOT, Trials, day, full_fit, pk, r3, reassign, runs, window, y_raw

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

(ROOT / "data/nuclear_checks.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
