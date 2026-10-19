"""Nuclear tests from the ground up: comparisons a reader can check by hand, each adding one thing to the last.

1. Raw: transients per plate on plates taken within a day of a US test (VASCO's window), against all other plates.
2. Within observing runs: each plate near a test compared only with the other plates of its own run (a run is a block
   of consecutive observing nights, a new one after more than 4 nights with no plates). The runs are pooled with the
   Mantel-Haenszel rate ratio, the textbook way of combining comparisons made within separate groups.
3. As 2, with each plate's count measured against what its properties predict (depth, sharpness, pointing, dust,
   exposure, airmass, the Moon, quality flag, ecliptic latitude; plate_properties.py), so plates taken in better
   conditions are not mistaken for an effect.
4. How often chance does as well: in each run, the nights counted as near a test reassigned at random among the run's
   nights; and the test calendar slid by whole years and up to 20 days.
5. Other windows: every run of days within 3, 7 and 14 days of a test, and windows in hours around each shot, each
   scored the same way; the global p is the share of slid calendars whose best window does as well as the data's best.

Each run is also shown on its own, with a range from the scatter between plates (not just counting).

6. A 95% upper limit on the ratio in 3, from effects injected into the trials of 4.
7. ABCD within plates: each plate's transients split into two classes (edge or central; crowded or isolated), and the
   ratio of class shares on plates in the window against the other plates of their runs. Anything that acts on a whole
   plate cancels; a sky signal would favour central, isolated objects, plate damage the edges and crowds.

python analysis/nuclear_basic.py -> data/nuclear_basic.json"""
import json

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.stats import norm

import plate_properties as PP
from nuclear_model import DAY0, GAP, P, ROOT, SHIFTS, day, hours, pk, r3, reassign, runs, window, y_raw

NREP = 20000
R = 14
HOURS = {"0 to 24 hours after a shot": (0, 24), "24 to 48 hours after": (24, 48), "24 hours before a shot": (-24, 0)}
run0 = runs(GAP)
S0 = window()
date = lambda d: str((DAY0 + pd.Timedelta(days=int(d))).date())

# each plate's expected count from its properties alone (no run levels, no window), and the scatter between plates
X = np.asarray(patsy.dmatrix(PP.FORMULA.split("~")[1].replace(" + bs(t, df=16)", ""), PP.standardised(P)))
a = 0.5
for _ in range(4):
    mu = sm.GLM(y_raw, X, family=sm.families.NegativeBinomial(alpha=a)).fit().fittedvalues
    a = float(np.sum(((y_raw - mu) ** 2 - mu) / mu ** 2) / (len(y_raw) - X.shape[1]))
e_props = np.asarray(mu)
G = (run0[None, :] == np.unique(run0)[:, None]).astype(float)           # runs x plates


def mh(S, e):
    """The Mantel-Haenszel rate ratio over runs, for each column of S: within each run, transients on plates in the
    window against those on the run's other plates, each over its expected count (e)."""
    S = np.atleast_2d(np.asarray(S, float).T).T
    a1, b1 = G @ (y_raw[:, None] * S), G @ (y_raw[:, None] * (1 - S))
    t1, t0 = G @ (e[:, None] * S), G @ (e[:, None] * (1 - S))
    t = t1 + t0
    num = np.sum(np.divide(a1 * t0, t, out=np.zeros_like(t), where=t > 0), 0)
    den = np.sum(np.divide(b1 * t1, t, out=np.zeros_like(t), where=t > 0), 0)
    return np.divide(num, den, out=np.full(S.shape[1], np.nan), where=den > 0)


ones = np.ones(len(y_raw))
out = {"plates_near": int(S0.sum()), "nights_near": int(len(np.unique(day[S0 == 1]))), "alpha": r3(a)}
# 1. raw
out["raw"] = r3(y_raw[S0 == 1].mean() / y_raw[S0 == 0].mean())
# 2. within runs, plates counted alike
out["within_runs"] = r3(mh(S0, ones)[0])
# 3. within runs, against each plate's expected count
r_hat = float(mh(S0, e_props)[0])
out["within_runs_expected"] = r3(r_hat)

# 4. chance: nights reassigned within runs; the calendar slid by whole years
T = reassign(S0, run0, NREP, 31)
A = mh(T, e_props)
SLIDES = [s for s in SHIFTS if window(shift=s).sum() >= 5 and (window(shift=s) * S0).sum() < 0.2 * window(shift=s).sum()]
B = mh(np.stack([window(shift=s) for s in SLIDES], 1), e_props)
q = lambda d: {k: r3(v) for k, v in zip(("q16", "median", "q84"), np.quantile(d[np.isfinite(d)], [0.16, 0.5, 0.84]))}
out["null_runs"] = {**q(A), "n": NREP, "p": r3(pk(A, r_hat)), "sigma": r3(norm.isf(pk(A, r_hat))),
                    "hist": np.histogram(A[np.isfinite(A)], np.linspace(0, 3, 61))[0].tolist()}
out["null_calendar"] = {**q(B), "n": len(SLIDES), "p": r3(pk(B, r_hat))}

# each informative run on its own, with a range from the scatter between plates: var(sum n) = sum(mu + a mu^2)
out["runs"] = []
for rn in np.unique(run0[S0 == 1]):
    k = run0 == rn; w, o = k & (S0 == 1), k & (S0 == 0)
    if not o.any():
        continue
    nw, ew, no, eo = y_raw[w].sum(), e_props[w].sum(), y_raw[o].sum(), e_props[o].sum()
    c = (nw / ew) / (no / eo)
    # each side's mean scaled to the run's own level, for its scatter
    lev = (nw + no) / (ew + eo)
    vw, vo = np.sum(lev * e_props[w] + a * (lev * e_props[w]) ** 2), np.sum(lev * e_props[o] + a * (lev * e_props[o]) ** 2)
    sl = np.sqrt(vw / (lev * ew) ** 2 + vo / (lev * eo) ** 2)
    out["runs"].append({"from": date(day[k].min()), "to": date(day[k].max()), "near": int(w.sum()), "others": int(o.sum()),
                        "n_near": int(nw), "n_others": int(no), "near_ratio": r3(nw / ew), "others_ratio": r3(no / eo),
                        "contrast": r3(c), "lo": r3(c * np.exp(-sl)), "hi": r3(c * np.exp(sl)),
                        "nights_near": [date(d) for d in np.unique(day[w])]})
cs = np.array([x["contrast"] for x in out["runs"]])
out["runs_summary"] = {"runs": len(cs), "above": int((cs > 1).sum()), "within_10pct": int((np.abs(cs - 1) < 0.1).sum()),
                       "median": r3(np.median(cs))}
# each run left out
out["leave_run_out"] = []
for rn in [rn for rn in np.unique(run0[S0 == 1]) if ((run0 == rn) & (S0 == 0)).any()]:
    S_ = S0 * (run0 != rn)
    out["leave_run_out"].append({"run": date(day[run0 == rn].min()), "r": r3(mh(S_, e_props)[0])})

# 5. other windows, by the calendar slides
WIN = [(a_, b_) for a_ in range(-R, R + 1) for b_ in range(a_, R + 1) if window(tuple(range(a_, b_ + 1))).sum() >= 5]
cols = lambda s: np.stack([window(tuple(range(a_, b_ + 1)), shift=s) for a_, b_ in WIN], 1)
real_w = mh(cols(0), e_props)
null_w = np.array([mh(cols(s), e_props) for s in SLIDES])
null_w = np.where(np.isfinite(null_w), null_w, -np.inf)
loc = np.array([pk(null_w[:, i], real_w[i]) for i in range(len(WIN))])
rank = (len(SLIDES) - np.argsort(np.argsort(null_w, axis=0), axis=0)) / (len(SLIDES) + 1)
span = np.array([max(abs(a_), abs(b_)) for a_, b_ in WIN]); iP = WIN.index((-1, 1))
out["scan"] = {"windows": len(WIN), "primary_p_local": r3(loc[iP]), "primary_rank": int(np.sum(loc <= loc[iP])),
               "grid": [{"a": a_, "b": b_, "r": r3(real_w[i]), "p": r3(loc[i]), "z": round(float(np.clip(norm.isf(loc[i]), -4, 4)), 2)}
                        for i, (a_, b_) in enumerate(WIN)], "by_size": []}
for days_ in (3, 7, 14):
    c_ = np.flatnonzero(span <= days_); b_ = c_[np.argmin(loc[c_])]
    pg = (np.sum(rank[:, c_].min(1) <= loc[b_]) + 1) / (len(SLIDES) + 1)
    out["scan"]["by_size"].append({"days": days_, "windows": len(c_), "best": list(WIN[b_]), "best_r": r3(real_w[b_]),
                                   "best_p_local": r3(loc[b_]), "p_global": r3(pg), "sigma": r3(norm.isf(pg))})
out["physical"] = []
for lab, (lo, hi) in HOURS.items():
    S_ = hours(lo, hi); rs = float(mh(S_, e_props)[0])
    As = mh(reassign(S_, run0, NREP, 41), e_props)
    out["physical"].append({"window": lab, "plates": int(S_.sum()), "r": r3(rs), "p_runs": r3(pk(As, rs))})

# 6. a 95% upper limit: effects injected into the within-run reassignments (r > 1 by thinning the other plates)
rng = np.random.default_rng(51)
def injected(rt, n=2000):
    T_ = reassign(S0, run0, n, int(rt * 1000)); Y = np.repeat(y_raw[:, None], n, 1).astype(int)
    Y = np.where(T_ == 0, rng.binomial(Y, 1 / rt), Y) if rt > 1 else np.where(T_ == 1, rng.binomial(Y, rt), Y)
    a1, b1 = G @ (Y * T_), G @ (Y * (1 - T_)); t1, t0 = G @ (e_props[:, None] * T_), G @ (e_props[:, None] * (1 - T_)); t = t1 + t0
    num = np.sum(np.divide(a1 * t0, t, out=np.zeros_like(t), where=t > 0), 0); den = np.sum(np.divide(b1 * t1, t, out=np.zeros_like(t), where=t > 0), 0)
    return num / den
grid = np.round(np.arange(1.0, 4.01, 0.05), 2)
low = np.array([np.mean(injected(rt) <= r_hat) for rt in grid])        # how often a true rt would give a value as low as ours
upper = float(grid[np.argmax(low < 0.05)]) if (low < 0.05).any() else None
out["upper95"] = upper
out["vasco_p_as_low"] = r3(np.interp(1.69, grid, low))

# 7. ABCD within plates
from scipy.spatial import cKDTree                                        # noqa: E402
from nuclear_model import cat                                            # noqa: E402
edge = (np.maximum(cat.xi.abs(), cat.eta.abs()) > 2.5).values
crowd = np.zeros(len(cat), bool)
for _, idx in cat.groupby("plate").indices.items():
    g = cat.iloc[idx]; xy = np.c_[g.xi, g.eta] * 60
    if len(xy) > 1:
        d, _ = cKDTree(xy).query(xy, k=2); crowd[idx] = d[:, 1] < 5
pidx = pd.Series(np.arange(len(P)), index=P.index)[cat.plate].values
out["abcd"] = []
for lab, cls in (("edge against central", edge), ("crowded against isolated", crowd)):
    nb = np.bincount(pidx, weights=cls, minlength=len(P)); ns = y_raw - nb    # background-like, signal-like per plate
    def share_ratio(S):
        """Mantel-Haenszel over runs of the odds that a transient is signal-like, plates in the window against the rest."""
        S = np.atleast_2d(np.asarray(S, float).T).T
        A_, B_ = G @ (ns[:, None] * S), G @ (nb[:, None] * S); C_, D_ = G @ (ns[:, None] * (1 - S)), G @ (nb[:, None] * (1 - S))
        tot = A_ + B_ + C_ + D_
        num = np.sum(np.divide(A_ * D_, tot, out=np.zeros_like(tot), where=tot > 0), 0)
        den = np.sum(np.divide(B_ * C_, tot, out=np.zeros_like(tot), where=tot > 0), 0)
        return np.divide(num, den, out=np.full(S.shape[1], np.nan), where=den > 0)
    rv = float(share_ratio(S0)[0]); Av = share_ratio(T)
    out["abcd"].append({"split": lab, "share_signal_like": r3(ns.sum() / y_raw.sum()), "odds_ratio": r3(rv),
                        "p_more_signal_like": r3(pk(Av, rv)), "p_less_signal_like": r3(pk(-Av, -rv)), "trials_median": r3(np.nanmedian(Av))})

(ROOT / "data/nuclear_basic.json").write_text(json.dumps(out, allow_nan=False))
print(json.dumps({k: v for k, v in out.items() if k not in ("scan",)} | {"scan": {k: v for k, v in out["scan"].items() if k != "grid"}},
                 default=str)[:5000])
