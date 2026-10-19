"""Nuclear tests, measured photograph by photograph, as set out in NUCLEAR_DESIGN.md.

Each plate's transients are compared with those of plates from the same observing run (a free level for each run),
with every plate property in the model (plate_properties.py, less the time trend, which the runs replace). r is the
ratio for plates in VASCO's window (a US test's UT date, the day before or the day after) against the other plates of
their runs; negative binomial counts, the scatter between plates re-estimated in every fit.

Nulls, each fitted exactly like the data: (a) within each run, its exposed nights reassigned at random among its
nights; (b) the test calendar slid by whole years and up to 20 days. A Neyman interval on r from effects injected into
(a). Look-elsewhere over every window of days within 3, 7 and 14 days of a test, by (b). Robustness: each run and test
series left out, other run definitions, a clump-robust count, the central 5 x 5 deg of each plate.

python analysis/nuclear_fit.py -> data/nuclear_fit.json"""
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.stats import gamma, norm

import plate_properties as PP
import sample

ROOT = Path(__file__).resolve().parents[1]
DAY0 = pd.Timestamp("1949-01-01")
GAP = 4                                                             # a new run after more than 4 nights with no plates
NREP, NINJ = 1000, 300
SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-20, 21)]
RTRUE = np.round(np.arange(0.5, 3.01, 0.1), 2)
VASCO = 1.69                                                        # their excess, on archive dates
R = 14
HOURS = {"0 to 24 hours after a shot": (0, 24), "24 to 48 hours after": (24, 48), "24 hours before a shot": (-24, 0)}
r3 = lambda x: round(float(x), 3)
rng = np.random.default_rng(11)

nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
TESTS = np.unique(((pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date) - DAY0).dt.days).values)
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
shots["h"] = (pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601") - DAY0).dt.total_seconds() / 3600
shots["day"] = np.floor(shots.h / 24).astype(int)
shots["series"] = shots.source_url.str.extract(r"Operation_([A-Za-z]+)")[0].fillna(shots.Shot_Name)

cat, plates = sample.ours()
P = PP.properties(cat, plates)
day = ((P.date - DAY0).dt.days).values
hour = (P.ut - DAY0).dt.total_seconds().values / 3600
y_raw = P.n.values.astype(float)
y_cells = cat.assign(i=np.floor(cat.xi * 30).astype(int), j=np.floor(cat.eta * 30).astype(int)).groupby("plate") \
             .apply(lambda g: len(set(zip(g.i, g.j)))).reindex(P.index).fillna(0).values.astype(float)
inner = cat[(cat.xi.abs() < 2.5) & (cat.eta.abs() < 2.5)]
y_inner = inner.groupby("plate").size().reindex(P.index).fillna(0).values.astype(float)


def runs(gap):
    nights = np.sort(np.unique(day))
    return pd.Series(np.cumsum(np.r_[1, np.diff(nights) > gap]), index=nights).reindex(day).values


def window(days_=(-1, 0, 1), shift=0, tests=TESTS):
    return np.isin(day, (tests[:, None] + shift + np.array(days_)[None, :]).ravel()).astype(float)


def hours(lo, hi, shift_days=0):
    d = hour[:, None] - (shots.h.values[None, :] + 24.0 * shift_days)
    return ((d >= lo) & (d < hi)).any(1).astype(float)


def design(run):
    Z = PP.standardised(P).assign(run=run.astype(str))
    return np.asarray(patsy.dmatrix(PP.FORMULA.split("~")[1].replace(" + bs(t, df=16)", "") + " + C(run)", Z))


class Model:
    """Plate counts against the plate properties and a level for each run; r for a window S."""
    def __init__(self, y, run, keep=None):
        self.k = np.ones(len(y), bool) if keep is None else keep
        self.y, self.X = y[self.k], design(run)[self.k]
        self.X = self.X[:, np.abs(self.X).sum(0) > 0]                 # a run left out leaves an empty column
        self.alpha = 0.5

    def fit(self, S, alpha=None, iters=3):
        """r and its counting error; the scatter (alpha) re-estimated from the fit's residuals unless given."""
        S = S[self.k]; X = np.c_[self.X, S]; a = self.alpha if alpha is None else alpha
        for _ in range(1 if alpha is not None else iters):
            m = sm.GLM(self.y, X, family=sm.families.NegativeBinomial(alpha=a)).fit()
            mu = m.fittedvalues
            if alpha is None:
                a = max(1e-4, float(np.sum(((self.y - mu) ** 2 - mu) / mu ** 2) / (len(self.y) - X.shape[1])))
        return float(np.exp(m.params[-1])), float(m.bse[-1]), a, np.exp(np.log(mu) - m.params[-1] * S)


run0 = runs(GAP)
S0 = window()
M = Model(y_raw, run0)
r_hat, se, alpha, mu_bkg = M.fit(S0)
M.alpha = alpha
out = {"plates_near": int(S0.sum()), "nights_near": int(len(np.unique(day[S0 == 1]))), "runs_near": int(len(np.unique(run0[S0 == 1]))),
       "fit": {"r": r3(r_hat), "count_lo": r3(r_hat * np.exp(-se)), "count_hi": r3(r_hat * np.exp(se)), "alpha": r3(alpha)}}

# the identification table: every run with a plate in the window
ident = []
for r_ in np.unique(run0[S0 == 1]):
    k = run0 == r_; e = k & (S0 == 1); u_ = k & (S0 == 0)
    ident.append({"from": str(P.date[k].min().date()), "to": str(P.date[k].max().date()), "nights": int(len(np.unique(day[k]))),
                  "exposed": int(e.sum()), "unexposed": int(u_.sum()), "per_plate_exposed": round(float(y_raw[e].mean()), 1),
                  "per_plate_unexposed": round(float(y_raw[u_].mean()), 1) if u_.any() else None, "informative": bool(u_.any())})
out["runs"] = ident


# null (a): within each run, its exposed nights reassigned at random among its nights
def reassign(S, run):
    T = np.zeros_like(S)
    for r_ in np.unique(run[S == 1]):
        nights = np.unique(day[run == r_]); k = len(np.unique(day[(run == r_) & (S == 1)]))
        T[np.isin(day, rng.choice(nights, k, replace=False)) & (run == r_)] = 1
    return T


# null (b): the calendar slid by whole years and up to 20 days
SLIDES = [s for s in SHIFTS if window(shift=s).sum() >= 5 and (window(shift=s) * S0).sum() < 0.2 * window(shift=s).sum()]
pk = lambda d, x: float((np.sum(np.asarray(d) >= x) + 1) / (len(d) + 1))
q5 = lambda d: {k: r3(v) for k, v in zip(("q05", "q16", "median", "q84", "q95"), np.quantile(d, [0.05, 0.16, 0.5, 0.84, 0.95]))}


def calibrate(model, S, run, n_rep=NREP):
    r = model.fit(S)[0]
    A = np.array([model.fit(reassign(S, run))[0] for _ in range(n_rep)])
    B = np.array([model.fit(window(shift=s))[0] for s in SLIDES])
    return r, A, B


r_hat, A, B = calibrate(M, S0, run0)
out["null_runs"] = {**q5(A), "n": len(A), "p_as_high": r3(pk(A, r_hat)), "hist": np.histogram(A, np.linspace(0, 3, 61))[0].tolist()}
out["null_calendar"] = {**q5(B), "n": len(B), "p_as_high": r3(pk(B, r_hat)), "hist": np.histogram(B, np.linspace(0, 3, 61))[0].tolist()}
out["fit"]["sigma"] = r3(norm.isf(pk(A, r_hat)))


# injected effects, on null (a): r_true < 1 by removing transients on the reassigned plates, > 1 by removing them
# from the others
def injected(rt, n_rep):
    vals = []
    for _ in range(n_rep):
        S = reassign(S0, run0); c = y_raw.copy()
        if rt < 1:
            c[S == 1] = rng.binomial(c[S == 1].astype(int), rt)
        elif rt > 1:
            c[S == 0] = rng.binomial(c[S == 0].astype(int), 1 / rt)
        Mi = copy.copy(M); Mi.y = c
        vals.append(Mi.fit(S, alpha=alpha)[0])
    return np.array(vals)


belt = np.array([[rt, *np.quantile(injected(rt, NINJ), [0.05, 0.16, 0.5, 0.84, 0.95])] for rt in RTRUE])
i68 = belt[(belt[:, 2] <= r_hat) & (belt[:, 4] >= r_hat), 0]; i90 = belt[(belt[:, 1] <= r_hat) & (belt[:, 5] >= r_hat), 0]
out["interval68"] = [r3(i68.min()), r3(i68.max())] if len(i68) else None
out["interval90"] = [r3(i90.min()), r3(i90.max())] if len(i90) else None
out["power_vasco"] = {"r_true": VASCO, "p_as_low": r3(np.mean(injected(VASCO, NREP) <= r_hat))}

# robustness: each informative run left out; each test series left out; run definitions; outcomes; the inner plate
loo_run = []
for r_ in [r_ for r_ in np.unique(run0[S0 == 1]) if ((run0 == r_) & (S0 == 0)).any()]:
    keep = run0 != r_
    loo_run.append({"run": str(P.date[run0 == r_].min().date()), "r": r3(Model(y_raw, run0, keep).fit(S0)[0])})
loo_series = []
for sname, g in shots.groupby("series"):
    S_ = window(tests=np.setdiff1d(TESTS, g.day.values))
    if S_.sum() < S0.sum():
        loo_series.append({"series": sname, "plates_dropped": int(S0.sum() - S_.sum()), "r": r3(M.fit(S_)[0])})
gaps = []
for g_ in (3, 4, 6, 10):
    run_ = runs(g_); Mg = Model(y_raw, run_)
    rg, Ag, _ = calibrate(Mg, S0, run_, n_rep=400)
    gaps.append({"gap": g_, "runs": int(len(np.unique(run_))), "r": r3(rg), "p_runs": r3(pk(Ag, rg))})
outcomes = []
for lab, yv in (("transients", y_raw), ("2' cells holding a transient", y_cells), ("central 5 x 5 deg only", y_inner)):
    Mo = copy.copy(M); Mo.y = yv; Mo.alpha = Mo.fit(S0)[2]; ro, Ao, Bo = calibrate(Mo, S0, run0, n_rep=400)
    outcomes.append({"outcome": lab, "r": r3(ro), "p_runs": r3(pk(Ao, ro)), "p_calendar": r3(pk(Bo, ro))})
out["robust"] = {"leave_run_out": loo_run, "leave_series_out": loo_series, "run_gap": gaps, "outcomes": outcomes}

# secondary: windows in hours from each shot
out["physical"] = []
for lab, (lo, hi) in HOURS.items():
    S_ = hours(lo, hi)
    rs = M.fit(S_)[0]; Bs = [M.fit(hours(lo, hi, s))[0] for s in SLIDES]; As = [M.fit(reassign(S_, run0))[0] for _ in range(400)]
    out["physical"].append({"window": lab, "plates": int(S_.sum()), "r": r3(rs), "p_runs": r3(pk(As, rs)), "p_calendar": r3(pk(Bs, rs))})

# look-elsewhere: every window of days [a, b] within R days of a test (at least 5 plates), on the calendar slides.
# Hundreds of windows on hundreds of slides need a fast fit: the conditional Poisson form of the same comparison
# (Armstrong, Gasparrini & Tobias 2014), in which each run's level drops out exactly and r is a one-number fit; the
# plate properties enter through each plate's expected count from the model above. Data and slides alike.
o = np.exp(np.log(mu_bkg) - np.log(mu_bkg).mean())                  # each plate's expected count, up to its run's level
G = (run0[None, :] == np.unique(run0)[:, None]).astype(float)        # runs x plates
Nrun = G @ y_raw


def cond_r(S, iters=30):
    """r for each column of S, comparing plates within runs (conditional Poisson)."""
    o1, o0, N1 = G @ (o[:, None] * S), G @ (o[:, None] * (1 - S)), y_raw @ S
    g = np.zeros(S.shape[1])
    for _ in range(iters):
        e1 = o1 * np.exp(g); den = o0 + e1
        f = np.divide(e1, den, out=np.zeros_like(den), where=den > 0)
        score = N1 - (Nrun[:, None] * f).sum(0); info = (Nrun[:, None] * f * (1 - f)).sum(0)
        g = g + np.divide(score, info, out=np.zeros_like(score), where=info > 0)
    return np.exp(g)


WINDOWS = [(a, b) for a in range(-R, R + 1) for b in range(a, R + 1)]
WIN = [w for w in WINDOWS if window(tuple(range(w[0], w[1] + 1))).sum() >= 5]
cols = lambda s: np.stack([window(tuple(range(a, b + 1)), shift=s) for a, b in WIN], 1)
real_w = cond_r(cols(0))
null_w = np.array([cond_r(cols(s)) for s in SLIDES])
loc = np.array([pk(null_w[:, i], real_w[i]) for i in range(len(WIN))])
rank = (len(SLIDES) - np.argsort(np.argsort(null_w, axis=0), axis=0)) / (len(SLIDES) + 1)
span = np.array([max(abs(a), abs(b)) for a, b in WIN]); iP = WIN.index((-1, 1))
out["scan"] = {"windows": len(WIN), "primary_r": r3(real_w[iP]), "primary_p_local": r3(loc[iP]),
               "primary_rank": int(np.sum(loc <= loc[iP])),
               "share_trials_rank_as_high": r3(np.mean([np.sum(rank[t] <= rank[t, iP]) <= np.sum(loc <= loc[iP]) for t in range(len(SLIDES))])),
               "grid": [{"a": a, "b": b, "r": r3(real_w[i]), "p": r3(loc[i]), "z": round(float(norm.isf(loc[i])), 2)} for i, (a, b) in enumerate(WIN)],
               "by_size": []}
for days_ in (3, 7, 14):
    cols = np.flatnonzero(span <= days_); b_ = cols[np.argmin(loc[cols])]
    pg = (np.sum(rank[:, cols].min(1) <= loc[b_]) + 1) / (len(SLIDES) + 1)
    out["scan"]["by_size"].append({"days": days_, "windows": len(cols), "best": list(WIN[b_]), "best_r": r3(real_w[b_]),
                                   "best_p_local": r3(loc[b_]), "p_global": r3(pg), "sigma": r3(norm.isf(pg))})

# each night in the window on its own: its plates' transients over what the model expects of them with no effect
out["nights"] = []
for dd in np.unique(day[S0 == 1]):
    k = day == dd; N_, E_ = y_raw[k].sum(), mu_bkg[k].sum()
    out["nights"].append({"date": str((DAY0 + pd.Timedelta(days=int(dd))).date()), "plates": list(P.index[k]), "n": int(N_), "E": round(float(E_), 1),
                          "r": r3(N_ / E_), "lo": r3(gamma.ppf(0.16, N_) / E_ if N_ > 0 else 0.0), "hi": r3(gamma.ppf(0.84, N_ + 1) / E_)})
(ROOT / "data/nuclear_fit.json").write_text(json.dumps(out))
print(json.dumps({k: v for k, v in out.items() if k not in ("nights",)}, default=str)[:6000])
