"""Nuclear tests, measured photograph by photograph, as set out in NUCLEAR_DESIGN.md.

Each plate's transients are compared with those of plates from the same observing run (a free level for each run),
with every plate property in the model (plate_properties.py, less the time trend, which the runs replace). r is the
ratio for plates in VASCO's window (a US test's UT date, the day before or the day after) against the other plates of
their runs; negative binomial counts.

The data are fitted once in full. Every trial then keeps the plate-property slopes of that fit and refits the run
levels, r and the scatter between plates: a negative binomial fit vectorised over trials, in which each step solves
for the run levels and r exactly (the equations couple each run only to r).

Nulls: (a) within each run, its exposed nights reassigned at random among its nights; (b) the test calendar slid by
whole years and up to 20 days. A Neyman interval on r from effects injected into (a). Look-elsewhere over every window
of days within 3, 7 and 14 days of a test, by (b). Robustness: each run and test series left out, other run
definitions, a clump-robust count, the central 5 x 5 deg of each plate.

python analysis/nuclear_fit.py -> data/nuclear_fit.json"""
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
NREP, NINJ = 2000, 500
SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-20, 21)]
RTRUE = np.round(np.arange(0.5, 3.01, 0.05), 2)
VASCO = 1.69                                                        # their excess, on archive dates
R = 14
HOURS = {"0 to 24 hours after a shot": (0, 24), "24 to 48 hours after": (24, 48), "24 hours before a shot": (-24, 0)}
r3 = lambda x: round(float(x), 3)

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
y_inner = cat[(cat.xi.abs() < 2.5) & (cat.eta.abs() < 2.5)].groupby("plate").size().reindex(P.index).fillna(0).values.astype(float)


def runs(gap):
    nights = np.sort(np.unique(day))
    return pd.Series(np.cumsum(np.r_[1, np.diff(nights) > gap]), index=nights).reindex(day).values


def window(days_=(-1, 0, 1), shift=0, tests=TESTS):
    return np.isin(day, (tests[:, None] + shift + np.array(days_)[None, :]).ravel()).astype(float)


def hours(lo, hi, shift_days=0):
    d = hour[:, None] - (shots.h.values[None, :] + 24.0 * shift_days)
    return ((d >= lo) & (d < hi)).any(1).astype(float)


def full_fit(y, run, S, iters=4):
    """The full fit on the data: plate properties, a level for each run, and r; the scatter alpha from the residuals.
    Returns r, its counting error, alpha, the model's degrees of freedom, and each plate's expected count with no effect."""
    Z = PP.standardised(P).assign(run=run.astype(str))
    X = np.c_[np.asarray(patsy.dmatrix(PP.FORMULA.split("~")[1].replace(" + bs(t, df=16)", "") + " + C(run)", Z)), S]
    a = 0.5
    for _ in range(iters):
        m = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=a)).fit()
        mu = m.fittedvalues
        a = max(1e-4, float(np.sum(((y - mu) ** 2 - mu) / mu ** 2) / (len(y) - X.shape[1])))
    return float(np.exp(m.params[-1])), float(m.bse[-1]), a, X.shape[1], np.exp(np.log(mu) - m.params[-1] * S)


class Trials:
    """r for many windows (columns of S) and counts (columns of Y) at once: counts ~ NB(o * exp(level[run] + g S), alpha),
    o each plate's expected count from the full fit, the run levels and g refitted by Newton steps (each run couples
    only to g, so each step is solved exactly), alpha re-estimated per column from the residuals and the fit repeated."""
    def __init__(self, o, run, alpha, dof, keep=None):
        k = np.ones(len(o), bool) if keep is None else keep
        self.k, self.o = k, o[k]
        u, self.ri = np.unique(run[k], return_inverse=True)
        self.G = (self.ri[None, :] == np.arange(len(u))[:, None]).astype(float)       # runs x plates
        self.alpha, self.dof = alpha, dof

    def r(self, S, Y, iters=30, rounds=2):
        S = np.atleast_2d(np.asarray(S, float).T).T[self.k]; Y = np.broadcast_to(np.atleast_2d(np.asarray(Y, float).T).T[self.k], S.shape)
        W = S.shape[1]; lam = np.zeros((self.G.shape[0], W)); g = np.zeros(W); a = np.full(W, self.alpha)
        logo = np.log(self.o)[:, None]
        for rnd in range(rounds):
            for _ in range(iters):
                mu = np.exp(logo + lam[self.ri] + g * S)
                u = (Y - mu) / (1 + a * mu); h = mu * (1 + a * Y) / (1 + a * mu) ** 2
                gl, Hll, Hlg = self.G @ u, self.G @ h, self.G @ (h * S)
                gg, Hgg = (u * S).sum(0), (h * S).sum(0)
                inv = np.divide(1.0, Hll, out=np.zeros_like(Hll), where=Hll > 0)
                den = Hgg - (Hlg ** 2 * inv).sum(0)
                dg = np.clip(np.divide(gg - (Hlg * gl * inv).sum(0), den, out=np.zeros_like(den), where=den > 1e-12), -1, 1)
                g = g + dg; lam = lam + (gl - Hlg * dg) * inv
            if rnd < rounds - 1:
                mu = np.exp(logo + lam[self.ri] + g * S)
                a = np.maximum(1e-4, (((Y - mu) ** 2 - mu) / mu ** 2).sum(0) / (len(self.o) - self.dof))
        ok = ((self.G @ S > 0) & (self.G @ (1 - S) > 0)).any(0)         # a run holding exposed and unexposed plates
        return np.where(ok, np.exp(g), np.nan)


fin = lambda d: np.asarray(d, float)[np.isfinite(np.asarray(d, float))]
pk = lambda d, x: float((np.sum(fin(d) >= x) + 1) / (len(fin(d)) + 1))
q5 = lambda d: {k: r3(v) for k, v in zip(("q05", "q16", "median", "q84", "q95"), np.quantile(fin(d), [0.05, 0.16, 0.5, 0.84, 0.95]))}


def reassign(S, run, n, seed):
    """n copies of S with each run's exposed nights reassigned at random among its nights (columns)."""
    g = np.random.default_rng(seed); T = np.zeros((len(S), n))
    for r_ in np.unique(run[S == 1]):
        inrun = run == r_; nights = np.unique(day[inrun]); k = len(np.unique(day[inrun & (S == 1)]))
        pick = np.argsort(g.random((n, len(nights))), 1)[:, :k]                      # k nights per copy
        chosen = nights[pick]                                                          # n x k
        T[inrun] = (day[inrun][:, None, None] == chosen[None, :, :]).any(2)
    return T


run0 = runs(GAP)
S0 = window()
r_hat, se, alpha, dof, mu_bkg = full_fit(y_raw, run0, S0)
TR = Trials(mu_bkg, run0, alpha, dof)
assert abs(TR.r(S0, y_raw)[0] - r_hat) < 0.01 * r_hat, (TR.r(S0, y_raw)[0], r_hat)     # the batched fit agrees with the full one
SLIDES = [s for s in SHIFTS if window(shift=s).sum() >= 5 and (window(shift=s) * S0).sum() < 0.2 * window(shift=s).sum()]
slides = lambda f: np.stack([f(s) for s in SLIDES], 1)


def calibrate(tr, S, y, run, seed, n=NREP, slide=lambda s: window(shift=s)):
    r = float(tr.r(S, y)[0])
    return r, tr.r(reassign(S, run, n, seed), y), tr.r(slides(slide), y)


r_hat, A, B = calibrate(TR, S0, y_raw, run0, 1)
out = {"plates_near": int(S0.sum()), "nights_near": int(len(np.unique(day[S0 == 1]))), "runs_near": int(len(np.unique(run0[S0 == 1]))),
       "fit": {"r": r3(r_hat), "count_lo": r3(r_hat * np.exp(-se)), "count_hi": r3(r_hat * np.exp(se)), "alpha": r3(alpha),
               "sigma": r3(norm.isf(pk(A, r_hat)))},
       "null_runs": {**q5(A), "n": len(fin(A)), "p_as_high": r3(pk(A, r_hat)), "hist": np.histogram(fin(A), np.linspace(0, 3, 61))[0].tolist()},
       "null_calendar": {**q5(B), "n": len(fin(B)), "p_as_high": r3(pk(B, r_hat)), "hist": np.histogram(fin(B), np.linspace(0, 3, 61))[0].tolist()}}

# the identification table: every run with a plate in the window
out["runs"] = []
for r_ in np.unique(run0[S0 == 1]):
    k = run0 == r_; e = k & (S0 == 1); u_ = k & (S0 == 0)
    out["runs"].append({"from": str(P.date[k].min().date()), "to": str(P.date[k].max().date()), "nights": int(len(np.unique(day[k]))),
                        "exposed": int(e.sum()), "unexposed": int(u_.sum()), "per_plate_exposed": round(float(y_raw[e].mean()), 1),
                        "per_plate_unexposed": round(float(y_raw[u_].mean()), 1) if u_.any() else None, "informative": bool(u_.any())})


# injected effects, on null (a): r_true < 1 by removing transients on the reassigned plates, > 1 by removing them
# from the others
def injected(rt, n, seed):
    S = reassign(S0, run0, n, seed); g = np.random.default_rng(seed + 1)
    Y = np.repeat(y_raw[:, None], n, 1).astype(int)
    if rt < 1:
        Y = np.where(S == 1, g.binomial(Y, rt), Y)
    elif rt > 1:
        Y = np.where(S == 0, g.binomial(Y, 1 / rt), Y)
    return TR.r(S, Y)


belt = np.array([[rt, *np.quantile(fin(injected(rt, NINJ, 100 + i)), [0.05, 0.16, 0.5, 0.84, 0.95])] for i, rt in enumerate(RTRUE)])
i68 = belt[(belt[:, 2] <= r_hat) & (belt[:, 4] >= r_hat), 0]; i90 = belt[(belt[:, 1] <= r_hat) & (belt[:, 5] >= r_hat), 0]
out["interval68"] = [r3(i68.min()), r3(i68.max())] if len(i68) else None
out["interval90"] = [r3(i90.min()), r3(i90.max())] if len(i90) else None
out["power_vasco"] = {"r_true": VASCO, "p_as_low": r3(np.mean(fin(injected(VASCO, NREP, 7)) <= r_hat))}

# robustness: each informative run left out; each test series left out; run definitions; outcomes
out["robust"] = {"leave_run_out": [], "leave_series_out": [], "run_gap": [], "outcomes": []}
for r_ in [r_ for r_ in np.unique(run0[S0 == 1]) if ((run0 == r_) & (S0 == 0)).any()]:
    tk = Trials(mu_bkg, run0, alpha, dof, keep=run0 != r_)
    out["robust"]["leave_run_out"].append({"run": str(P.date[run0 == r_].min().date()), "r": r3(tk.r(S0, y_raw)[0])})
for sname, g_ in shots.groupby("series"):
    S_ = window(tests=np.setdiff1d(TESTS, g_.day.values))
    if S_.sum() < S0.sum():
        out["robust"]["leave_series_out"].append({"series": sname, "plates_dropped": int(S0.sum() - S_.sum()), "r": r3(TR.r(S_, y_raw)[0])})
for gp in (3, 4, 6, 10):
    run_ = runs(gp); rg, _, ag, dg, mg = full_fit(y_raw, run_, S0)
    rg, Ag, Bg = calibrate(Trials(mg, run_, ag, dg), S0, y_raw, run_, 20 + gp)
    out["robust"]["run_gap"].append({"gap": gp, "runs": int(len(np.unique(run_))), "r": r3(rg), "p_runs": r3(pk(Ag, rg)), "p_calendar": r3(pk(Bg, rg))})
for i, (lab, yv) in enumerate((("transients", y_raw), ("2' cells holding a transient", y_cells), ("central 5 x 5 deg only", y_inner))):
    _, _, ao, do, mo = full_fit(yv, run0, S0)
    ro, Ao, Bo = calibrate(Trials(mo, run0, ao, do), S0, yv, run0, 40 + i)
    out["robust"]["outcomes"].append({"outcome": lab, "r": r3(ro), "p_runs": r3(pk(Ao, ro)), "p_calendar": r3(pk(Bo, ro))})

# secondary: windows in hours from each shot
out["physical"] = []
for i, (lab, (lo, hi)) in enumerate(HOURS.items()):
    rs, As, Bs = calibrate(TR, hours(lo, hi), y_raw, run0, 60 + i, slide=lambda s, lo=lo, hi=hi: hours(lo, hi, s))
    out["physical"].append({"window": lab, "plates": int(hours(lo, hi).sum()), "r": r3(rs), "p_runs": r3(pk(As, rs)), "p_calendar": r3(pk(Bs, rs))})

# look-elsewhere: every window of days [a, b] within R days of a test (at least 5 plates), on the calendar slides
WIN = [(a, b) for a in range(-R, R + 1) for b in range(a, R + 1) if window(tuple(range(a, b + 1))).sum() >= 5]
cols = lambda s: np.stack([window(tuple(range(a, b + 1)), shift=s) for a, b in WIN], 1)
real_w = TR.r(cols(0), y_raw)
null_w = np.array([TR.r(cols(s), y_raw) for s in SLIDES])
null_w = np.where(np.isfinite(null_w), null_w, -np.inf)                 # an undetermined window never beats the data
loc = np.array([pk(null_w[:, i], real_w[i]) for i in range(len(WIN))])
rank = (len(SLIDES) - np.argsort(np.argsort(null_w, axis=0), axis=0)) / (len(SLIDES) + 1)
span = np.array([max(abs(a), abs(b)) for a, b in WIN]); iP = WIN.index((-1, 1))
out["scan"] = {"windows": len(WIN), "primary_r": r3(real_w[iP]), "primary_p_local": r3(loc[iP]), "primary_rank": int(np.sum(loc <= loc[iP])),
               "share_trials_rank_as_high": r3(np.mean([np.sum(rank[t] <= rank[t, iP]) <= np.sum(loc <= loc[iP]) for t in range(len(SLIDES))])),
               "grid": [{"a": a, "b": b, "r": r3(real_w[i]), "p": r3(loc[i]), "z": round(float(norm.isf(loc[i])), 2)} for i, (a, b) in enumerate(WIN)],
               "by_size": []}
for days_ in (3, 7, 14):
    c_ = np.flatnonzero(span <= days_); b_ = c_[np.argmin(loc[c_])]
    pg = (np.sum(rank[:, c_].min(1) <= loc[b_]) + 1) / (len(SLIDES) + 1)
    out["scan"]["by_size"].append({"days": days_, "windows": len(c_), "best": list(WIN[b_]), "best_r": r3(real_w[b_]),
                                   "best_p_local": r3(loc[b_]), "p_global": r3(pg), "sigma": r3(norm.isf(pg))})

# each night in the window on its own: its plates' transients over what the model expects of them with no effect
out["nights"] = []
for dd in np.unique(day[S0 == 1]):
    k = day == dd; N_, E_ = y_raw[k].sum(), mu_bkg[k].sum()
    out["nights"].append({"date": str((DAY0 + pd.Timedelta(days=int(dd))).date()), "plates": list(P.index[k]), "n": int(N_), "E": round(float(E_), 1),
                          "r": r3(N_ / E_), "lo": r3(gamma.ppf(0.16, N_) / E_ if N_ > 0 else 0.0), "hi": r3(gamma.ppf(0.84, N_ + 1) / E_)})
(ROOT / "data/nuclear_fit.json").write_text(json.dumps(out))
print("written")
