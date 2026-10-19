"""The data and model behind the nuclear-test analysis (NUCLEAR_DESIGN.md), shared by nuclear_fit.py and
nuclear_checks.py: the plates and their properties, observing runs, windows around tests, the full within-run fit on
the data, and the negative binomial fit vectorised over trials."""
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm

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


def full_fit(y, run, S, iters=4, keep=None):
    """The full fit on the data: plate properties, a level for each run, and r; the scatter alpha from the residuals.
    Returns r, its counting error, alpha, the model's degrees of freedom, and each plate's expected count with no effect
    (for the plates in keep, if given)."""
    k = np.ones(len(y), bool) if keep is None else keep
    Z = PP.standardised(P).assign(run=run.astype(str))[k]
    X = np.c_[np.asarray(patsy.dmatrix(PP.FORMULA.split("~")[1].replace(" + bs(t, df=16)", "") + " + C(run)", Z)), S[k]]
    X = X[:, np.abs(X).sum(0) > 0]; y, S = y[k], S[k]
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

    def r(self, S, Y, iters=30, rounds=2, levels=False):
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
        if levels:
            return np.exp(g), lam, a
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


