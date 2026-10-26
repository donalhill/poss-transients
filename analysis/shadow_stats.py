"""The shadow test's statistics, for one set of shadowed plates (real, replicated, or with permuted times).
Each plate enters as its transients' shares in shadow s (the share of the exposure their location spent in Earth's
shadow), the map's density a at them (the hosts' map A, mean 1 over the plate: where transients lie on the plates the
shadow never reaches), their weights u (ratings), and its shadow on the grid's sub-points. With F the share of a
plate's transients expected in shadow with no effect (s averaged over the plate, weighted by A):
  O_j = sum_i u_i s_i,  E_j = (sum_i u_i) F_j          observed and expected in shadow, the plate's own total held fixed
  Z = sum_j (E_j - O_j) / sqrt(sum_j (sum_i u_i^2)_j (E_A[s^2]_j - F_j^2))
Z is the score test of "fewer transients where s is higher", conditional on each plate's total (Berman's Z1 summed over
plates; with a sharp shadow, a Mantel-Haenszel rate ratio): positive for a shortfall in shadow. It is quoted as the
ratio O / E = sum O / sum E. Also, as checks: the ratio fit (density proportional to A (1 - (1 - r) s), the plate's
total conditioned on; r the density in shadow over that in sunlight), the glint-share fit (glint_fit.py: density
A - g s), and two pooled counts that do not hold each plate's total fixed: E1, VASCO's (the share of all transients in
shadow over the share of all plates' area), and E2, the count in shadow over the pointing model's expectation. The
hybrid H (B. Villarroel's suggestion): the ratio fit's likelihood on the plates the test uses, plus every crossed plate's
total, a negative binomial about the pointing model's expectation reduced in shadow, mu (1 - (1 - r) F), with the
model's dispersion alpha; the plates wholly or barely in shadow enter through their totals."""
import numpy as np
from numba import njit

import glint_fit as GF

NB = 1000


def geometry(sub, A):
    """A shadow's geometry: sub (cells x sub-points) the share in shadow at the grid's sub-points, A the map on the
    cells. S (the share of the area), F and E_A[s^2] (over the map), and glint_fit's QCAP and I."""
    w = A / A.sum(); sc = sub.mean(1); q = sub / A[:, None]
    return dict(S=float(sub.mean()), F=float(sc @ w), Es2=float((sub ** 2).mean(1) @ w), QCAP=float(q.max()),
                I=float(np.mean(q * sub) - sub.mean() ** 2))


@njit(cache=True)
def _llr(lr, S, W, F, C):
    r = np.exp(lr); L = 0.0
    for i in range(len(S)):
        L += W[i] * np.log(1.0 - (1.0 - r) * S[i])
    for j in range(len(F)):
        L -= C[j] * np.log(1.0 - (1.0 - r) * F[j])
    return L


@njit(cache=True)
def _fit_r(S, W, F, C, lo, hi):
    gr = 0.6180339887498949; a, b = lo, hi
    x1, x2 = b - gr * (b - a), a + gr * (b - a); f1, f2 = _llr(x1, S, W, F, C), _llr(x2, S, W, F, C)
    while b - a > 1e-8:
        if f1 < f2:
            a, x1, f1 = x1, x2, f2; x2 = a + gr * (b - a); f2 = _llr(x2, S, W, F, C)
        else:
            b, x2, f2 = x2, x1, f1; x1 = b - gr * (b - a); f1 = _llr(x1, S, W, F, C)
    return 0.5 * (a + b)


def ratio_inputs(s, u):
    """The ratio fit's inputs for one plate: shares s > 0 binned (NB bins on [0, 1], centres), weights k u summed, and
    k n; k = sum u / sum u^2 scales the weighted likelihood to the plate's effective number of transients."""
    k = float(u.sum() / (u ** 2).sum()) if len(u) else 1.0
    m = s > 0; b = np.minimum((s[m] * NB).astype(int), NB - 1)
    acc = np.bincount(b, weights=k * u[m], minlength=NB); nz = acc > 0
    return (np.flatnonzero(nz) + 0.5) / NB, acc[nz], k * float(u.sum())


def fit_r(P, interval=False):
    """The ratio fit over plates P (each with "s", "u", "F"); with interval, its 68% likelihood range (counting only)."""
    R = [ratio_inputs(p["s"], p["u"]) for p in P]
    S, W = np.concatenate([x[0] for x in R]), np.concatenate([x[1] for x in R])
    F, C = np.array([p["F"] for p in P]), np.array([x[2] for x in R])
    lo, hi = np.log(0.02), np.log(50.0)
    x = _fit_r(S, W, F, C, lo, hi)
    if not interval:
        return float(np.exp(x))
    Lm = _llr(x, S, W, F, C)

    def edge(end):
        if 2 * (Lm - _llr(end, S, W, F, C)) <= 1:
            return end
        a, b = x, end
        for _ in range(60):
            m = 0.5 * (a + b)
            a, b = (m, b) if 2 * (Lm - _llr(m, S, W, F, C)) <= 1 else (a, m)
        return 0.5 * (a + b)
    return float(np.exp(x)), float(np.exp(edge(lo))), float(np.exp(edge(hi)))


def _ll_tot(lr, n, mu, F, alpha):
    """The plates' totals: negative binomial log-likelihood (the terms depending on r) at ratio exp(lr)."""
    r = np.exp(lr); m = mu * (1 - (1 - r) * F); ia = 1 / alpha
    return float(np.sum(n * np.log(m / (m + ia)) + ia * np.log(ia / (m + ia))))


def fit_hybrid(P, Pall, alpha):
    """The hybrid's ratio: the ratio fit's likelihood on plates P plus the totals of every crossed plate in Pall (each
    with "u", "mu", "F")."""
    R = [ratio_inputs(p["s"], p["u"]) for p in P]
    S, W = np.concatenate([x[0] for x in R]), np.concatenate([x[1] for x in R])
    F, C = np.array([p["F"] for p in P]), np.array([x[2] for x in R])
    n = np.array([float(p["u"].sum()) for p in Pall]); mu = np.array([p["mu"] for p in Pall]); Ft = np.array([p["F"] for p in Pall])
    from scipy.optimize import minimize_scalar
    f = lambda lr: -(_llr(lr, S, W, F, C) + _ll_tot(lr, n, mu, Ft, alpha))
    return float(np.exp(minimize_scalar(f, bounds=(np.log(0.02), np.log(50.0)), method="bounded", options={"xatol": 1e-6}).x))


def fit_g(P, qmax):
    """The glint-share fit over plates P (each with "s", "a", "u", and the geometry)."""
    G = []
    for p in P:
        k = float(p["u"].sum() / (p["u"] ** 2).sum()); q = p["s"] / p["a"]
        Q, W = GF.binned(q, k * p["u"], qmax)
        G.append(dict(Q=Q, W=W, W0=float(k * p["u"][q <= 0].sum()), K=k, S=p["S"], QCAP=p["QCAP"]))
    return GF.fit(G)


def score(P):
    """O, E, Var and Z over plates P (each with "s", "u", "F", "Es2")."""
    O = sum(float(p["u"] @ p["s"]) for p in P); E = sum(float(p["u"].sum()) * p["F"] for p in P)
    V = sum(float((p["u"] ** 2).sum()) * (p["Es2"] - p["F"] ** 2) for p in P)
    return O, E, (E - O) / np.sqrt(V) if V > 0 else np.nan


def stats(P, totals=None, fits=False, qmax=None):
    """The statistics for plates P. totals = (N, J, mu): every plate's summed weight and number of plates, and each
    plate's pointing-model expectation (p["mu"]) for the pooled counts."""
    if not P:
        return {}
    O, E, Z = score(P)
    out = {"Z": Z, "OE": O / E if E > 0 else np.nan, "O": O, "E": E}
    if totals is not None:
        N, J = totals
        out["E1"] = (O / N) / (sum(p["S"] for p in P) / J)
        out["E2"] = O / sum(p["mu"] * p["F"] for p in P)
    if fits:
        out["r"] = fit_r(P); out["g"] = fit_g(P, qmax)
    return out


def pointing(pl, cat, k_fold=10):
    """Every plate's expected summed weight from the pointing model (expected.py), cross-fitted over ten folds by plate
    number (no plate sets its own), the plates' summed weights, and the model's dispersion: (mu, y, alpha)."""
    import pandas as pd
    import expected as EX
    y = cat.groupby("plate").w.sum().reindex(pl.index).fillna(0.0)
    P = pl.copy(); X, off = EX.design(P)
    fold = np.array([int(q[2:]) % k_fold for q in P.index]); mu = np.empty(len(P))
    for k in range(k_fold):
        Ek, _ = EX.fit(y.values, X, off, fold != k); mu[fold == k] = Ek[fold == k]
    _, alpha = EX.fit(y.values, X, off, np.ones(len(P), bool))
    return pd.Series(mu, index=pl.index), y, alpha
