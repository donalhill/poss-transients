"""The shadow test's model and likelihood: transients as artefacts plus glints, fitted transient by transient.
On a plate, the density of transients at a location x, relative to the plate's mean, is
  A(x) - g s(x)
A is the average density over the face of the plates no shadow reaches (the hosts' map, mean 1 over the plate): every
location there is sunlit, so A holds the artefacts, which gather towards the edges, and any glints in sunlight. Glints
are spread evenly over a plate, a share g of its transients, and vanish in Earth's shadow, so a location in shadow for a
share s of the exposure loses g s of its density; the artefacts are unaffected. g cannot exceed the smallest value of A
under the shadow, where glints would be every transient. Given the plate's number of transients n, the likelihood of
where they lie is
  ln L(g) = sum_i ln(1 - g q_i) - n ln(1 - g S),    q = s / A at each transient, S the plate's mean share in shadow
so how rich or poor a plate is does not enter. One g is shared by the plates the shadow partly covers.
Transients count by their weight u (a rating), with every plate's terms scaled by k = sum u / sum u^2 so the
likelihood carries the plate's effective number of transients. The shares q are binned (NB bins on [0, QMAX]) so a
fit takes well under a millisecond: the bin width moves ln(1 - g q) by under 1e-3.
Plate inputs (inputs()): q binned (Q, W: bin centres and summed k u), W0 (k u summed over the transients never in
shadow), K = k, S, QCAP (the largest s / A on the plate's grid, which sets g's ceiling), I (the information on g per
transient at g = 0), and nw (the summed weight)."""
import numpy as np
from numba import njit

NB = 1000
G_LO = -20.0                                                                 # g's lower search bound (an excess in shadow)


def binned(q, w, qmax):
    """Shares q > 0 and weights w summed into NB bins on [0, qmax] (centres)."""
    m = q > 0; b = np.minimum((q[m] / qmax * NB).astype(int), NB - 1)
    acc = np.bincount(b, weights=w[m], minlength=NB); nz = acc > 0
    return (np.flatnonzero(nz) + 0.5) / NB * qmax, acc[nz]


def inputs(s, a, u, sub_s, A_cells, qmax):
    """One plate's inputs: s, a, u at its transients (share in shadow, the map's density there, weight); sub_s, the share
    in shadow at the grid's sub-points (cells x sub-points) and A_cells, the map on the cells."""
    k = float(u.sum() / (u ** 2).sum()) if len(u) else 1.0
    q = s / a; Q, W = binned(q, k * u, qmax)
    qs = sub_s / A_cells[:, None]
    return dict(Q=Q, W=W, W0=float(k * u[q <= 0].sum()), K=k, S=float(sub_s.mean()), QCAP=float(qs.max()),
                I=float(np.mean(qs * sub_s) - sub_s.mean() ** 2), nw=float(u.sum()))


def pack(P):
    """Plates' inputs as flat arrays for the numba functions."""
    off = np.zeros(len(P) + 1, np.int64); off[1:] = np.cumsum([len(p["Q"]) for p in P])
    cat = lambda key: np.concatenate([p[key] for p in P]) if P else np.zeros(0)
    arr = lambda key: np.array([p[key] for p in P], float)
    return cat("Q"), cat("W"), off, arr("W0"), arr("S")


@njit(cache=True)
def ll(g, Q, W, off, W0, S):
    """The log-likelihood of the shared g, summed over the plates."""
    L = 0.0
    for j in range(len(S)):
        n = W0[j]
        for i in range(off[j], off[j + 1]):
            L += W[i] * np.log(1.0 - g * Q[i]); n += W[i]
        L -= n * np.log(1.0 - g * S[j])
    return L


@njit(cache=True)
def _golden(lo, hi, Q, W, off, W0, S):
    """The maximum of the (concave) log-likelihood on [lo, hi], by golden-section search."""
    gr = 0.6180339887498949
    a, b = lo, hi
    x1, x2 = b - gr * (b - a), a + gr * (b - a)
    f1, f2 = ll(x1, Q, W, off, W0, S), ll(x2, Q, W, off, W0, S)
    while b - a > 1e-7:
        if f1 < f2:
            a, x1, f1 = x1, x2, f2; x2 = a + gr * (b - a); f2 = ll(x2, Q, W, off, W0, S)
        else:
            b, x2, f2 = x2, x1, f1; x1 = b - gr * (b - a); f1 = ll(x1, Q, W, off, W0, S)
    return 0.5 * (a + b)


def ceiling(P):
    """g's ceiling over the plates P: the density stays positive everywhere under the shadow."""
    return 1.0 / max(p["QCAP"] for p in P) if P else 1.0


def fit(P):
    """The shared g of highest likelihood over the plates P (inputs() records); NaN with no plates."""
    if not P:
        return np.nan
    return float(_golden(G_LO, 0.999 * ceiling(P), *pack(P)))


def loglik(P, g):
    return float(ll(g, *pack(P)))


def interval(P, level=1.0):
    """The estimate and its profile-likelihood range (2 dlnL = level: 1 for 68%), each end found by bisection."""
    a = pack(P); x = fit(P); Lm = ll(x, *a)

    def edge(end):
        if 2 * (Lm - ll(end, *a)) <= level:
            return end
        u, v = x, end
        for _ in range(80):
            m = 0.5 * (u + v)
            u, v = (m, v) if 2 * (Lm - ll(m, *a)) <= level else (u, m)
        return 0.5 * (u + v)
    return x, edge(G_LO), edge(0.999 * ceiling(P))
