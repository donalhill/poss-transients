"""Bounded interval for phi, the fraction of transients that vanish inside Earth's shadow (0 <= phi <= 1).

Sampling model: the real straddling shadow footprints pasted onto random plates (shadow_paste_null.py), with a true
phi imposed by removing each in-footprint transient with probability phi. For every phi on a grid this gives the
distribution of the within-plate estimator phi_hat. From these we build
  - a Feldman-Cousins (unified) confidence interval: Neyman construction, likelihood-ratio ordering, phi in [0, 1];
  - a Bayesian credible interval with a flat prior on [0, 1].
Writes data/phi_interval_{GEOM}.json."""
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).parent))
from shadow import ROOT

GEOM = os.environ.get("LL_GEOM", "geo")
rng = np.random.default_rng(11)
d = np.load(ROOT / f"data/paste_phi_{GEOM}.npz")
O, N, f = d["O"], d["N"], d["f"]
keep = np.sum((N - O) * f, axis=1) > 0
O, N, f = O[keep], N[keep], f[keep]
den = np.sum((N - O) * f, axis=1)


def phi_hat(Oi):
    return 1 - np.sum(Oi * (1 - f), axis=1) / den


x0 = float(1 - np.sum(d["obs_O"] * (1 - d["obs_f"])) / np.sum((d["obs_N"] - d["obs_O"]) * d["obs_f"]))
phis = np.round(np.arange(0, 1.0001, 0.0025), 4)
edges = np.linspace(-2.5, 1.0, 701)
xc, dx = 0.5 * (edges[1:] + edges[:-1]), edges[1] - edges[0]
P = np.empty((len(phis), len(xc)))      # P[i, j] = probability that phi_hat falls in bin j when the true value is phis[i]
# common random numbers: the same uniforms decide which transients survive at every phi, so curves in phi are smooth
U = [rng.random(O.shape) for _ in range(4)]
for i, ph in enumerate(phis):
    xs = np.concatenate([phi_hat(binom.ppf(u, O.astype(int), 1 - ph)) for u in U])
    h, _ = np.histogram(np.clip(xs, edges[0], edges[-1] - 1e-9), edges)
    h = gaussian_filter1d(h.astype(float), 3.0, mode="constant") if ph < 1 else h.astype(float)
    P[i] = h / h.sum()

j0 = int(np.clip(np.searchsorted(edges, x0) - 1, 0, len(xc) - 1))
best = P.max(axis=0)                     # likelihood at the best physical phi, for each possible measurement


def fc(cl):
    accept = np.zeros_like(P, dtype=bool)
    for i in range(len(phis)):
        r = P[i] / np.where(best > 0, best, 1)
        order = np.argsort(-r)
        n = np.searchsorted(np.cumsum(P[i][order]), cl) + 1
        accept[i, order[:n]] = True
    inside = phis[accept[:, j0]]
    step = phis[1] - phis[0]   # the edge lies between the last accepted grid point and the next one
    return float(inside.min()), float(min(inside.max() + step / 2, 1))


like = P[:, j0]
post = like / np.trapezoid(like, phis)
cdf = np.concatenate([[0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(phis))])
q = lambda p: float(np.interp(p, cdf, phis))
null = P[0]
res = {
    "phi_hat": round(x0, 3),
    "mode": float(phis[np.argmax(like)]),
    "null_sd": round(float(np.sqrt(np.sum(null * xc ** 2) - np.sum(null * xc) ** 2)), 3),
    "fc": {str(c): [round(v, 3) for v in fc(c)] for c in (0.683, 0.9, 0.95)},
    "bayes": {"upper90": round(q(0.9), 3), "upper95": round(q(0.95), 3), "median": round(q(0.5), 3)},
    "p_as_low": {str(ph): round(float(P[np.argmin(abs(phis - ph))][: j0 + 1].sum()), 3) for ph in (0.32, 0.55)},
    "posterior": {"phi": phis[::4].tolist(), "density": [round(v, 4) for v in post[::4]]},
}
print(json.dumps({k: v for k, v in res.items() if k != "posterior"}, indent=1))
(ROOT / f"data/phi_interval_{GEOM}.json").write_text(json.dumps(res))
