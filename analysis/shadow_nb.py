"""Between-plate test of avoidance with an overdispersed model: N_p ~ NegBin(mu_p, alpha),
log mu_p = b0 + b1 x + b2 x^2 + log(1 - phi f_p), x = |galactic latitude|/30 (strong plate-level covariate).
Profile likelihood for phi (the fraction of transients that would vanish inside the shadow). Corrected dates."""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import minimize
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parents[1]
geom = os.environ.get("LL_GEOM", "geo")
pt = pd.read_parquet(ROOT / f"data/plate_shadow_{geom}_corr.parquet")
pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet").loc[pt.index]
x = np.abs(SkyCoord(pl.ra_c.values * u.deg, pl.dec_c.values * u.deg).galactic.b.deg) / 30
N, f = pt.N.values.astype(float), pt.f_tmpl.values


def nll(theta, phi):
    b0, b1, b2, la = theta
    mu = np.exp(b0 + b1 * x + b2 * x**2) * (1 - phi * f)
    a = np.exp(la); r = 1 / a
    return -np.sum(gammaln(N + r) - gammaln(r) - gammaln(N + 1) + r * np.log(r / (r + mu)) + N * np.log(mu / (r + mu)))


def prof(phi):
    return minimize(nll, [np.log(N.mean()), 0, 0, 0], args=(phi,), method="Nelder-Mead",
                    options=dict(maxiter=20000, xatol=1e-6, fatol=1e-8)).fun


grid = np.linspace(-1.5, 0.99, 250)
L = np.array([prof(p) for p in grid])
best = grid[L.argmin()]
ci = grid[L - L.min() < 1.92]
print(f"[{geom}] NB between-plate: phi_hat = {best:+.2f}, 95% profile CI {ci.min():+.2f} .. {ci.max():+.2f}; "
      f"-2dlnL(phi=0) = {2 * (prof(0) - L.min()):.2f}; NB alpha = {np.exp(minimize(nll, [np.log(N.mean()), 0, 0, 0], args=(0,), method='Nelder-Mead').x[3]):.2f}")
