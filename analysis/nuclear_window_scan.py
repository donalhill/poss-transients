"""Look-elsewhere check for the nuclear-test excess: every window, with its start and end left free.

VASCO's window (one day either side of a US test) is one choice among many. Here every run of consecutive days [a, b]
within R days of a test is tried, for R = 3, 7 and 14 (28, 120 and 435 windows). The test is VASCO's random-night test
generalised to all windows at once: shuffling nights shuffles each night's membership of every window together. For
each window the local p is the share of shuffles that match or beat the real excess; the global p is the share of
shuffles whose best window (smallest local p) does at least as well as the real best window.

Run on transient counts as they are, and with the Milky Way effect removed (each night's transients measured against
the number expected from its plates' distances to the Milky Way, as in build_site.py).

Results (5,000 shuffles): with the Milky Way effect removed, the global significance is about 0.1, 0.9 and 0.7
standard deviations for R = 3, 7 and 14 (the exact values move a little with the random shuffles). Without it, the best windows are broad ones in the days after tests, and for R = 14
no shuffle matches the real best window (global p < 1 in 5,000): an excess that removing the Milky Way effect erases."""
import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
NSHUF = 5000
g = {"__file__": str(ROOT / "analysis/nuclear_lee.py"), "__name__": "scan"}
with contextlib.redirect_stdout(io.StringIO()):
    exec((ROOT / "analysis/nuclear_lee.py").read_text().split("windows = {")[0], g)
nights, values, N, nt = g["nights"], g["values"], g["N"], g["nt"]
cat, pl = g["g"]["cat"], g["g"]["pl"]
day0 = pd.Timestamp("1949-01-01")
us = np.unique(((pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date) - day0).dt.days).values)

# expected transients per night from the plates' distances to the Milky Way (same curve as build_site.py)
pln = pl[pl.cor_date.isin(nights)].copy()
pln["n"] = cat.groupby("plate").size().reindex(pln.index).fillna(0)
pln["b"] = np.abs(SkyCoord(pln.ra_c.values * u.deg, pln.dec_c.values * u.deg).galactic.b.deg)
sat = lambda x, A, b0, C: C + A * (1 - np.exp(-x / b0))
psat, _ = curve_fit(sat, pln.b, pln.n, p0=[200, 20, 30], maxfev=20000)
expected = pln.assign(e=sat(pln.b, *psat)).groupby("cor_date").e.sum().reindex(nights).fillna(0).values

rng = np.random.default_rng(2026)


def excess(M, e):
    a, b = values @ M, e @ M
    return (a / b) / ((values.sum() - a) / (e.sum() - b))


for R in (3, 7, 14):
    wins = [(a, b) for a in range(-R, R + 1) for b in range(a, R + 1)]
    M = np.array([np.isin(N, (us[:, None] + np.arange(a, b + 1)[None, :]).ravel()) for a, b in wins]).T.astype(float)
    for name, e in (("counts as they are", np.ones(len(values))), ("Milky Way effect removed", expected)):
        real = excess(M, e)
        P = np.array([excess(M[rng.permutation(len(values))], e) for _ in range(NSHUF)])
        local = np.array([np.mean(P[:, i] >= real[i]) for i in range(len(wins))])
        rank = (np.argsort(np.argsort(-P, axis=0), axis=0) + 1) / len(P)
        glob = np.mean(rank.min(axis=1) <= local.min())
        ib = int(local.argmin())
        g_sd = f"{norm.isf(glob):.2f} sd" if glob > 0 else f"> {norm.isf(1 / NSHUF):.1f} sd"
        print(f"+/-{R:2d} days ({len(wins):3d} windows), {name}: best window {wins[ib]} days, local p {local[ib]:.4f}; "
              f"VASCO's +/-1 local p {local[wins.index((-1, 1))]:.3f}; global p {glob:.3f} ({g_sd})")
