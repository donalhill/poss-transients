"""Earth-shadow deficit: reproduce the published per-transient binomial, then redo it with plates as the unit.

Geometries: 'geo' = geocentric anti-sun (J2000 axes, plate UT), as in VASCO/Doherty up to precession;
            'topo' = line of sight from Palomar to the GEO sphere, angle from the shadow axis.
Expected fractions per plate from a dense grid over the plate footprint, weighted either uniformly
(area, as the published control sample) or by the stacked plate-frame density of transients.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).parent))
from shadow import ROOT, SHADOW_DEG, antisun_geocentric, sep_topocentric
import os

CORR = os.environ.get("LL_DATES") == "corrected"  # use plate dates corrected from labels/asteroids
SUFFIX = "_corr" if CORR else ""
rng = np.random.default_rng(42)
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
plates = pd.read_parquet(ROOT / ("data/plates_corrected.parquet" if CORR else "data/plates.parquet"))
strict = set(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv")["obs_date"].astype(str).str[:10])
cat = cat[cat.obs_date.astype(str).str[:10].isin(strict)].copy()
plates = plates[plates.index.isin(cat.plate.unique())].copy()
print(f"strict-date sample: {len(cat)} transients on {len(plates)} plates, {cat.obs_date.nunique()} nights")

HALF, STEP = 3.0, 0.1   # the central 6 x 6 deg of each plate, as in shadow.py
cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].copy()
g = np.arange(-HALF + STEP / 2, HALF, STEP)
GX, GY = [a.ravel() for a in np.meshgrid(g, g)]

# stacked plate-frame density template from all transients (shadow holds <1%, negligible)
tmpl, _, _ = np.histogram2d(cat.eta.clip(-HALF + 1e-6, HALF - 1e-6), cat.xi.clip(-HALF + 1e-6, HALF - 1e-6),
                            bins=[np.append(g - STEP / 2, HALF), np.append(g - STEP / 2, HALF)])
W_TMPL = tmpl.ravel() / tmpl.sum()  # rows: eta (GY), cols: xi (GX) -> matches meshgrid ravel order


def deproject(xi, eta, ra0, dec0):
    xi, eta, ra0, dec0 = map(np.radians, (xi, eta, ra0, dec0))
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def separations(ra, dec, ut, geom, dra=0.0):
    t = Time(ut)
    if geom == "geo":
        return antisun_geocentric((ra - dra) % 360, dec, t)
    return sep_topocentric((ra - dra) % 360, dec, t)  # rotating the sky by dra == moving the shadow by +dra


def plate_table(geom, dra=0.0):
    rows = []
    for name, p in plates.iterrows():
        ra, dec = deproject(GX, GY, p.ra_c, p.dec_c)
        ins = separations(ra, dec, np.full(len(ra), p.ut), geom, dra) < SHADOW_DEG
        rows.append((name, ins.mean(), np.sum(ins * W_TMPL)))
    pt = pd.DataFrame(rows, columns=["plate", "f_area", "f_tmpl"]).set_index("plate")
    sub = cat[["plate", "ra", "dec"]].join(plates[["ut"]], on="plate")
    ins = separations(sub.ra.values, sub.dec.values, sub.ut.values, geom, dra) < SHADOW_DEG
    pt["N"] = cat.groupby("plate").size()
    pt["n_in"] = pd.Series(ins, index=sub.index).groupby(sub.plate).sum()
    pt["night"] = plates.ut.dt.strftime("%Y-%m-%d")
    return pt


def summary(pt, label, nboot=4000):
    N, n_in = pt.N.sum(), pt.n_in.sum()
    e_area = pt.f_area.mean()  # published style: equal weight per plate
    z_binom = (n_in - N * e_area) / np.sqrt(N * e_area * (1 - e_area))
    E_w = (pt.N * pt.f_tmpl).sum()  # within-plate expectation
    strad = pt[(pt.f_area > 0) & (pt.f_area < 1)]
    z_within = (n_in - E_w) / np.sqrt((pt.N * pt.f_tmpl * (1 - pt.f_tmpl)).sum())
    # cluster bootstrap over nights of the published statistic r = (n_in/N)/e_area
    nights = pt.night.unique()
    by = {k: v for k, v in pt.groupby("night")}
    r = np.empty(nboot)
    for b in range(nboot):
        s = pd.concat([by[k] for k in rng.choice(nights, len(nights))])
        r[b] = (s.n_in.sum() / s.N.sum()) / s.f_area.mean() if s.f_area.mean() > 0 else np.nan
    r_obs = (n_in / N) / e_area
    se = np.nanstd(r)
    print(f"\n[{label}] plates touching shadow {np.sum(pt.f_area > 0)}, straddling {len(strad)}, fully inside {np.sum(pt.f_area == 1)}")
    print(f"  observed {n_in}/{N} = {n_in/N*100:.3f}%   area-expected {e_area*100:.3f}% -> {N*e_area:.0f}   ratio {r_obs:.2f}")
    print(f"  per-transient binomial z = {z_binom:.2f}")
    print(f"  night-bootstrap: ratio {r_obs:.2f} +/- {se:.2f}  ->  z = {(r_obs-1)/se:.2f}   (P[ratio>=1] = {np.mean(r >= 1):.3f})")
    print(f"  within-plate (template-weighted) expected {E_w:.0f}: ratio {n_in/E_w:.2f}, binomial z = {z_within:.2f}")
    print(f"  transient-weighted vs area-weighted expectation: {E_w/N*100:.3f}% vs {e_area*100:.3f}%")
    return r_obs


if __name__ == "__main__":  # noqa
    for geom in ("geo", "topo"):
        pt = plate_table(geom)
        pt.to_parquet(ROOT / f"data/plate_shadow_{geom}{SUFFIX}.parquet")
        summary(pt, geom + SUFFIX)
