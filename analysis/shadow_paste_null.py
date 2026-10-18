"""Plate-level checks of the shadow deficit (corrected dates by default; LL_DATES=header for header dates, LL_GEOM=topo for topocentric).

1. Template variants for the within-plate density estimator phi (wider grid, template excluding shadow plates, uniform).
2. Watters-style footprint paste on V: real shadow footprints (plate frame) pasted onto random plates, with random
   dihedral orientation -> null for the pooled published-style statistic AND for phi_hat (calibrates the bootstrap).
3. VASCO's r < 2 deg centre-of-plate test at plate level.
4. Coverage-hole analogue: N of shadow plates vs plates within 20 deg; out-of-shadow density of straddlers.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from shadow import SHADOW_DEG, antisun_geocentric, sep_topocentric, unit

HDR = os.environ.get("LL_DATES") == "header"
GEOM = os.environ.get("LL_GEOM", "geo")
rng = np.random.default_rng(11)
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
plates = pd.read_parquet(ROOT / ("data/plates.parquet" if HDR else "data/plates_corrected.parquet"))
strict = set(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv")["obs_date"].astype(str).str[:10])
cat = cat[cat.obs_date.astype(str).str[:10].isin(strict)].copy()
plates = plates[plates.index.isin(cat.plate.unique())].copy()
print(f"dates={'header' if HDR else 'corrected'} geom={GEOM}: {len(cat)} transients, {len(plates)} plates")

HALF, STEP = 3.6, 0.1
g = np.arange(-HALF + STEP / 2, HALF, STEP)
GX, GY = [a.ravel() for a in np.meshgrid(g, g)]
NG = len(g)
SQ31 = (np.abs(GX) < 3.1) & (np.abs(GY) < 3.1)  # control-sample footprint (published area weighting)


def deproject(xi, eta, ra0, dec0):
    xi, eta, ra0, dec0 = map(np.radians, (xi, eta, ra0, dec0))
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def sep(ra, dec, ut):
    t = Time(np.full(len(ra), ut) if np.ndim(ut) == 0 else ut)
    return antisun_geocentric(ra, dec, t) if GEOM == "geo" else sep_topocentric(ra, dec, t)


def cell(xi, eta):
    i = np.floor((eta + HALF) / STEP).astype(int)
    j = np.floor((xi + HALF) / STEP).astype(int)
    ok = (i >= 0) & (i < NG) & (j >= 0) & (j < NG)
    return np.where(ok, i * NG + j, -1)


cat["cell"] = cell(cat.xi.values, cat.eta.values)
print(f"transients outside the {HALF} deg grid: {(cat.cell < 0).sum()}")
sub = cat.join(plates[["ut"]], on="plate")
cat["in_sh"] = sep(sub.ra.values, sub.dec.values, sub.ut.values) < SHADOW_DEG
N = cat.groupby("plate").size().reindex(plates.index)

# shadow masks in the plate frame (grid cells), for plates touching the shadow
masks = {}
for name, p in plates.iterrows():
    if np.degrees(np.arccos(np.clip(unit(p.ra_c, p.dec_c) @ -unit(*[x for x in _sun]) if False else 1, -1, 1))) > 99:
        pass
for name, p in plates.iterrows():
    ra, dec = deproject(GX, GY, p.ra_c, p.dec_c)
    m = sep(ra, dec, p.ut) < SHADOW_DEG
    if m[SQ31].any() or m.any():
        masks[name] = m
sh = list(masks)
print(f"plates with shadow inside the {HALF} deg grid: {len(sh)}; inside the 3.1 deg square: {sum(masks[k][SQ31].any() for k in sh)}")


def template(exclude=()):
    c = cat[(cat.cell >= 0) & ~cat.plate.isin(exclude)]
    w = np.bincount(c.cell, minlength=NG * NG).astype(float)
    return w / w.sum()


T = {"wide (all)": template(), "wide (excl. shadow plates)": template(sh), "uniform 3.1 sq": SQ31 / SQ31.sum()}


def phi_hat(O, Nn, f):
    return 1 - np.sum(O * (1 - f)) / np.sum((Nn - O) * f)


rows = []
for k in sh:
    m = masks[k]
    rows.append(dict(plate=k, N=int(N[k]), O_exact=int(cat.in_sh[cat.plate == k].sum()),
                     O_grid=int(m[cat.cell[(cat.plate == k) & (cat.cell >= 0)]].sum()),
                     f_area31=m[SQ31].mean(), **{f"f[{n}]": np.sum(m * w) for n, w in T.items()}))
tab = pd.DataFrame(rows).set_index("plate").join(plates[["ra_c", "dec_c"]])
print(tab.round(3).to_string())

print("\n== 1. phi_hat (within-plate density ratio), straddlers 0 < f < 0.97, plate-bootstrap sd ==")
for n in T:
    st = tab[(tab[f"f[{n}]"] > 0) & (tab[f"f[{n}]"] < 0.97)]
    O, Nn, f = st.O_exact.values.astype(float), st.N.values.astype(float), st[f"f[{n}]"].values
    bs = [phi_hat(O[i], Nn[i], f[i]) for i in (rng.integers(0, len(st), len(st)) for _ in range(5000))]
    loo = [phi_hat(np.delete(O, i), np.delete(Nn, i), np.delete(f, i)) for i in range(len(st))]
    print(f"  {n:28s} {len(st)} plates: phi = {phi_hat(O, Nn, f):+.3f}  boot sd {np.std(bs):.3f}  "
          f"leave-one-out range {min(loo):+.2f}..{max(loo):+.2f}")
    per = (1 - (O / np.maximum(Nn - O, 1e-9)) * (1 - f) / f)
    print("     per-plate phi:", dict(zip(st.index, per.round(2))))
# binomial heterogeneity across straddlers (is within-plate scatter > binomial?)
n0 = "wide (all)"
st = tab[(tab[f"f[{n0}]"] > 0) & (tab[f"f[{n0}]"] < 0.97)]
E = st.N * st[f"f[{n0}]"]
chi2 = np.sum((st.O_exact - E) ** 2 / (E * (1 - st[f"f[{n0}]"])))
from scipy import stats
print(f"  heterogeneity under no avoidance: chi2 = {chi2:.1f} on {len(st)} dof, p = {stats.chi2.sf(chi2, len(st)):.2g}")

print("\n== 2. Footprint paste (Watters-style) on V ==")
cells_by_plate = {k: v.cell.values[v.cell.values >= 0] for k, v in cat.groupby("plate")}
pnames = np.array(plates.index)
if os.environ.get("LL_PASTE") == "near":  # sky-matched null: only plates 10-45 deg from the anti-sun (opposition region)
    from astropy.coordinates import get_sun
    sun = get_sun(Time(plates.ut.values))
    dA = np.degrees(np.arccos(np.clip(np.sum(unit(plates.ra_c.values, plates.dec_c.values) * -unit(sun.ra.deg, sun.dec.deg), 1), -1, 1)))
    pnames = np.array(plates.index[(dA > 10) & (dA < 45)])
    print(f"  paste candidates restricted to {len(pnames)} plates 10-45 deg from the anti-sun")
W = T["wide (all)"].reshape(NG, NG)


def dihedral(m, k):
    a = m.reshape(NG, NG)
    a = np.rot90(a, k % 4)
    return (a.T if k >= 4 else a).ravel()


obs_sum = tab.O_grid.sum()
mask_list = [masks[k] for k in sh]
strad_keys = [k for k in sh if 0 < np.sum(masks[k] * T[n0]) < 0.97]
NREP = 20000
sums, phis = np.empty(NREP), np.empty(NREP)
draws = {"O": np.empty((NREP, len(strad_keys))), "N": np.empty((NREP, len(strad_keys))), "f": np.empty((NREP, len(strad_keys)))}
for r in range(NREP):
    pick = rng.choice(pnames, len(mask_list), replace=False)
    s = 0
    for m, q in zip(mask_list, pick):
        mm = dihedral(m, rng.integers(8))
        s += mm[cells_by_plate[q]].sum()
    sums[r] = s
    pick = rng.choice(pnames, len(strad_keys), replace=False)
    O_, N_, f_ = [], [], []
    for k, q in zip(strad_keys, pick):
        mm = dihedral(masks[k], rng.integers(8))
        O_.append(mm[cells_by_plate[q]].sum()); N_.append(len(cells_by_plate[q])); f_.append(np.sum(mm * T[n0]))
    O_, N_, f_ = map(np.array, (O_, N_, f_))
    draws["O"][r], draws["N"][r], draws["f"][r] = O_, N_, f_
    phis[r] = phi_hat(O_, N_, f_) if np.sum((N_ - O_) * f_) > 0 else np.nan
E_area = len(cat) * np.mean([masks[k][SQ31].mean() if k in masks else 0 for k in pnames]) if False else \
    len(cat) * sum(masks[k][SQ31].mean() for k in sh) / len(plates)
print(f"  observed in-footprint (grid) {obs_sum}; area-expected {E_area:.0f}; ratio {obs_sum / E_area:.2f}")
print(f"  pasted: mean {sums.mean():.0f}, median {np.median(sums):.0f}, 5-95% {np.percentile(sums, 5):.0f}-{np.percentile(sums, 95):.0f}")
print(f"  P(pasted <= observed) = {np.mean(sums <= obs_sum):.3f}   (one-sided p for the pooled deficit)")
print(f"  equivalent z = {stats.norm.ppf(max(np.mean(sums <= obs_sum), 1 / NREP)):.2f}")
ph = phis[np.isfinite(phis)]
print(f"  phi_hat under paste null ({len(strad_keys)} straddler footprints): mean {ph.mean():+.3f}, sd {ph.std():.3f}, "
      f"5-95% {np.percentile(ph, 5):+.2f}..{np.percentile(ph, 95):+.2f}")
obs_phi = phi_hat(tab.loc[strad_keys].O_exact.values, tab.loc[strad_keys].N.values, tab.loc[strad_keys][f"f[{n0}]"].values)
print(f"  observed phi {obs_phi:+.3f}: P(null phi >= obs) = {np.mean(ph >= obs_phi):.3f}")
# injection on the paste null: how often would a true phi=0.3 be detected at 90% one-sided?
thr = np.percentile(ph, 90)
print(f"  90% one-sided null threshold phi > {thr:+.2f}; 1-(1-null)*(1-0.3) median {np.median(1 - (1 - ph) * 0.7):+.2f} "
      f"-> power vs phi=0.3 ~ {np.mean(1 - (1 - ph) * 0.7 > thr):.2f}, vs 0.5 ~ {np.mean(1 - (1 - ph) * 0.5 > thr):.2f}")
np.save(ROOT / f"data/paste_sums_{GEOM}{os.environ.get('LL_PASTE', '')}{'_hdr' if HDR else ''}.npy", sums)
# straddler footprints on random plates, with the observed counts, for the bounded interval (shadow_phi_interval.py)
np.savez(ROOT / f"data/paste_phi_{GEOM}{os.environ.get('LL_PASTE', '')}{'_hdr' if HDR else ''}.npz", **draws,
         obs_O=tab.loc[strad_keys].O_exact.values, obs_N=tab.loc[strad_keys].N.values, obs_f=tab.loc[strad_keys][f"f[{n0}]"].values)

print("\n== 3. VASCO centre-of-plate test (r < 2 deg) ==")
R2 = np.hypot(GX, GY) < 2
c2 = cat[np.hypot(cat.xi, cat.eta) < 2]
f2 = {k: masks[k][R2].mean() for k in sh}
n2 = c2.in_sh.sum(); e2 = len(c2) * sum(f2.values()) / len(plates)
print(f"  r<2: {len(c2)} transients, {n2} in shadow vs area-expected {e2:.1f} (ratio {n2 / e2:.2f}); plates with r<2 overlap "
      f"{sum(v > 0 for v in f2.values())}: " + ", ".join(f"{k}:{v:.2f}" for k, v in f2.items() if v > 0))
z_binom = (n2 - e2) / np.sqrt(e2)
print(f"  per-transient Poisson z = {z_binom:.2f}")
cells2 = {k: v.cell.values[(v.cell.values >= 0) & (np.hypot(v.xi, v.eta).values < 2)] for k, v in cat.groupby("plate")}
obs2 = sum(masks[k][cells2[k]].sum() for k in sh)
s2 = np.empty(NREP)
for r in range(NREP):
    pick = rng.choice(pnames, len(sh), replace=False)
    s2[r] = sum(dihedral(masks[k] & R2, rng.integers(8))[cells2[q]].sum() for k, q in zip(sh, pick))
print(f"  paste null for r<2: observed {obs2}, null median {np.median(s2):.0f}, P(null <= obs) = {np.mean(s2 <= obs2):.3f}")

print("\n== 4. Coverage hole / low-count shadow plates ==")
v = unit(plates.ra_c.values, plates.dec_c.values)
for k in sh:
    d = np.degrees(np.arccos(np.clip(v @ unit(plates.loc[k, "ra_c"], plates.loc[k, "dec_c"]), -1, 1)))
    nb = (d < 20) & (d > 0)
    f = tab.loc[k, f"f[{n0}]"]
    out_density = (N[k] - tab.loc[k, "O_exact"]) / (1 - f) if f < 0.999 else np.nan
    print(f"  {k}: N={int(N[k]):4d} (pct {int((N < N[k]).mean() * 100):2d}), neighbours<20deg n={nb.sum():2d} median N={np.median(N.values[nb]):.0f}; "
          f"N implied from out-of-shadow density = {out_density:.0f}")
tab.to_csv(ROOT / f"data/shadow_plates_{GEOM}{'_hdr' if HDR else ''}.csv")
