"""Earth's shadow as a ratio, measured plate by plate: the density of transients inside the shadow divided by the density
in sunlight on the same plate (1 = no effect). Uses nothing beyond the shadow plates themselves.

On each plate partly in shadow, the plate is split into rings by distance from its centre (Chebyshev, as the plate is
square); in each ring the sunlit density times the shadowed area gives the count expected in shadow, so the rise of
transients towards plate edges cancels within the plate. r = observed in shadow / expected.
Uncertainty per plate: block bootstrap, resampling the plate's 0.5 deg patches with replacement (keeps the clumping of
transients within a plate). Combined: total observed / total expected over plates, with every plate resampled
independently; a chi-square says whether the plates agree with one another.
Sets: VASCO's catalogue on its nine shadow plates; our own extraction on the southern plates, on VASCO's nine plates, and
on both together (shadow_ours.py writes data/shadow_ours_grid.npz).
python analysis/shadow_ratio.py   -> data/shadow_ratio.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from shadow import SHADOW_DEG, antisun_geocentric   # noqa: E402

rng = np.random.default_rng(5)
NBOOT, STEP, BLOCK, NRING = 4000, 0.05, 0.5, 10


def deproject(xi, eta, ra0, dec0):
    """Plate-frame offsets (deg) -> RA, Dec (deg), gnomonic, as in shadow_paste_null.py."""
    xi, eta, ra0, dec0 = map(np.radians, (xi, eta, ra0, dec0))
    rho = np.hypot(xi, eta); c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def vasco_grid(half=3.6):
    """VASCO's catalogue on its nine shadow plates (corrected dates, nights in VASCO's list)."""
    cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
    strict = set(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date.astype(str).str[:10])
    cat = cat[cat.obs_date.astype(str).str[:10].isin(strict)]
    pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
    sh = pd.read_parquet(ROOT / "data/plate_shadow_geo_corr.parquet")
    names = sh.index[sh.f_area > 0].tolist()
    g = np.arange(-half + STEP / 2, half, STEP); GX, GY = np.meshgrid(g, g); ng = len(g)
    masks, cells = [], []
    for q in names:
        p = pl.loc[q]
        ra, dec = deproject(GX.ravel(), GY.ravel(), p.ra_c, p.dec_c)
        masks.append((antisun_geocentric(ra, dec, Time(np.full(ra.size, pd.Timestamp(p.ut).to_datetime64()))) < SHADOW_DEG).reshape(ng, ng))
        c = cat[cat.plate == q]
        ix, iy = np.floor((c.xi.values + half) / STEP).astype(int), np.floor((c.eta.values + half) / STEP).astype(int)
        ok = (ix >= 0) & (ix < ng) & (iy >= 0) & (iy < ng)
        cells.append(np.ravel_multi_index((iy[ok], ix[ok]), (ng, ng)))
    return dict(plates=np.array(names), sets=np.array(["vasco"] * len(names)), masks=np.array(masks),
                usable=np.ones((len(names), ng, ng), bool), cells=cells, half=half)


def plate_tables(G, q):
    """Per 0.5 deg block and ring: transients and cells, in shadow and in sunlight -> array (blocks, rings, 4)."""
    m, u = G["masks"][q], G["usable"][q]
    ng = m.shape[0]; half = ng * STEP / 2
    g = np.arange(-half + STEP / 2, half, STEP); GX, GY = np.meshgrid(g, g)
    ring = np.clip((np.maximum(np.abs(GX), np.abs(GY)) / half * NRING).astype(int), 0, NRING - 1)
    nb = int(np.ceil(ng * STEP / BLOCK)); per = int(round(BLOCK / STEP))
    block = (np.arange(ng)[:, None] // per) * nb + (np.arange(ng)[None, :] // per)
    count = np.bincount(np.asarray(G["cells"][q], dtype=int), minlength=ng * ng).reshape(ng, ng)
    T = np.zeros((nb * nb, NRING, 4))
    for k, (sel, val) in enumerate([(m & u, count), (~m & u, count), (m & u, 1), (~m & u, 1)]):
        np.add.at(T[:, :, k], (block[sel], ring[sel]), val[sel] if np.ndim(val) else 1)
    return T[T[:, :, 2:].sum(axis=(1, 2)) > 0]          # blocks with usable sky


def obs_exp(T):
    """T summed over blocks -> (observed in shadow, expected from the sunlit density ring by ring)."""
    n_sh, n_sun, a_sh, a_sun = (T[..., k] for k in range(4))
    ok = a_sun > 0
    return np.where(ok, n_sh, 0).sum(-1), np.where(ok, n_sun * a_sh / np.where(ok, a_sun, 1), 0).sum(-1)


def measure(G, which, label):
    rows, boots = [], []
    for q in range(len(G["plates"])):
        if G["sets"][q] not in which:
            continue
        T = plate_tables(G, q)
        o, e = obs_exp(T.sum(0))
        frac = T[..., 2].sum() / T[..., 2:].sum()
        if e <= 0 or frac >= 0.97:                       # no sunlit sky to compare with
            continue
        b = obs_exp(T[rng.integers(0, len(T), (NBOOT, len(T)))].sum(1))
        boots.append(np.c_[b[0], b[1]])
        rb = b[0] / np.where(b[1] > 0, b[1], np.nan)
        rows.append(dict(plate=str(G["plates"][q]), set=str(G["sets"][q]), shadow=round(float(frac), 3), obs=int(o),
                         exp=round(float(e), 1), r=round(float(o / e), 3), sd=round(float(np.nanstd(rb)), 3)))
    t = pd.DataFrame(rows)
    O, E = t.obs.sum(), t.exp.sum()
    bO = sum(b[:, 0] for b in boots); bE = sum(b[:, 1] for b in boots)
    r, sd = O / E, float(np.std(bO / bE))
    # agreement between plates: each plate's count against the combined ratio, with its bootstrap variance
    var = np.array([np.var(b[:, 0] - r * b[:, 1]) for b in boots])
    chi2 = float(np.sum((t.obs - r * t.exp) ** 2 / np.maximum(var, 1)))
    res = {"label": label, "plates": rows, "obs": int(O), "exp": round(float(E), 1), "r": round(float(r), 3),
           "sd": round(sd, 3), "sigma_from_1": round(float((1 - r) / sd), 2),
           "sigma_from_068": round(float((r - 0.68) / sd), 2), "chi2": round(chi2, 1), "ndof": len(t) - 1}
    print(f"\n{label}: {len(t)} plates partly in shadow, observed {O}, expected {E:.0f}, r = {r:.3f} +- {sd:.3f} "
          f"({(r - 1) / sd:+.1f} sigma from 1, {(r - 0.68) / sd:+.1f} sigma from 0.68); chi2 {chi2:.1f} / {len(t) - 1}")
    print(t.to_string(index=False))
    return res


s = np.load(ROOT / "data/shadow_ours_grid.npz", allow_pickle=True)
OURS = {k: s[k] for k in ("plates", "sets", "masks", "usable", "cells")}
out = {"vasco": measure(vasco_grid(), ["vasco"], "VASCO's plates, VASCO's catalogue"),
       "south": measure(OURS, ["south"], "Southern plates, our transients")}
if "north" in set(OURS["sets"]):
    out["north"] = measure(OURS, ["north"], "VASCO's plates, our transients")
out["combined"] = measure(OURS, ["south", "north"], "All plates, our transients")
(ROOT / "data/shadow_ratio.json").write_text(json.dumps(out))
