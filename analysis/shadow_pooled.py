"""Earth's shadow, pooled over the plates it crosses: VASCO's count of transients in shadow against an average plate, and
how far that count moves when the real shadow outlines are laid on plates the shadow never crossed.

Each plate is its central 6 x 6 deg (shadow.py). For each sample and geometry: the share of each plate in shadow, its
transients and how many lie in shadow; the expected number in shadow (the shadowed plates' worth of sky times the
average plate); the spread VASCO's binomial model gives; and the share of 20,000 random placements (each shadow outline
on another plate, turned and flipped at random) with as few transients in shadow. Also on plates near the shadow on
the sky only, and within 2 deg of each plate's centre (VASCO's own check, away from the edges).

python analysis/shadow_pooled.py -> data/shadow_pooled.json, data/shadow_plates.parquet, data/shadow_paste_sums.npy"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import get_sun
from astropy.time import Time

import sample
from shadow import CELL, GX, GY, HALF, SHADOW_DEG, antisun_geocentric, deproject, sep_topocentric, unit

ROOT = Path(__file__).resolve().parents[1]
NREP = 20000
NG = int(round(2 * HALF / CELL))
R2 = np.hypot(GX, GY) < 2


def sep(ra, dec, ut, geom):
    t = Time(np.broadcast_to(np.asarray(ut, dtype="datetime64[ns]"), len(ra)))
    return antisun_geocentric(ra, dec, t) if geom == "geo" else sep_topocentric(ra, dec, t)


def cells(xi, eta):
    i, j = np.floor((eta + HALF) / CELL).astype(int), np.floor((xi + HALF) / CELL).astype(int)
    return i * NG + j


def dihedral(m, k):
    a = np.rot90(m.reshape(NG, NG), k % 4)
    return (a.T if k >= 4 else a).ravel()


def pooled(cat, plates, geom):
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].copy()
    ut = cat.plate.map(plates.ut)
    cat["dark"] = sep(cat.ra.values, cat.dec.values, ut.values, geom) < SHADOW_DEG
    cat["cell"] = cells(cat.xi.values, cat.eta.values)
    masks = {}
    for q, p in plates.iterrows():
        ra, dec = deproject(GX, GY, p.ra_c, p.dec_c)
        m = sep(ra, dec, p.ut, geom) < SHADOW_DEG
        if m.any():
            masks[q] = m
    t = pd.DataFrame({"f_area": pd.Series({q: m.mean() for q, m in masks.items()}),
                      "N": cat.groupby("plate").size(), "n_in": cat.groupby("plate").dark.sum()}).reindex(plates.index)
    t = t.fillna({"f_area": 0.0, "N": 0, "n_in": 0}).astype({"N": int, "n_in": int})
    avg = t.N.mean(); expected = t.f_area.sum() * avg; observed = int(t.n_in[t.f_area > 0].sum())
    share = expected / t.N.sum()
    res = {"plates": len(masks), "whole": [q for q, m in masks.items() if m.all()], "plates_worth": round(float(t.f_area.sum()), 2),
           "avg_per_plate": round(float(avg), 1), "observed": observed, "expected": round(float(expected), 1),
           "binom_sd": round(float(np.sqrt(t.N.sum() * share * (1 - share))), 1)}
    res["binom_z"] = round((observed - expected) / res["binom_sd"], 1)
    return res, t, cat, masks


def paste(cat, plates, masks, observed, pick_from, rng, inner=None):
    """Lay each shadow outline on a random plate (turned and flipped at random) and count the transients it covers."""
    by = {q: g.cell.values for q, g in cat.groupby("plate")}
    sums = np.empty(NREP)
    for r in range(NREP):
        s = 0
        for k, m in masks.items():
            q = rng.choice(pick_from(k))
            mm = dihedral(m if inner is None else m & inner, rng.integers(8))
            s += mm[by.get(q, np.empty(0, int))].sum()
        sums[r] = s
    return sums, float(np.mean(sums <= observed))


out = {}
rng = np.random.default_rng(11)

# VASCO's construction: their plates, on archive dates
vcat, vplates = sample.vasco()
out["vasco_geo"], *_ = pooled(vcat, vplates, "geo")

cat, plates = sample.ours()
for geom in ("geo", "topo"):
    res, t, c, masks = pooled(cat, plates, geom)
    allp = np.array(plates.index)
    sums, res["p_random"] = paste(c, plates, masks, res["observed"], lambda k: allp, rng)
    if geom == "geo":
        # random plates only from those 10-45 deg from the anti-sun point at their own time, as the shadow plates are
        sun = get_sun(Time(plates.ut.values))
        dA = np.degrees(np.arccos(np.clip(np.sum(unit(plates.ra_c.values, plates.dec_c.values) * -unit(sun.ra.deg, sun.dec.deg), 1), -1, 1)))
        nearp = allp[(dA > 10) & (dA < 45)]
        _, res["p_near"] = paste(c, plates, masks, res["observed"], lambda k: nearp, rng)
        # VASCO's check within 2 deg of the plate centre
        c2 = c[np.hypot(c.xi, c.eta) < 2]
        obs2 = int(sum(masks[k][c2.cell[c2.plate == k].values].sum() for k in masks))
        exp2 = len(c2) * sum(m[R2].sum() / R2.sum() for m in masks.values()) / len(plates)
        _, p2 = paste(c2, plates, masks, obs2, lambda k: allp, rng, inner=R2)
        res["centre"] = {"observed": obs2, "expected": round(float(exp2), 1), "p_random": round(p2, 3)}
        t.to_parquet(ROOT / "data/shadow_plates.parquet")
        np.save(ROOT / "data/shadow_paste_sums.npy", sums)
    res["p_random"] = round(res["p_random"], 3); res["p_near"] = round(res.get("p_near", np.nan), 3) if geom == "geo" else None
    out[f"ours_{geom}"] = res

(ROOT / "data/shadow_pooled.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
