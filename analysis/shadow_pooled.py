"""Earth's shadow, pooled over the plates it crosses: VASCO's count of transients in shadow against an average plate, and
how far that count moves when the real shadow outlines are laid on plates the shadow never crossed.

Each plate is its central 6 x 6 deg (shadow.py). For each sample and geometry: the share of each plate in shadow, its
transients and how many lie in shadow; the expected number in shadow (the shadowed plates' worth of sky times the
average plate); the spread VASCO's binomial model gives; and the share of 20,000 random placements (each shadow outline
on a plate neither shadow, geocentric or topocentric, crosses, turned and flipped at random) with as few transients in shadow. Also on plates near the shadow on
the sky only, and within 2 deg of each plate's centre (VASCO's own check, away from the edges). Also how close the
host plates come to the shadow.

Also the random placements with each count scaled by the pointing curve (plate_properties.pointing): what the host plate would
hold had it pointed where the shadow plate did.

python analysis/shadow_pooled.py -> data/shadow_pooled.json, data/shadow_plates.parquet, data/shadow_paste_sums.npy,
    data/shadow_paste_scaled.npy"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import get_sun
from astropy.time import Time

import plate_properties as PP
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
    masks, clear = {}, {}
    for q, p in plates.iterrows():
        ra, dec = deproject(GX, GY, p.ra_c, p.dec_c)
        d = sep(ra, dec, p.ut, geom); m = d < SHADOW_DEG
        clear[q] = float(d.min() - SHADOW_DEG)                                  # how far the plate stays from the shadow
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
    return res, t.assign(clear_deg=pd.Series(clear)), cat, masks


def counts(cat, plates, masks, inner=None):
    """Every shadow outline, in each of its eight turns and flips, laid on every plate: the transients it covers, as an
    array outline x orientation x plate (plates in the order of plates.index)."""
    qi = pd.Series(np.arange(len(plates)), plates.index)
    H = np.zeros((len(plates), NG * NG)); np.add.at(H, (qi[cat.plate].values, cat.cell.values), 1)
    M = np.array([[dihedral(m if inner is None else m & inner, d) for d in range(8)] for m in masks.values()], float)
    return M @ H.T


def paste(cat, plates, masks, observed, pick_from, rng, inner=None, E=None):
    """Lay each shadow outline on a random plate (turned and flipped at random) and count the transients it covers. With
    E (each plate's count predicted by its pointing), also each count scaled to what that plate would hold had it pointed
    where the shadow plate did."""
    C, qi, K = counts(cat, plates, masks, inner), pd.Series(np.arange(len(plates)), plates.index), list(masks)
    q = np.stack([qi[pick_from(k)].values[rng.integers(0, len(pick_from(k)), NREP)] for k in K], 1)   # trials x outlines
    d = rng.integers(0, 8, (NREP, len(K))); n = C[np.arange(len(K))[None, :], d, q]
    sums = n.sum(1)
    if E is None:
        return sums, float(np.mean(sums <= observed))
    e = E.reindex(plates.index).values; scaled = (n * E[K].values[None, :] / e[q]).sum(1)
    return sums, float(np.mean(sums <= observed)), scaled, float(np.mean(scaled <= observed))


if __name__ == "__main__":
    out = {}
    rng = np.random.default_rng(11)

    # VASCO's construction: their plates, on archive dates
    vcat, vplates = sample.vasco()
    out["vasco_geo"], *_ = pooled(vcat, vplates, "geo")

    cat, plates = sample.ours()
    both = {geom: pooled(cat, plates, geom) for geom in ("geo", "topo")}
    # each plate's count predicted by its pointing, from the curve fitted without the plates either shadow crosses, so they
    # cannot shape their own expected counts
    P_ = PP.properties(cat, plates)
    E_pt = pd.Series(PP.pointing(P_, use=~P_.index.isin(set(both["geo"][3]) | set(both["topo"][3])))[0], P_.index)
    # outlines go only on plates neither shadow (geocentric or topocentric) crosses, so no placement can hold a real effect
    crossed = set(both["geo"][3]) | set(both["topo"][3])
    allp = np.array([q for q in plates.index if q not in crossed])
    clear = np.minimum(both["geo"][1].clear_deg, both["topo"][1].clear_deg).reindex(allp)
    assert (clear > 0).all()
    for geom in ("geo", "topo"):
        res, t, c, masks = both[geom]
        res["hosts"] = {"plates": len(allp), "closest_deg": round(float(clear.min()), 2), "within_2deg": int((clear < 2).sum())}
        if geom == "geo":
            sums, res["p_random"], scaled, res["p_pointing"] = paste(c, plates, masks, res["observed"], lambda k: allp, rng, E=E_pt)
            res["pointing"] = {"expected": round(float(sum(t.f_area[k] * E_pt[k] for k in masks)), 1), "null_mean": round(float(scaled.mean()), 1),
                               "p": round(res.pop("p_pointing"), 3), "plates": {k: [int(t.N[k]), round(float(E_pt[k]), 1)] for k in masks}}
        else:
            sums, res["p_random"] = paste(c, plates, masks, res["observed"], lambda k: allp, rng)
        if geom == "geo":
            # random plates only from those 10-45 deg from the anti-sun point at their own time, as the shadow plates are
            sun = get_sun(Time(plates.ut.values))
            dA = pd.Series(np.degrees(np.arccos(np.clip(np.sum(unit(plates.ra_c.values, plates.dec_c.values) * -unit(sun.ra.deg, sun.dec.deg), 1), -1, 1))), plates.index)
            nearp = allp[((dA > 10) & (dA < 45)).reindex(allp).values]
            _, res["p_near"] = paste(c, plates, masks, res["observed"], lambda k: nearp, rng)
            # VASCO's check within 2 deg of the plate centre
            c2 = c[np.hypot(c.xi, c.eta) < 2]
            obs2 = int(sum(masks[k][c2.cell[c2.plate == k].values].sum() for k in masks))
            exp2 = len(c2) * sum(m[R2].sum() / R2.sum() for m in masks.values()) / len(plates)
            _, p2 = paste(c2, plates, masks, obs2, lambda k: allp, rng, inner=R2)
            res["centre"] = {"observed": obs2, "expected": round(float(exp2), 1), "p_random": round(p2, 3)}
            t.assign(f_area_topo=both["topo"][1].f_area).to_parquet(ROOT / "data/shadow_plates.parquet")
            np.save(ROOT / "data/shadow_paste_sums.npy", sums); np.save(ROOT / "data/shadow_paste_scaled.npy", scaled)
        res["p_random"] = round(res["p_random"], 3); res["p_near"] = round(res.get("p_near", np.nan), 3) if geom == "geo" else None
        out[f"ours_{geom}"] = res

    (ROOT / "data/shadow_pooled.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
