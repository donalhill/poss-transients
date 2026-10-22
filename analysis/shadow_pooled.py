"""Earth's shadow, pooled over the plates it crosses: VASCO's count of transients in shadow against an average plate, and
how far that count moves when the real shadow outlines are laid on plates the shadow never crossed.

Each plate is its central 6 x 6 deg (shadow.py). Each transient counts by VASCO-ML's rating of how likely it is to be
real (sample.py), as Bruehl et al. (2026) weight theirs; VASCO's own 2025 construction, reproduced, counts each once.
For each sample and geometry: the share of each plate in shadow, its transients and how many lie in shadow; the expected number in shadow (the shadowed plates' worth of sky times the
average plate); the spread VASCO's binomial model gives; and the share of 20,000 random placements (each shadow outline
on a plate neither shadow, geocentric or topocentric, crosses, turned and flipped at random) with as few transients in shadow. Also on plates near the shadow on
the sky only, and within 2 deg of each plate's centre (VASCO's own check, away from the edges). Also how close the
host plates come to the shadow.

Also the random placements with each count scaled by the pointing curve (plate_properties.pointing): what the host plate would
hold had it pointed where the shadow plate did.

Also each tenth of the ratings alone (VASCO-ML's decile edges; their shadow analysis, results/shadow_deficit_deciles,
finds the largest shortfall in the top tenth), each transient counted once and tested in the same way, with how often
chance makes the most extreme of the ten bands look as extreme as the top one; and the page's test at VASCO-ML's shadow
radius (shadow.SHADOW_DEG_ML).

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
from shadow import CELL, GX, GY, HALF, SHADOW_DEG, SHADOW_DEG_ML, antisun_geocentric, deproject, sep_topocentric, unit

ROOT = Path(__file__).resolve().parents[1]
NREP = 20000
NEAR = 20.0                          # deg: "near" host plates lie within this of the anti-sun point
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


def pooled(cat, plates, geom, radius=SHADOW_DEG, weighted=True):
    """Each transient counted by its rating (weighted) or once (VASCO's 2025 construction). VASCO's binomial spread treats
    every transient as independent: with ratings, the spread of a sum of independent weights."""
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].copy()
    cat["u"] = cat.w if weighted else 1.0
    ut = cat.plate.map(plates.ut)
    cat["dark"] = sep(cat.ra.values, cat.dec.values, ut.values, geom) < radius
    cat["cell"] = cells(cat.xi.values, cat.eta.values)
    masks, clear = {}, {}
    for q, p in plates.iterrows():
        ra, dec = deproject(GX, GY, p.ra_c, p.dec_c)
        d = sep(ra, dec, p.ut, geom); m = d < radius
        clear[q] = float(d.min() - radius)                                  # how far the plate stays from the shadow
        if m.any():
            masks[q] = m
    t = pd.DataFrame({"f_area": pd.Series({q: m.mean() for q, m in masks.items()}),
                      "N": cat.groupby("plate").size(), "n_in": cat.groupby("plate").dark.sum(),
                      "W": cat.groupby("plate").u.sum(), "w_in": cat[cat.dark].groupby("plate").u.sum()}).reindex(plates.index)
    t = t.fillna({"f_area": 0.0, "N": 0, "n_in": 0, "W": 0.0, "w_in": 0.0}).astype({"N": int, "n_in": int})
    avg = t.W.mean(); expected = t.f_area.sum() * avg; observed = float(t.w_in[t.f_area > 0].sum())
    share = expected / t.W.sum()
    res = {"weighted": weighted, "plates": len(masks), "whole": [q for q, m in masks.items() if m.all()], "plates_worth": round(float(t.f_area.sum()), 2),
           "avg_per_plate": round(float(avg), 1), "observed": round(observed, 1), "expected": round(float(expected), 1),
           "binom_sd": round(float(np.sqrt((cat.u ** 2).sum() * share * (1 - share))), 1),
           "transients": int(t.n_in[t.f_area > 0].sum())}                       # the transients in shadow, each counted once
    res["binom_z"] = round((observed - expected) / res["binom_sd"], 1)
    return res, t.assign(clear_deg=pd.Series(clear)), cat, masks


def counts(cat, plates, masks, inner=None):
    """Every shadow outline, in each of its eight turns and flips, laid on every plate: the transients it covers (each by
    its weight u, from pooled), as an array outline x orientation x plate (plates in the order of plates.index)."""
    qi = pd.Series(np.arange(len(plates)), plates.index)
    H = np.zeros((len(plates), NG * NG)); np.add.at(H, (qi[cat.plate].values, cat.cell.values), cat.u.values)
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

    # VASCO's 2025 construction: their plates, on archive dates, each transient counted once
    vcat, vplates = sample.vasco()
    out["vasco_geo"], *_ = pooled(vcat, vplates, "geo", weighted=False)

    # ours: every transient weighted by its rating
    cat, plates = sample.ours()
    both = {geom: pooled(cat, plates, geom) for geom in ("geo", "topo")}
    # outlines go only on plates neither shadow (geocentric or topocentric) crosses, so no placement can hold a real effect
    crossed = set(both["geo"][3]) | set(both["topo"][3])
    allp = np.array([q for q in plates.index if q not in crossed])
    clear = np.minimum(both["geo"][1].clear_deg, both["topo"][1].clear_deg).reindex(allp)
    assert (clear > 0).all()
    # each plate's summed ratings predicted by its pointing, from the curve fitted without the plates either shadow
    # crosses, so they cannot shape their own expected counts
    P_ = PP.properties(cat, plates)
    E_pt = pd.Series(PP.pointing(P_, use=~P_.index.isin(crossed), y="nw")[0], P_.index)
    for geom in ("geo", "topo"):
        res, t, c, masks = both[geom]
        res["hosts"] = {"plates": len(allp), "closest_deg": round(float(clear.min()), 2), "within_2deg": int((clear < 2).sum())}
        if geom == "geo":
            sums, res["p_random"], scaled, p_pt = paste(c, plates, masks, res["observed"], lambda k: allp, rng, E=E_pt)
            res["pointing"] = {"expected": round(float(sum(t.f_area[k] * E_pt[k] for k in masks)), 1), "null_mean": round(float(scaled.mean()), 1),
                               "p": round(p_pt, 3), "plates": {k: [round(float(t.W[k]), 1), round(float(E_pt[k]), 1)] for k in masks}}
            # random plates only from those pointed close to the anti-sun point at their own time, as the shadow plates
            # are (their centres lie 4-12 deg from it), with each count scaled by pointing as above
            sun = get_sun(Time(plates.ut.values))
            dA = pd.Series(np.degrees(np.arccos(np.clip(np.sum(unit(plates.ra_c.values, plates.dec_c.values) * -unit(sun.ra.deg, sun.dec.deg), 1), -1, 1))), plates.index)
            nearp = allp[(dA < NEAR).reindex(allp).values]
            *_, res["p_near"] = paste(c, plates, masks, res["observed"], lambda k: nearp, rng, E=E_pt)
            res["near"] = {"max_deg": NEAR, "hosts": len(nearp), "shadow_deg": [round(float(dA[list(masks)].min()), 1), round(float(dA[list(masks)].max()), 1)]}
            # VASCO's check within 2 deg of the plate centre
            c2 = c[np.hypot(c.xi, c.eta) < 2]
            obs2 = float(sum((masks[k][c2.cell[c2.plate == k].values] * c2.u[c2.plate == k].values).sum() for k in masks))
            exp2 = c2.u.sum() * sum(m[R2].sum() / R2.sum() for m in masks.values()) / len(plates)
            _, p2 = paste(c2, plates, masks, obs2, lambda k: allp, rng, inner=R2)
            res["centre"] = {"observed": round(obs2, 1), "expected": round(float(exp2), 1), "p_random": round(p2, 3)}
            # each tenth of the ratings (VASCO-ML's edges), each transient counted once, as their deciles count them, with its
            # own pointing curve. Every band is laid on the same random plates in each trial, so the trials keep how the bands
            # move together: in each trial, the lowest of the ten bands' chances (each from its own trials) gives how often
            # some band looks as extreme as the one picked out
            band = sample.decile(cat.w); edges = np.r_[0.0, sample.decile_edges(), 1.0]
            K = list(masks); qi_ = pd.Series(np.arange(len(plates)), plates.index)
            q = qi_[allp].values[rng.integers(0, len(allp), (NREP, len(K)))]; d = rng.integers(0, 8, (NREP, len(K)))
            res["bands"], P_null = [], []
            for k in range(10):
                r_k, t_k, c_k, _ = pooled(cat[band == k], plates, "geo", weighted=False)
                E_k = pd.Series(PP.pointing(P_.assign(n=t_k.N.reindex(P_.index)), use=~P_.index.isin(crossed))[0], P_.index)
                e_k = E_k.reindex(plates.index).values
                n_k = counts(c_k, plates, masks)[np.arange(len(K))[None, :], d, q]
                s_k = (n_k * E_k[K].values[None, :] / e_k[q]).sum(1)
                P_null.append(np.searchsorted(np.sort(s_k), s_k, side="right") / NREP)    # each trial's own chance in this band
                res["bands"].append({"band": k, "lo": round(float(edges[k]), 3), "hi": round(float(edges[k + 1]), 3),
                                     "transients": int(t_k.N.sum()), "observed": int(r_k["observed"]), "expected": r_k["expected"],
                                     "binom_z": r_k["binom_z"], "pointing_expected": round(float(sum(t_k.f_area[j] * E_k[j] for j in K)), 1),
                                     "p_pointing": round(float(np.mean(s_k <= r_k["observed"])), 4)})
            lowest = np.min(P_null, 0)
            top = res["bands"][-1]
            res["top_tenth"] = {**top, "cut": top["lo"], "p_any_band": round(float(np.mean(lowest <= top["p_pointing"])), 3),
                                "lowest_band": int(np.argmin([b_["p_pointing"] for b_ in res["bands"]]))}
            # the same test at VASCO-ML's radius: its outlines lie inside the page's, so the same host plates are clear of them
            r_, t_, c_, m_ = pooled(cat, plates, "geo", SHADOW_DEG_ML)
            *_, p_ml = paste(c_, plates, m_, r_["observed"], lambda k: allp, rng, E=E_pt)
            res["radius_ml"] = {"deg": SHADOW_DEG_ML, "plates": r_["plates"], "observed": r_["observed"], "expected": r_["expected"],
                                "p_pointing": round(p_ml, 3)}
            t.assign(f_area_topo=both["topo"][1].f_area).to_parquet(ROOT / "data/shadow_plates.parquet")
            np.save(ROOT / "data/shadow_paste_sums.npy", sums); np.save(ROOT / "data/shadow_paste_scaled.npy", scaled)
        else:
            sums, res["p_random"] = paste(c, plates, masks, res["observed"], lambda k: allp, rng)
        res["p_random"] = round(res["p_random"], 3); res["p_near"] = round(res.get("p_near", np.nan), 3) if geom == "geo" else None
        out[f"ours_{geom}"] = res

    (ROOT / "data/shadow_pooled.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
