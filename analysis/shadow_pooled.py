"""Earth's shadow, pooled over the plates it crosses: VASCO's count of transients in shadow against an average plate, and
how far that count moves when the real shadow outlines are laid on plates the shadow never crossed.

The shadow follows shadow.py's model (shadow.MAIN): seen from Palomar, where it fell on the plates, each transient
counted by the share of its plate's exposure that its place spent in shadow (the shadow slides across a plate during an
exposure), with the edge from the share of the Sun's disc the Earth hides. VASCO's 2025 construction is reproduced as
they made it (their plates, archive dates and times, each transient counted once, a sharp 8.69 deg edge at one
instant) with their binomial, seen from Palomar and, as they place it, from the centre of the Earth;
ours, each transient counted by VASCO-ML's rating (sample.py), on the central 6 x 6 deg of each plate: the count in shadow
against an average plate and VASCO's binomial spread; the share of 20,000 random placements (each outline on a plate no
shadow reaches, turned and flipped at random) with as few in shadow, without and with each count scaled by the pointing
model (plate_properties.pointing) to the shadow plate's pointing; the test's power for VASCO's shortfall measured from
Palomar; random plates near the anti-sun point only; the central 2 deg (VASCO's check); each tenth of the ratings,
with how often the most extreme of ten looks as extreme as the top one; VASCO-ML's 8.5 deg radius; and the model's
systematics (shadow.variants). The host plates and the pointing model are shared by every variant.

python analysis/shadow_pooled.py -> data/shadow_pooled.json, data/shadow_plates.parquet,
    data/shadow_paste_sums.npy, data/shadow_paste_scaled.npy"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import get_sun
from astropy.time import Time

import plate_properties as PP
import sample
import shadow as SH
from shadow import GX, GY, HALF, SHADOW_DEG_ML, Shadow, turn, unit

ROOT = Path(__file__).resolve().parents[1]
NREP = 20000
NEAR = 20.0                          # deg: "near" host plates lie within this of the anti-sun point
R2 = np.hypot(GX, GY) < 2


class Outlines(dict):
    """Each shadow plate's share map (cells), keyed by plate, with the shadow model S that made them."""
    def __init__(self, S, maps):
        super().__init__(maps); self.S = S


def pooled(cat, plates, S, weighted=True, exact=False):
    """Each transient counted by its rating (weighted) or once (VASCO's 2025 construction), times the share of its
    plate's exposure that its place spent in shadow under model S; a plate's shadowed area is that share integrated over
    its face. exact: VASCO's construction as they made it, each cell's share from its centre (with S a sharp edge at one
    instant). VASCO's binomial spread treats every transient as independent: with ratings, the spread of a sum of
    independent weights."""
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].copy()
    cat["u"] = cat.w if weighted else 1.0
    masks = {}
    for q, p in plates.iterrows():
        m = SH.cell_dark(p, S, centres=exact)
        if m.any():
            masks[q] = m
    s = np.zeros(len(cat))
    for q in masks:
        i = (cat.plate == q).values; p = plates.loc[q]
        s[i] = SH.dark_dirs(unit(cat.ra.values[i], cat.dec.values[i]), p, S) if exact else SH.dark(cat.xi.values[i], cat.eta.values[i], p, S)
    cat["s"] = s; cat["us"] = cat.u * s
    g = cat.groupby("plate")
    t = pd.DataFrame({"f_area": pd.Series({q: m.mean() for q, m in masks.items()}, dtype=float), "N": g.size(), "n_in": g.s.sum(),
                      "W": g.u.sum(), "w_in": g.us.sum()}).reindex(plates.index)
    t = t.fillna({"f_area": 0.0, "N": 0, "n_in": 0.0, "W": 0.0, "w_in": 0.0}).astype({"N": int})
    avg = t.W.mean(); expected = t.f_area.sum() * avg; observed = float(t.w_in[t.f_area > 0].sum())
    share = expected / t.W.sum()
    res = {"weighted": weighted, "plates": len(masks), "whole": [q for q, m in masks.items() if (m > 1 - 1e-9).all()],
           "plates_worth": round(float(t.f_area.sum()), 2), "avg_per_plate": round(float(avg), 1), "observed": round(observed, 1),
           "expected": round(float(expected), 1), "binom_sd": round(float(np.sqrt((cat.u ** 2).sum() * share * (1 - share))), 1),
           "transients": int(round(t.n_in[t.f_area > 0].sum()))}            # the transients in shadow, each counted once
    res["binom_z"] = round((observed - expected) / res["binom_sd"], 1)
    return res, t, cat, Outlines(S, masks)


_LAID = {}


def laid_one(cat, p, S):
    """One shadow plate's shadow (model S) laid on every transient of cat in each of the eight turns and flips: the share
    in shadow at the transient's place on the shadow plate (orientation x transient)."""
    xy = [turn(cat.xi.values, cat.eta.values, d) for d in range(8)]
    return SH.dark(np.concatenate([x for x, _ in xy]), np.concatenate([y for _, y in xy]), p, S).reshape(8, len(cat)).astype(np.float32)


def counts(cat, plates, masks):
    """Every shadow, in each of its eight turns and flips, laid on every plate: the transients it covers, each by its
    weight u (from pooled) times the share in shadow at its place, as an array outline x orientation x plate (plates in
    the order of plates.index). The laid shares are kept for later subsets of the same transients (the rating bands) when
    there are few outlines."""
    K = list(masks); key = (masks.S, tuple(K))
    qi = pd.Series(np.arange(len(plates)), plates.index)[cat.plate].values
    keep = key in _LAID and cat.index.isin(_LAID[key][1]).all()
    if not keep and len(K) <= 40:
        _LAID.clear(); _LAID[key] = (np.stack([laid_one(cat, plates.loc[q], masks.S) for q in K]), cat.index); keep = True
    C = np.zeros((len(K), 8, len(plates)))
    for k, q in enumerate(K):
        A = _LAID[key][0][k][:, _LAID[key][1].get_indexer(cat.index)] if keep else laid_one(cat, plates.loc[q], masks.S)
        for d in range(8):
            C[k, d] = np.bincount(qi, weights=A[d] * cat.u.values, minlength=len(plates))
    return C


def paste(cat, plates, masks, observed, pick_from, rng, E=None):
    """Lay each shadow outline on a random plate (turned and flipped at random) and count the transients it covers. With
    E (each plate's count predicted by its pointing), also each count scaled to what that plate would hold had it pointed
    where the shadow plate did."""
    C, qi, K = counts(cat, plates, masks), pd.Series(np.arange(len(plates)), plates.index), list(masks)
    q = np.stack([qi[pick_from(k)].values[rng.integers(0, len(pick_from(k)), NREP)] for k in K], 1)   # trials x outlines
    d = rng.integers(0, 8, (NREP, len(K))); n = C[np.arange(len(K))[None, :], d, q]
    sums = n.sum(1)
    if E is None:
        return sums, float(np.mean(sums <= observed))
    e = E.reindex(plates.index).values; scaled = (n * E[K].values[None, :] / e[q]).sum(1)
    return sums, float(np.mean(sums <= observed)), scaled, float(np.mean(scaled <= observed))


def run(S, cat, plates, vcat, vplates, P_, E_pt, allp, crossed, rng):
    """Every step of the pooled test under shadow model S: VASCO's construction in S's frame, and ours with its trials."""
    vasco, *_ = pooled(vcat, vplates, Shadow(S.frame, edge="sharp", instant=True), weighted=False, exact=True)
    r_v = vasco["observed"] / vasco["expected"]                     # VASCO's shortfall, measured in this frame
    res, t, c, masks = pooled(cat, plates, S)
    sums, res["p_random"], scaled, p_pt = paste(c, plates, masks, res["observed"], lambda k: allp, rng, E=E_pt)
    # power: the share of no-effect trials that, scaled down by VASCO's shortfall in this frame, fall below the no-effect
    # trials' 5th percentile
    res["pointing"] = {"expected": round(float(sum(t.f_area[k] * E_pt[k] for k in masks)), 1), "null_mean": round(float(scaled.mean()), 1),
                       "p": round(p_pt, 3), "plates": {k: [round(float(t.W[k]), 1), round(float(E_pt[k]), 1)] for k in masks},
                       "parts": {k: {"shadow": round(float(t.w_in[k] / (E_pt[k] * t.f_area[k])), 3),
                                     "sunlit": round(float((t.W[k] - t.w_in[k]) / (E_pt[k] * (1 - t.f_area[k]))), 3) if t.f_area[k] < 1 else None} for k in masks},
                       "power_vasco": round(float(np.mean(r_v * scaled <= np.quantile(scaled, 0.05))), 3), "vasco_ratio": round(r_v, 3)}
    # random plates only from those pointed close to the anti-sun point at their own time, as the shadow plates are, with
    # each count scaled by pointing as above
    sun = get_sun(Time(plates.ut.values))
    dA = pd.Series(np.degrees(np.arccos(np.clip(np.sum(unit(plates.ra_c.values, plates.dec_c.values) * -unit(sun.ra.deg, sun.dec.deg), 1), -1, 1))), plates.index)
    nearp = allp[(dA < NEAR).reindex(allp).values]
    *_, p_near = paste(c, plates, masks, res["observed"], lambda k: nearp, rng, E=E_pt)
    res["near"] = {"max_deg": NEAR, "hosts": len(nearp), "p": round(p_near, 3),
                   "shadow_deg": [round(float(dA[list(masks)].min()), 1), round(float(dA[list(masks)].max()), 1)]}
    # random plates only from those near each shadow plate on the sky, so that a region of the catalogue the pointing model
    # misses (sparse or rich) is shared by the shadow plate and its hosts; counts scaled by pointing as above
    xyz = pd.DataFrame(unit(plates.ra_c.values, plates.dec_c.values), plates.index)
    res["local"] = {}
    for R in (12, 20):
        loc = {k: allp[np.degrees(np.arccos(np.clip(xyz.loc[allp].values @ xyz.loc[k].values, -1, 1))) < R] for k in masks}
        *_, p_loc = paste(c, plates, masks, res["observed"], lambda k: loc[k], rng, E=E_pt)
        res["local"][str(R)] = {"p": round(p_loc, 3), "fewest_hosts": int(min(len(v) for v in loc.values()))}
    # VASCO's check within 2 deg of the plate centre
    c2 = c[np.hypot(c.xi, c.eta) < 2]
    # (a turn or flip keeps a transient's distance from the centre, so the outlines laid on these count only within 2 deg)
    obs2 = float(c2.us[c2.plate.isin(list(masks))].sum())
    exp2 = c2.u.sum() * sum(m[R2].sum() / R2.sum() for m in masks.values()) / len(plates)
    _, p2 = paste(c2, plates, masks, obs2, lambda k: allp, rng)
    res["centre"] = {"observed": round(obs2, 1), "expected": round(float(exp2), 1), "p_random": round(p2, 3)}
    # each tenth of the ratings (VASCO-ML's edges), each transient counted once, as their deciles count them, with its own
    # pointing model. Every band is laid on the same random plates in each trial, so the trials keep how the bands move
    # together: in each trial, the lowest of the ten bands' chances gives how often some band looks as extreme as the one
    # picked out
    band = sample.decile(cat.w); edges = np.r_[0.0, sample.decile_edges(), 1.0]
    K = list(masks); qi_ = pd.Series(np.arange(len(plates)), plates.index)
    q = qi_[allp].values[rng.integers(0, len(allp), (NREP, len(K)))]; d = rng.integers(0, 8, (NREP, len(K)))
    res["bands"], P_null = [], []
    for k in range(10):
        r_k, t_k, c_k, _ = pooled(cat[band == k], plates, S, weighted=False)
        E_k = pd.Series(PP.pointing(P_.assign(n=t_k.N.reindex(P_.index)), use=~P_.index.isin(crossed))[0], P_.index)
        e_k = E_k.reindex(plates.index).values
        n_k = counts(c_k, plates, masks)[np.arange(len(K))[None, :], d, q]
        s_k = (n_k * E_k[K].values[None, :] / e_k[q]).sum(1)
        P_null.append(np.searchsorted(np.sort(s_k), s_k, side="right") / NREP)        # each trial's own chance in this band
        res["bands"].append({"band": k, "lo": round(float(edges[k]), 3), "hi": round(float(edges[k + 1]), 3),
                             "transients": int(t_k.N.sum()), "observed": int(r_k["observed"]), "expected": r_k["expected"],
                             "binom_z": r_k["binom_z"], "pointing_expected": round(float(sum(t_k.f_area[j] * E_k[j] for j in K)), 1),
                             "p_pointing": round(float(np.mean(s_k <= r_k["observed"])), 4)})
    lowest = np.min(P_null, 0); top = res["bands"][-1]
    res["top_tenth"] = {**top, "cut": top["lo"], "p_any_band": round(float(np.mean(lowest <= top["p_pointing"])), 3),
                        "lowest_band": int(np.argmin([b_["p_pointing"] for b_ in res["bands"]]))}
    # the same test with a sharp edge at VASCO-ML's radius, over the exposure: inside the widest shadow, so the same hosts
    r_, t_, c_, m_ = pooled(cat, plates, S.but(edge="sharp", radius=SHADOW_DEG_ML))
    *_, p_ml = paste(c_, plates, m_, r_["observed"], lambda k: allp, rng, E=E_pt)
    res["radius_ml"] = {"deg": SHADOW_DEG_ML, "plates": r_["plates"], "observed": r_["observed"], "expected": r_["expected"], "p_pointing": round(p_ml, 3)}
    # the model's systematics, each tested as above (the same hosts and pointing model)
    res["systematics"] = {}
    for name, S2 in SH.variants(S).items():
        r2, t2, c2, m2 = pooled(cat, plates, S2)
        *_, p2 = paste(c2, plates, m2, r2["observed"], lambda k: allp, rng, E=E_pt)
        res["systematics"][name] = {"plates": r2["plates"], "observed": r2["observed"],
                                    "pointing_expected": round(float(sum(t2.f_area[k] * E_pt[k] for k in m2)), 1), "p_pointing": round(p2, 3)}
    res["p_random"] = round(res["p_random"], 3)
    np.save(ROOT / "data/shadow_paste_sums.npy", sums); np.save(ROOT / "data/shadow_paste_scaled.npy", scaled)
    return {"vasco": vasco, "ours": res}, t


if __name__ == "__main__":
    rng = np.random.default_rng(11)
    vcat, vplates = sample.vasco()                 # VASCO's 2025 construction: their plates, on archive dates, counted once
    cat, plates = sample.ours()                    # ours: every transient weighted by its rating
    # outlines go only on plates the shadow never reaches (at its widest, at any time of the exposure), so no placement can
    # hold a real effect; the same host plates, and the same pointing model, serve every variant
    crossed = {q for q, p in plates.iterrows() if SH.cell_dark(p, SH.wide("topo")).any()}
    allp = np.array([q for q in plates.index if q not in crossed])
    clear = pd.Series({q: SH.clearance(plates.loc[q], SH.wide("topo")) for q in allp})
    assert (clear > 0).all()
    out = {"hosts": {"plates": len(allp), "crossed": sorted(crossed), "closest_deg": round(float(clear.min()), 2), "within_2deg": int((clear < 2).sum())}}
    # each plate's summed ratings predicted by its pointing, fitted without the plates the shadow reaches
    P_ = PP.properties(cat, plates)
    E_pt = pd.Series(PP.pointing(P_, use=~P_.index.isin(crossed), y="nw")[0], P_.index)
    out["topo"], tab = run(SH.MAIN, cat, plates, vcat, vplates, P_, E_pt, allp, crossed, rng)
    print({k: v for k, v in out["topo"]["ours"].items() if k not in ("bands", "pointing")}, flush=True)
    # VASCO's construction with the shadow placed as they place it, from the centre of the Earth: their published result
    out["geo"] = {"vasco": pooled(vcat, vplates, Shadow("geo", edge="sharp", instant=True), weighted=False, exact=True)[0]}
    print("VASCO's construction, from the centre of the Earth:", out["geo"]["vasco"], flush=True)
    # each plate's share, transients and summed ratings in shadow
    tab[["N", "W", "f_area", "n_in", "w_in"]].to_parquet(ROOT / "data/shadow_plates.parquet")
    tmp = ROOT / "data/shadow_pooled.json.tmp"; tmp.write_text(json.dumps(out, indent=1))
    tmp.replace(ROOT / "data/shadow_pooled.json")                              # whole, so a reader never sees half of it
