"""Earth's shadow, measured plate by plate on VASCO's catalogue.

Hypothesis: on each plate, transients in the shadowed part are r times as dense as in the sunlit part (r = 1: no effect).
Fit: given a plate's n transients and the shadowed share f of its transients with no effect (its shadowed area, each
part weighted by the average density of transients there on plates the shadow never crosses), the count in shadow is binomial,
O ~ Bin(n, p), p = r f / (r f + 1 - f). Conditioning on n removes how rich or poor each plate is. One common r is fitted
over the plates partly in shadow by maximum likelihood; each plate's own r_i is shown with its 68% counting range
(Clopper-Pearson; an upper bound when none lie in shadow).

The binomial treats transients as spread evenly and independently within a plate; they are not (clumps, busy edges,
holes). So the fit is calibrated with toys:
  - from the fitted model (parametric): checks that the fit returns what went in, with the coverage it claims;
  - from the data: each real shadow footprint (its exact shape and place on the plate) laid on another catalogue plate the
    shadow never crosses, with a similar number of transients; the toy's counts are that plate's real transients inside
    and outside the footprint. An effect r < 1 is injected by removing transients inside at random, r > 1 by removing
    sunlit ones. Fitted exactly like the data: the no-effect distribution, and a Neyman interval on r.
Cross-checks: host plates chosen differently (any plate; the nearest on the sky); the shadow moved at random over its
own plate with the transients fixed (shadow-sized discs, similar coverage); each plate left out in turn.
Geocentric shadow (8.5 deg about the anti-sun point at mid-exposure); every catalogue plate on its corrected date
(sample.ours); the central 6 x 6 deg of each plate.
python analysis/shadow_fit.py -> data/shadow_fit.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
from shadow import (HALF, ROOT, SHADOW_DEG, SUB, antisun_axis, cell_counts, deproject, shadow_share,   # noqa: E402
                    unit, vasco_plate)

rng = np.random.default_rng(11)
T, TBELT = 20000, 5000
GRID = np.exp(np.linspace(np.log(0.05), np.log(20), 1201))
RBELT = np.round(np.arange(0.3, 2.501, 0.025), 3)


def counting_range(k, n, f):
    """68% Clopper-Pearson range on O/n, as a ratio r; None for an open upper end."""
    to_r = lambda x: float(x * (1 - f) / ((1 - x) * f))
    lo = to_r(beta.ppf(0.16, k, n - k + 1)) if k > 0 else 0.0
    hi = to_r(beta.ppf(0.84, k + 1, n - k)) if k < n else None
    return lo, hi


def make_fit(f):
    """Grid maximum likelihood of the common r for plates with shadowed shares f (rows of O, n: one experiment each)."""
    pg = GRID[:, None] * f[None, :] / (GRID[:, None] * f[None, :] + 1 - f[None, :])
    lp, lq = np.log(pg), np.log(1 - pg)
    def fit(O, n, interval=False):
        L = O @ lp.T + (n - O) @ lq.T
        best = GRID[L.argmax(1)]
        if not interval:
            return best
        ok = 2 * (L.max(1, keepdims=True) - L) <= 1
        return best, np.where(ok, GRID, np.inf).min(1), np.where(ok, GRID, -np.inf).max(1)
    return fit


def toys_from_hosts(P, pick, r_true, size):
    """Data-driven toys: for each plate, a host drawn by pick(p) (indices into the host list) per toy."""
    O = np.zeros((size, len(P))); n = np.zeros((size, len(P)))
    for j, p in enumerate(P):
        h = rng.choice(pick(p), size); o = p["hO"][h]; s = p["hN"][h] - o
        if r_true < 1:
            o = rng.binomial(o, r_true)
        elif r_true > 1:
            s = rng.binomial(s, 1 / r_true)
        O[:, j] = o; n[:, j] = o + s
    return O, n


def summary(d):
    q = np.quantile(d, [0.05, 0.16, 0.5, 0.84, 0.95])
    return {k: round(float(v), 3) for k, v in zip(("q05", "q16", "median", "q84", "q95"), q)}


def neyman(P, pick, fit, r_obs):
    """Central Neyman intervals (68%, 90%) on the true r from data-driven toys at each injected r."""
    rows = []
    for rt in RBELT:
        d = fit(*toys_from_hosts(P, pick, rt, TBELT))
        rows.append([rt, *np.quantile(d, [0.05, 0.16, 0.5, 0.84, 0.95])])
    B = np.array(rows)
    in68 = B[(B[:, 2] <= r_obs) & (B[:, 4] >= r_obs), 0]; in90 = B[(B[:, 1] <= r_obs) & (B[:, 5] >= r_obs), 0]
    return B, (float(in68.min()), float(in68.max())), (float(in90.min()), float(in90.max()))


def placements(P_, want=6000, batch=512):
    """Shadow-sized discs at random positions on one plate (shares of each 0.1 deg cell), cached in data/placements/."""
    cache = ROOT / f"data/placements/placements_vasco_{P_['plate']}.npz"
    if cache.exists():
        d = np.load(cache); return d["S"].astype(np.float32) / 255, d["f"]
    g = np.random.default_rng(int(P_["plate"][2:]))
    c = P_["sub"].reshape(-1, len(SUB), 3).mean(1).mean(0); c /= np.linalg.norm(c)
    e1 = np.cross(c, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(c, e1)
    cosr = np.cos(np.radians(SHADOW_DEG)); keep, fs = [], []
    for _ in range(200):
        th, d = np.radians(g.uniform(0, 360, batch)), np.radians(g.uniform(0, 13, batch))
        ax = np.cos(d)[:, None] * c + np.sin(d)[:, None] * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
        S = ((P_["sub"] @ ax.T) > cosr).reshape(-1, len(SUB), batch).mean(1).T.astype(np.float32)
        f = S.mean(1); ok = (f > 0.02) & (f < 0.97)
        keep.append(np.round(S[ok] * 255).astype(np.uint8)); fs.append(f[ok])
        if sum(len(k) for k in keep) >= want:
            break
    S = np.concatenate(keep)[:want]; f = np.concatenate(fs)[:want]
    cache.parent.mkdir(exist_ok=True); np.savez_compressed(cache, S=S, f=f)
    return S.astype(np.float32) / 255, f


if __name__ == "__main__":
    cat, pl = sample.ours()
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)]
    shadow9 = sorted(pd.read_parquet(ROOT / "data/shadow_plates.parquet").query("f_area > 0").index)
    hosts = {q: gq for q, gq in cat[~cat.plate.isin(shadow9)].groupby("plate")}
    # acceptance: the average density of transients over the face of a plate (each plate other than the shadow plates
    # counted equally); a cell's shadowed share is weighted by it, so edges and corners count for what they hold
    ACC = np.mean([cell_counts(gq.xi.values, gq.eta.values) / len(gq) for gq in hosts.values()], 0)
    ACC = ACC / ACC.mean()
    share = lambda sh: float(np.sum(sh * ACC) / ACC.sum())
    H = list(hosts); Hn = np.array([len(hosts[h]) for h in H])
    Hxyz = np.array([[np.cos(np.radians(pl.loc[h].dec_c)) * np.cos(np.radians(pl.loc[h].ra_c)),
                      np.cos(np.radians(pl.loc[h].dec_c)) * np.sin(np.radians(pl.loc[h].ra_c)), np.sin(np.radians(pl.loc[h].dec_c))] for h in H])

    P, whole = [], []
    for q in shadow9:
        p = pl.loc[q]; V = vasco_plate(q, p, cat); s = shadow_share(V["sub"], p.ut); area = float(s.mean()); f = share(s)
        if not 0.02 < area < 0.97:
            whole.append(q); continue
        ax = np.asarray(antisun_axis(Time(pd.Timestamp(p.ut).to_datetime64()))).reshape(3)
        inside = lambda xi, eta: (unit(*deproject(xi, eta, p.ra_c, p.dec_c)) @ ax) > np.cos(np.radians(SHADOW_DEG))
        c = cat[cat.plate == q]; n = len(c); O = int(inside(c.xi.values, c.eta.values).sum())
        cen = Hxyz @ np.asarray(unit(p.ra_c, p.dec_c)).reshape(3)
        P.append(dict(plate=q, f=f, area=area, n=n, O=O, s=s, V=V,
                      hO=np.array([inside(hosts[h].xi.values, hosts[h].eta.values).sum() for h in H]), hN=Hn,
                      similar=np.flatnonzero((Hn > n / 1.5) & (Hn < n * 1.5)), near=np.argsort(-cen)[:12]))
    f = np.array([p["f"] for p in P]); nR = np.array([p["n"] for p in P], float); OR = np.array([p["O"] for p in P], float)
    fit = make_fit(f)
    r_hat, c_lo, c_hi = (float(v[0]) for v in fit(OR[None], nR[None], interval=True))
    out = {"plates": [], "whole": whole, "fit": {"r": round(r_hat, 3), "count_lo": round(c_lo, 3), "count_hi": round(c_hi, 3),
                                                 "O": int(OR.sum()), "n": int(nR.sum())}}
    for p in P:
        lo, hi = counting_range(p["O"], p["n"], p["f"]); E = (p["n"] - p["O"]) * p["f"] / (1 - p["f"])
        out["plates"].append({"plate": p["plate"], "shadow": round(p["area"], 3), "weighted": round(p["f"], 3), "n": p["n"], "O": p["O"], "E": round(E, 1),
                              "r": round(p["O"] / E, 3) if E > 0 else None, "lo": round(lo, 3), "hi": None if hi is None else round(hi, 3)})
    print(pd.DataFrame(out["plates"]).to_string(index=False))
    print(f"common r {r_hat:.3f}, counting-only 68% {c_lo:.3f}-{c_hi:.3f}")

    # 1. the fit on its own model: generate from the binomial at the fitted r, fit back
    Op = rng.binomial(nR.astype(int), r_hat * f / (r_hat * f + 1 - f), (T, len(P))).astype(float)
    rb, lo_b, hi_b = fit(Op, np.broadcast_to(nR, Op.shape), interval=True)
    out["closure"] = {"median": round(float(np.median(rb)), 3), "coverage68": round(float(np.mean((lo_b <= r_hat) & (hi_b >= r_hat))), 3),
                      **{k: v for k, v in summary(rb).items() if k in ("q16", "q84")}}
    print("closure on the binomial model:", out["closure"])

    # 2. toys from the data: footprints on other plates with a similar number of transients
    pick = lambda p: p["similar"]
    d0 = fit(*toys_from_hosts(P, pick, 1.0, T))
    B, i68, i90 = neyman(P, pick, fit, r_hat)
    out["toys"] = {"null": summary(d0), "p_as_high": round(float(np.mean(d0 >= r_hat)), 3),
                   "hist_null": np.histogram(d0, np.linspace(0, 3, 61))[0].tolist(),
                   "interval68": [round(v, 3) for v in i68], "interval90": [round(v, 3) for v in i90],
                   "belt": [[round(float(x), 3) for x in row] for row in B]}
    print("toys from the data, no effect:", out["toys"]["null"], "share as high as the data", out["toys"]["p_as_high"])
    print("Neyman interval on r: 68%", i68, "90%", i90)
    # a shortfall as large as VASCO's (their pooled count over expected, on their construction): how often the fit would
    # come back as high as it does on the real plates
    V_ = json.loads((ROOT / "data/shadow_pooled.json").read_text())["vasco_geo"]
    r_v = round(V_["observed"] / V_["expected"], 2)
    out["injected"] = {"r_true": r_v, "p_as_high": round(float(np.mean(fit(*toys_from_hosts(P, pick, r_v, T)) >= r_hat)), 3)}
    print("injected", out["injected"])

    # 3. cross-checks: other host choices; the shadow moved over its own plate; each plate left out
    out["checks"] = {}
    for key, pk in (("any_plate", lambda p: np.arange(len(H))), ("near_on_sky", lambda p: p["near"])):
        d = fit(*toys_from_hosts(P, pk, 1.0, T)); _, a68, a90 = neyman(P, pk, fit, r_hat)
        out["checks"][key] = {"null": summary(d), "p_as_high": round(float(np.mean(d >= r_hat)), 3),
                              "interval68": [round(v, 3) for v in a68], "interval90": [round(v, 3) for v in a90]}
        print(f"hosts {key}:", out["checks"][key])
    # same plate, transients fixed: shadow-sized discs of similar coverage at random positions; each toy fitted with
    # its own shadowed shares
    Os, ns, fs = [], [], []
    for p in P:
        S, fa = placements(p["V"]); fa = S @ ACC / ACC.sum(); m = (fa > p["f"] / 1.5) & (fa < min(0.97, p["f"] * 1.5)); S, fa = S[m], fa[m]
        i = rng.integers(0, len(S), T); cnt = cell_counts(cat[cat.plate == p["plate"]].xi.values, cat[cat.plate == p["plate"]].eta.values)
        Os.append(S[i] @ cnt); ns.append(np.full(T, cnt.sum())); fs.append(fa[i])
    Os, ns, fs = np.array(Os).T, np.array(ns).T, np.array(fs).T
    best = np.empty(T)
    for k0 in range(0, T, 2000):
        o, n_, ff = Os[k0:k0 + 2000], ns[k0:k0 + 2000], fs[k0:k0 + 2000]
        pg = GRID[None, None, :] * ff[..., None] / (GRID[None, None, :] * ff[..., None] + 1 - ff[..., None])
        L = (o[..., None] * np.log(pg) + (n_ - o)[..., None] * np.log(1 - pg)).sum(1)
        best[k0:k0 + 2000] = GRID[L.argmax(1)]
    out["checks"]["same_plate"] = {"null": summary(best), "p_as_high": round(float(np.mean(best >= r_hat)), 3)}
    print("shadow moved over its own plate:", out["checks"]["same_plate"])
    out["checks"]["leave_one_out"] = {}
    for j, p in enumerate(P):
        k = np.arange(len(P)) != j
        out["checks"]["leave_one_out"][p["plate"]] = round(float(make_fit(f[k])(OR[None, k], nR[None, k])[0]), 3)
    print("each plate left out:", out["checks"]["leave_one_out"])
    (ROOT / "data/shadow_fit.json").write_text(json.dumps(out))
