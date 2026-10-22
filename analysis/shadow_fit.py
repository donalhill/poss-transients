"""Earth's shadow, measured plate by plate on VASCO's catalogue.

Hypothesis: on each plate, transients in the shadowed part are r times as dense as in the sunlit part (r = 1: no effect).
Each transient counts by VASCO-ML's rating of how likely it is to be real (sample.py), in every share, density and toy
below; the likelihood uses each plate's effective number of transients (see measure). Also the fit in each tenth of the
ratings (VASCO-ML's decile edges, sample.decile), each transient counted once, as their deciles count them; and every
band at once, with a free signal purity in each (joint_bands).
Fit: given a plate's n transients and the shadowed share f of its transients with no effect (its shadowed area, each
part weighted by the average density of transients there on plates the shadow never crosses), the count in shadow is binomial,
O ~ Bin(n, p), p = r f / (r f + 1 - f). Conditioning on n removes how rich or poor each plate is. One common r is fitted
over the plates partly in shadow by maximum likelihood; each plate's own r_i is shown with its 68% counting range
(Clopper-Pearson; an upper bound when none lie in shadow).

The binomial treats transients as spread evenly and independently within a plate; they are not (clumps, busy edges,
holes). So the fit is calibrated with toys:
  - from the fitted model (parametric): checks that the fit returns what went in, with the coverage it claims;
  - from the data: each real shadow footprint (its exact shape and place on the plate) laid on another catalogue plate the
    shadow never crosses (geocentric or topocentric), with a similar number of transients; the toy's counts are that plate's real transients inside
    and outside the footprint. An effect r < 1 is injected by removing transients inside at random, r > 1 by removing
    sunlit ones. Fitted exactly like the data: the no-effect distribution, and a Neyman interval on r.
Cross-checks: host plates chosen differently (any plate; the nearest on the sky); the shadow moved at random over its
own plate with the transients fixed (shadow-sized discs, similar coverage); each plate left out in turn.
Geocentric shadow (8.69 deg about the anti-sun point at mid-exposure); every catalogue plate on its corrected date
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
RBELT = np.round(np.arange(0.1, 2.501, 0.025), 3)


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
    """Data-driven toys: for each plate, a host drawn by pick(p) (indices into the host list) per toy. An effect is put in
    by removing transients at random, inside (r < 1) or outside (r > 1), each taking its weight with it (the host's
    weights inside or outside scaled by the share of transients kept). Returned as effective counts (the host's scale k)."""
    O = np.zeros((size, len(P))); n = np.zeros((size, len(P)))
    for j, p in enumerate(P):
        h = rng.choice(pick(p), size); o = p["hO"][h]; s = p["hN"][h] - o; oc = p["hOc"][h]; sc = p["hNc"][h] - oc
        if r_true < 1:
            o = o * rng.binomial(oc, r_true) / np.maximum(oc, 1)
        elif r_true > 1:
            s = s * rng.binomial(sc, 1 / r_true) / np.maximum(sc, 1)
        O[:, j] = o * p["hK"][h]; n[:, j] = (o + s) * p["hK"][h]
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
    # each end where the toys' quantile crosses the observed r, interpolated between the two injected values either side
    # (the lower end at the first crossing of the upper quantile, the upper end at the last crossing of the lower one)

    def edge(q, first):
        s = np.flatnonzero(np.diff(np.sign(B[:, q] - r_obs)) != 0)
        if not len(s):                                                       # no crossing: the end of the injected range
            return float(B[0, 0] if first else B[-1, 0])
        i = s[0] if first else s[-1]
        return float(np.interp(r_obs, B[i:i + 2, q], B[i:i + 2, 0])) if B[i, q] < B[i + 1, q] else float(B[i, 0])
    return B, (edge(4, True), edge(2, False)), (edge(5, True), edge(1, False))


def placements(P_, want=6000, batch=512):
    """Shadow-sized discs at random positions on one plate (shares of each 0.1 deg cell), cached in data/placements/."""
    cache = ROOT / f"data/placements/placements_vasco_{P_['plate']}_{SHADOW_DEG}.npz"
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


def measure(cat, pl, weighted, full):
    """The fit on the plates the shadow partly covers, each transient counted by its rating (weighted) or once. Shares,
    the density map and every toy use the same weights. The likelihood takes each plate's effective number of transients,
    k = sum(w) / sum(w^2) times its summed ratings, so counting statistics are those of the ratings; with no weights k = 1.
    full: also the cross-checks (other hosts, the shadow moved over its own plate, each plate left out)."""
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].assign(u=lambda d: d.w if weighted else 1.0)
    SPL = pd.read_parquet(ROOT / "data/shadow_plates.parquet")
    shadow9 = sorted(SPL.query("f_area > 0").index)
    # hosts: plates neither shadow (geocentric or topocentric) crosses, so no toy can hold a real effect
    hosts = {q: gq for q, gq in cat[~cat.plate.isin(SPL.query("f_area > 0 or f_area_topo > 0").index)].groupby("plate")}
    # acceptance: the average density of transients (by weight) over the face of a plate, each plate other than the shadow
    # plates counted equally; a cell's shadowed share is weighted by it, so edges and corners count for what they hold
    ACC = np.mean([cell_counts(gq.xi.values, gq.eta.values, gq.u.values) / gq.u.sum() for gq in hosts.values()], 0)
    ACC = ACC / ACC.mean()
    share = lambda sh: float(np.sum(sh * ACC) / ACC.sum())
    H = list(hosts); Hc = np.array([len(hosts[h]) for h in H]); Hn = np.array([hosts[h].u.sum() for h in H])
    HK = np.array([hosts[h].u.sum() / (hosts[h].u ** 2).sum() for h in H])
    Hxyz = np.array([[np.cos(np.radians(pl.loc[h].dec_c)) * np.cos(np.radians(pl.loc[h].ra_c)),
                      np.cos(np.radians(pl.loc[h].dec_c)) * np.sin(np.radians(pl.loc[h].ra_c)), np.sin(np.radians(pl.loc[h].dec_c))] for h in H])

    P, whole = [], []
    for q in shadow9:
        p = pl.loc[q]; V = vasco_plate(q, p, cat); s = shadow_share(V["sub"], p.ut); area = float(s.mean()); f = share(s)
        if not 0.02 < area < 0.97:
            whole.append(q); continue
        c = cat[cat.plate == q]
        if c.empty:
            continue
        ax = np.asarray(antisun_axis(Time(pd.Timestamp(p.ut).to_datetime64()))).reshape(3)
        inside = lambda xi, eta: (unit(*deproject(xi, eta, p.ra_c, p.dec_c)) @ ax) > np.cos(np.radians(SHADOW_DEG))
        m = inside(c.xi.values, c.eta.values); hin = [inside(hosts[h].xi.values, hosts[h].eta.values) for h in H]
        cen = Hxyz @ np.asarray(unit(p.ra_c, p.dec_c)).reshape(3)
        n_w = c.u.sum(); sim = np.flatnonzero((Hn > n_w / 1.5) & (Hn < n_w * 1.5))
        P.append(dict(plate=q, f=f, area=area, n=len(c), O=int(m.sum()), nw=float(n_w), Ow=float(c.u[m].sum()),
                      k=float(n_w / (c.u ** 2).sum()), s=s, V=V,
                      hO=np.array([hosts[h].u.values[i].sum() for h, i in zip(H, hin)]), hN=Hn,
                      hOc=np.array([i.sum() for i in hin]), hNc=Hc, hK=HK,
                      similar=sim if len(sim) >= 20 else np.argsort(np.abs(Hn - n_w))[:20], near=np.argsort(-cen)[:12]))
    f = np.array([p["f"] for p in P]); K_ = np.array([p["k"] for p in P])
    nR = np.array([p["nw"] for p in P]) * K_; OR = np.array([p["Ow"] for p in P]) * K_
    fit = make_fit(f)
    r_hat, c_lo, c_hi = (float(v[0]) for v in fit(OR[None], nR[None], interval=True))
    out = {"weighted": weighted, "plates": [], "whole": whole,
           "fit": {"r": round(r_hat, 3), "count_lo": round(c_lo, 3), "count_hi": round(c_hi, 3), "O": round(float(sum(p["Ow"] for p in P)), 1),
                   "n": round(float(sum(p["nw"] for p in P)), 1), "transients": int(sum(p["n"] for p in P))}}
    for p in P:
        lo, hi = counting_range(p["Ow"] * p["k"], p["nw"] * p["k"], p["f"]); E = (p["nw"] - p["Ow"]) * p["f"] / (1 - p["f"])
        out["plates"].append({"plate": p["plate"], "shadow": round(p["area"], 3), "weighted": round(p["f"], 3), "n": p["n"], "O": p["O"],
                              "nw": round(p["nw"], 1), "Ow": round(p["Ow"], 1), "E": round(E, 1),
                              "r": round(p["Ow"] / E, 3) if E > 0 else None, "lo": round(lo, 3), "hi": None if hi is None else round(hi, 3)})
    print(pd.DataFrame(out["plates"]).to_string(index=False))
    print(f"common r {r_hat:.3f}, counting-only 68% {c_lo:.3f}-{c_hi:.3f}")

    # 1. the fit on its own model: generate from the binomial at the fitted r, fit back
    Op = rng.binomial(np.round(nR).astype(int), r_hat * f / (r_hat * f + 1 - f), (T, len(P))).astype(float)
    rb, lo_b, hi_b = fit(Op, np.broadcast_to(np.round(nR), Op.shape), interval=True)
    out["closure"] = {"median": round(float(np.median(rb)), 3), "coverage68": round(float(np.mean((lo_b <= r_hat) & (hi_b >= r_hat))), 3),
                      **{k: v for k, v in summary(rb).items() if k in ("q16", "q84")}}
    print("closure on the binomial model:", out["closure"])

    # 2. toys from the data: footprints on other plates with a similar total
    pick = lambda p: p["similar"]
    d0 = fit(*toys_from_hosts(P, pick, 1.0, T))
    B, i68, i90 = neyman(P, pick, fit, r_hat)
    # a 95% upper limit on the share of transients that vanish in shadow (1 - r): the lowest true r whose toys reach the
    # fitted r at least 5% of the time, from the belt's 95th percentile; and the limit expected with no effect, from the
    # no-effect toys in place of the data
    hi95 = np.maximum.accumulate(B[:, 5])
    up = lambda r_: 1 - np.interp(r_, hi95, B[:, 0])
    lim = {"share95": round(float(up(r_hat)), 3), "expected": [round(float(np.quantile(up(d0), q)), 3) for q in (0.025, 0.16, 0.5, 0.84, 0.975)]}
    out["toys"] = {"trials": T, "null": summary(d0), "p_as_high": round(float(np.mean(d0 >= r_hat)), 3), "limit": lim,
                   "p_as_low": round(float(np.mean(d0 <= r_hat)), 3),
                   "hist_null": np.histogram(d0, np.linspace(0, 3, 61))[0].tolist(),
                   "interval68": [round(v, 3) for v in i68], "interval90": [round(v, 3) for v in i90],
                   "belt": [[round(float(x), 3) for x in row] for row in B]}
    print("toys from the data, no effect:", out["toys"]["null"], "share as high as the data", out["toys"]["p_as_high"])
    print("Neyman interval on r: 68%", i68, "90%", i90)
    # a shortfall as large as VASCO's (their pooled count over expected, on their construction): how often the fit would
    # come back as high as it does on the real plates
    V_ = json.loads((ROOT / "data/shadow_pooled.json").read_text())["vasco_geo"]
    r_v = round(V_["observed"] / V_["expected"], 2)
    d_v = fit(*toys_from_hosts(P, pick, r_v, T))
    out["injected"] = {"r_true": r_v, "p_as_high": round(float(np.mean(d_v >= r_hat)), 3),
                       "found": round(float(np.mean(d_v <= np.quantile(d0, 0.05))), 3),    # flagged at the 95% level
                       "hist": np.histogram(d_v, np.linspace(0, 3, 61))[0].tolist()}
    print("injected", out["injected"])
    if not full:
        return out

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
        c = cat[cat.plate == p["plate"]]
        i = rng.integers(0, len(S), T); cnt = cell_counts(c.xi.values, c.eta.values, c.u.values) * p["k"]
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
    return out


def joint_bands(cat, pl, band, size=T):
    """Signal purity in every band of the ratings at once. In band k a share pi_k of transients vanish in shadow (glints),
    the rest do not (flaws, missed stars), so its ratio is r_k = 1 - pi_k; pi_k is free in each band, as nothing says
    beforehand how purity varies with rating. Each band's count in shadow on each plate partly covered is binomial, with
    that band's own density map (ratings vary across a plate). The test of no signal in any band, fixed before looking:
    q = sum over bands of 2 (ln L(pi hat) - ln L(0)), with each pi_k held at 0 or above (signal only removes transients).
    Its null: in each trial every shadow outline is laid on one host plate for all ten bands at once (a plate the shadow
    never crosses, with a similar number of transients), so the bands keep how they move together."""
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].assign(k=band[(cat.xi.abs() < HALF).values & (cat.eta.abs() < HALF).values])
    SPL = pd.read_parquet(ROOT / "data/shadow_plates.parquet")
    hosts = {q: gq for q, gq in cat[~cat.plate.isin(SPL.query("f_area > 0 or f_area_topo > 0").index)].groupby("plate")}
    H = list(hosts); NB = 10
    ACC = []
    for k in range(NB):
        m_ = [cell_counts(g.xi.values[g.k.values == k], g.eta.values[g.k.values == k]) / (g.k.values == k).sum() for g in hosts.values() if (g.k.values == k).any()]
        a = np.mean(m_, 0); ACC.append(a / a.mean())
    Hn = np.array([len(hosts[h]) for h in H]); Hk = np.array([np.bincount(hosts[h].k.values, minlength=NB) for h in H])   # hosts x bands
    f, O, n, hO, sim = [], [], [], [], []
    for q in sorted(SPL.query("f_area > 0").index):
        p = pl.loc[q]; V = vasco_plate(q, p, cat); sh = shadow_share(V["sub"], p.ut)
        if not 0.02 < sh.mean() < 0.97:
            continue
        ax = np.asarray(antisun_axis(Time(pd.Timestamp(p.ut).to_datetime64()))).reshape(3)
        inside = lambda xi, eta: (unit(*deproject(xi, eta, p.ra_c, p.dec_c)) @ ax) > np.cos(np.radians(SHADOW_DEG))
        c = cat[cat.plate == q]; m = inside(c.xi.values, c.eta.values)
        f.append([float(np.sum(sh * a) / a.sum()) for a in ACC])
        O.append(np.bincount(c.k.values[m], minlength=NB)); n.append(np.bincount(c.k.values, minlength=NB))
        hO.append(np.array([np.bincount(hosts[h].k.values[inside(hosts[h].xi.values, hosts[h].eta.values)], minlength=NB) for h in H]))
        s_ = np.flatnonzero((Hn > len(c) / 1.5) & (Hn < len(c) * 1.5)); sim.append(s_ if len(s_) >= 20 else np.argsort(np.abs(Hn - len(c)))[:20])
    f, O, n = np.array(f), np.array(O, float), np.array(n, float)                     # plates x bands
    R = np.linspace(0.0, 1.0, 401)[1:]                                                  # r = 1 - purity, from 1 (none) down

    # benchmark D, VASCO's labels at face value: the signal spread over the bands as their labelled-real examples are, so
    # band k's purity is mu t_k / u_k (t: the labelled-real examples' share in each band; u: each band's share of the
    # transients on these plates), one strength mu for every band, up to the value at which a band's purity reaches 1
    lab = sample.labelled(); lr = lab[(lab.label == 1) & (lab.row >= 0)]
    t_ = np.bincount(sample.decile(lr.prob.values), minlength=NB) / len(lr); u_ = n.sum(0) / n.sum()
    MU = np.linspace(0, float(np.min(u_[t_ > 0] / t_[t_ > 0])), 201)
    RI = np.abs(R[None, None, :] - (1 - MU[:, None] * t_[None, :] / u_[None, :])[..., None]).argmin(2)          # mu x bands -> r index

    def q_stat(O_, n_):
        """O_, n_: trials x plates x bands. Per band, 2 (ln L at the best r <= 1 - ln L at r = 1); the best r; and benchmark
        D's 2 (ln L at the best mu - ln L at mu = 0), with the best mu."""
        pg = R[:, None, None] * f[None] / (R[:, None, None] * f[None] + 1 - f[None])                # r x plates x bands
        L = np.einsum("tpb,rpb->trb", O_, np.log(pg)) + np.einsum("tpb,rpb->trb", n_ - O_, np.log(1 - pg))
        LD = L[:, RI, np.arange(NB)[None, :]].sum(2)                                                 # trials x mu
        return 2 * (L.max(1) - L[:, -1, :]), R[L.argmax(1)], 2 * (LD.max(1) - LD[:, 0]), MU[LD.argmax(1)]
    qb, rb, qd, mud = q_stat(O[None], n[None]); q_obs = float(qb.sum())

    def trials(purity, size):
        """q on trials from host plates, with a share purity[k] of each band's transients inside the outline removed."""
        Q = []
        for t0 in range(0, size, 1000):
            t = min(1000, size - t0)
            h = np.stack([rng.choice(s_, t) for s_ in sim], 1)                                     # trials x plates
            Ot = np.stack([hO[j][h[:, j]] for j in range(len(sim))], 1).astype(float)
            nt = Hk[h].astype(float)
            gone = Ot - rng.binomial(Ot.astype(int), 1 - np.asarray(purity)[None, None, :])
            q_ = q_stat(Ot - gone, nt - gone); Q.append(np.c_[q_[0].sum(1), q_[2]])
        return np.concatenate(Q)
    QQ = trials(np.zeros(NB), size); Q, QD = QQ[:, 0], QQ[:, 1]; q95 = float(np.quantile(Q, 0.95))
    # a trend with rating, as a signal among genuine transients would give: the slope of log r on each band's mean rating
    # (r unconstrained), against the same slope on trials from host plates, every band on the same hosts
    Rw = np.exp(np.linspace(np.log(0.02), np.log(8), 1201))
    def best(O_, n_):
        pg = Rw[:, None, None] * f[None] / (Rw[:, None, None] * f[None] + 1 - f[None])
        L = np.einsum("tpb,rpb->trb", O_, np.log(pg)) + np.einsum("tpb,rpb->trb", n_ - O_, np.log(1 - pg))
        return Rw[L.argmax(1)]
    mr = np.array([cat.w[cat.k == k].mean() for k in range(NB)])
    slope = lambda r_: np.polyfit(mr, np.log(r_), 1)[0]
    s_obs = slope(best(O[None], n[None])[0]); S_ = []
    for t0 in range(0, size, 1000):
        t = min(1000, size - t0); h = np.stack([rng.choice(s_, t) for s_ in sim], 1)
        Ot = np.stack([hO[j][h[:, j]] for j in range(len(sim))], 1).astype(float)
        S_ += [slope(r_) for r_ in best(Ot, Hk[h].astype(float))]
    S_ = np.array(S_)
    # the test's reach: the top band's fitted purity in the top band alone, and VASCO's pooled shortfall in every band
    V_ = json.loads((ROOT / "data/shadow_pooled.json").read_text())["vasco_geo"]; v_ = round(1 - V_["observed"] / V_["expected"], 2)
    inj = {"top_band": [0.0] * (NB - 1) + [round(float(1 - rb[0][-1]), 2)], "every_band": [v_] * NB}
    reach = {k: {"purity": v, "found": round(float(np.mean(trials(v, 5000)[:, 0] >= q95)), 3)} for k, v in inj.items()}
    return {"q": round(q_obs, 2), "p": round(float((np.sum(Q >= q_obs) + 1) / (len(Q) + 1)), 3), "trials": int(len(Q)),
            "q_by_band": [round(float(v), 2) for v in qb[0]], "purity": [round(float(1 - v), 3) for v in rb[0]],
            "null_q": {"median": round(float(np.median(Q)), 2), "q95": round(q95, 2)}, "reach": reach,
            "labels_template": {"template": [round(float(v), 3) for v in t_], "mu": round(float(mud[0]), 3), "q": round(float(qd[0]), 2),
                                "p": round(float((np.sum(QD >= qd[0]) + 1) / (len(QD) + 1)), 3), "top_band_purity_at_mu": round(float(mud[0] * t_[-1] / u_[-1]), 3)},
            "trend": {"slope": round(float(s_obs), 3), "null_median": round(float(np.median(S_)), 3),
                      "null_90": [round(float(np.quantile(S_, 0.05)), 3), round(float(np.quantile(S_, 0.95)), 3)],
                      "p": round(float(np.mean(S_ <= s_obs)), 3), "mean_rating": [round(float(v), 3) for v in mr]}}


if __name__ == "__main__":
    cat, pl = sample.ours()
    out = measure(cat, pl, weighted=True, full=True)
    # the same fit in each tenth of the ratings (VASCO-ML's edges), each transient counted once, as their deciles count them
    edges = np.r_[0.0, sample.decile_edges(), 1.0]; band = sample.decile(cat.w)
    out["bands"] = []
    for k in range(10):
        print(f"\nrating band {k} ({edges[k]:.3f} to {edges[k + 1]:.3f})")
        m = measure(cat[band == k], pl, weighted=False, full=False)
        out["bands"].append({"band": k, "lo": round(float(edges[k]), 3), "hi": round(float(edges[k + 1]), 3), "fit": m["fit"],
                             "plates": [{q: p_[q] for q in ("plate", "n", "O", "E", "r")} for p_ in m["plates"]],
                             "toys": {q: m["toys"][q] for q in ("null", "p_as_low", "p_as_high", "interval68", "interval90", "limit")}})
    # every band at once, each with its own free purity
    out["bands_joint"] = joint_bands(cat, pl, band)
    print("all bands at once:", out["bands_joint"])
    (ROOT / "data/shadow_fit.json").write_text(json.dumps(out))
