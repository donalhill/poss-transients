"""The shadow test's null, p-value, and limit on glints, from replicas of the real measurement: each of the real
shadows (shadow_fit.py's plates) moved, as it lies on its plate, onto a host plate (one the shadow does not reach at
its own time) chosen at random among those holding a similar number of transients (0.7 to 1.4 times, plus 2) at a
similar declination (within DEC_MATCH deg), no host used twice in a replica. Each replica keeps the real geometry and
gives it another plate's real transients, with their clumping and edge pile-up, so the spread of the statistics over
replicas is their spread with no effect, for this geometry (the geometry is ancillary: it carries no information about
glints, so the test is conditioned on it). Two other choices of host are kept as checks: on number alone, and on the
pointing model's expectation (within a factor 1.25).
Glints are put in by adding them: on each host plate, a Poisson number of glints with mean g0 / (1 - g0) times its
summed weight, spread evenly over the plate, each kept unless its location is in shadow (kept with chance 1 - s), and
each weighted by a rating drawn from the host's transients. g0 is then the share of a plate's transients that are
glints where it is sunlit throughout; glints on every plate raise the survey's total (E1's N) by 1 / (1 - g0) too. The page's observable is O / E, the transients observed in shadow over those
expected (shadow_stats.py); the 95% upper limit on g0 is the smallest g0 at which replicas give an O / E at or above the
real one only 5% of the time, and the expected limit is the median over the replicas with no glints. The same for Z.
Also, for the same replicas: the ratio fit, the glint-share fit, the pooled counts E1 and E2, and the hybrid H
(shadow_stats.py),
their spread with no effect under each choice of host, and how often each detects glints at a share of POWER_G0 (beyond
the 95th percentile of its own spread with no effect).
python analysis/replicas.py [REPLICAS] -> $POSS_SAMPLE/replicas.json, replica_draws.npz"""
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402
import shadow_stats as ST                                                    # noqa: E402
from shadow import HALF                                                      # noqa: E402
from shadow_fit import acceptance, cell_of, crossed_plates, real_plates       # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
M = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
DEC_MATCH, POWER_G0 = 15.0, 0.3
G0 = np.round(np.arange(0, 0.7001, 0.025), 3)                                # glint shares for the limit
SEED = int(os.environ.get("POSS_REPLICA_SEED", "20261026"))


def subpoint(xi, eta):
    """Each location's nearest sub-point of the grid (cell * 9 + k, SH.SUB's order)."""
    c = cell_of(xi, eta)
    ix = np.clip(np.rint((xi - SH.GX[c]) / (SH.CELL / 3)), -1, 1).astype(int)
    iy = np.clip(np.rint((eta - SH.GY[c]) / (SH.CELL / 3)), -1, 1).astype(int)
    return c * 9 + (ix + 1) * 3 + (iy + 1)


def setup():
    cat, pl = sample.load()
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)]
    A, qmax = acceptance(pl)
    P, _, _ = real_plates(cat, pl, aside=True)
    mu, y, alpha = ST.pointing(pl, cat)
    own = crossed_plates()
    hosts = [q for q in pl.index if q not in own and (cat.plate == q).any()]
    G = {q: g for q, g in cat.groupby("plate")}
    H = {q: (subpoint(G[q].xi.values, G[q].eta.values), G[q].w.values.astype(float)) for q in hosts}
    return cat, pl, A, qmax, P, mu, y, alpha, hosts, H


def candidates(P, pl, mu, hosts, H, match):
    out = []
    for p in P:
        q = p["plate"]; n0 = p["n"]
        if match == "count_dec":
            c = [h for h in hosts if 0.7 * n0 <= len(H[h][0]) <= 1.4 * n0 + 2 and abs(pl.dec_c[h] - pl.dec_c[q]) < DEC_MATCH]
        elif match == "count":
            c = [h for h in hosts if 0.7 * n0 <= len(H[h][0]) <= 1.4 * n0 + 2]
        else:
            c = [h for h in hosts if 0.8 * mu[q] <= mu[h] <= 1.25 * mu[q]]
        out.append(c)
    return out


def replica(P, cand, H, A9, mu, rng, g0s, fits=False, qmax=None, totals=None):
    """One replica: hosts drawn for each real shadow, then the statistics at each glint share in g0s (the same hosts and
    random numbers throughout): Z and the fits over the plates the test uses, the pooled counts over every plate the
    shadow reaches (with totals)."""
    used, hosts = set(), []
    for c in cand:
        free = [h for h in c if h not in used]
        h = free[rng.integers(len(free))]; used.add(h); hosts.append(h)
    seed = int(rng.integers(1 << 62)); res = []
    for g0 in g0s:
        r = np.random.default_rng(seed); Q = []
        for p, h in zip(P, hosts):
            sp, u = H[h]; s = p["s_sub"][sp]; a = A9[sp]
            if g0 > 0:
                k = r.poisson(g0 / (1 - g0) * u.sum() / u.mean()); idx = r.integers(0, A9.size, k)
                keep = r.random(k) >= p["s_sub"][idx]; idx = idx[keep]
                s = np.r_[s, p["s_sub"][idx]]; a = np.r_[a, A9[idx]]; u = np.r_[u, u[r.integers(0, len(u), len(idx))]]
            Q.append({**{k_: p[k_] for k_ in ("S", "F", "Es2", "QCAP", "aside")}, "s": s, "a": a, "u": u, "mu": mu[h] / (1 - g0)})
        st = ST.stats([x for x in Q if not x["aside"]], fits=fits, qmax=qmax)
        if totals is not None:                                               # glints on every plate raise the survey's total too
            tg = (totals[0] / (1 - g0), totals[1])
            st.update({k_: v for k_, v in ST.stats(Q, totals=tg).items() if k_ in ("E1", "E2")})
            if fits:
                st["H"] = ST.fit_hybrid([x for x in Q if not x["aside"]], Q, ALPHA)
        res.append(st)
    return res


def q5(x):
    return {k: round(float(v), 4) for k, v in zip(("q05", "q16", "median", "q84", "q95"), np.nanquantile(x, [0.05, 0.16, 0.5, 0.84, 0.95]))}


CHUNK = 500


def _null(args):
    """A chunk of replicas with no glints and at POWER_G0, every statistic (run in a worker: the setup is inherited)."""
    match, c = args
    rng = np.random.default_rng([SEED, c]); cand = CAND[match]
    return [replica(P, cand, H, A9, mu, rng, (0.0, POWER_G0), fits=True, qmax=qmax, totals=totals) for _ in range(CHUNK)]


def _at(args):
    """A chunk of replicas at one glint share, O / E alone (the limit's own distribution, for the chart)."""
    match, c, g0 = args
    rng = np.random.default_rng([SEED + 2, c]); cand = CAND[match]
    return [replica(P, cand, H, A9, mu, rng, (g0,))[0]["OE"] for _ in range(CHUNK)]


def _lim(args):
    """A chunk of replicas over the glint shares G0, O / E and Z alone."""
    match, c = args
    rng = np.random.default_rng([SEED + 1, c]); cand = CAND[match]
    return [[(t["OE"], t["Z"]) for t in replica(P, cand, H, A9, mu, rng, G0)] for _ in range(CHUNK)]


if __name__ == "__main__":
    import multiprocessing as mp
    cat, pl, A, qmax, P, mu, y, alpha, hosts, H = setup()
    A9 = np.repeat(A, 9)
    for p in P:
        p["s_sub"] = p["sub"].ravel()
    totals = (float(y.sum()), len(pl))
    for p in P:
        p["mu"] = float(mu[p["plate"]])
    real = ST.stats([p for p in P if not p["aside"]], fits=True, qmax=qmax)
    real.update({k: v for k, v in ST.stats(P, totals=totals).items() if k in ("E1", "E2")})
    ALPHA = alpha; real["H"] = ST.fit_hybrid([p for p in P if not p["aside"]], P, ALPHA)
    out = {"replicas": M, "real": {k: round(float(v), 4) for k, v in real.items()}, "hosts": len(hosts), "dec_match": DEC_MATCH,
           "matches": {}, "power_g0": POWER_G0}
    MATCHES = ("count_dec", "count", "pointing")
    CAND = {m_: candidates(P, pl, mu, hosts, H, m_) for m_ in MATCHES}
    nch, ML = max(M // CHUNK, 1), max(M // 4, 2000); nl = max(ML // CHUNK, 1)
    with mp.get_context("fork").Pool(int(os.environ.get("POSS_REPLICA_WORKERS", "10"))) as pool:
        NULL = {m_: sum(pool.map(_null, [(m_, c) for c in range(nch)]), []) for m_ in MATCHES}
        LIM = {m_: np.array(sum(pool.map(_lim, [(m_, c) for c in range(nl)]), [])) for m_ in MATCHES}
    draws = {}
    for match in MATCHES:
        cand = CAND[match]; two = NULL[match]; OEg, Zg = LIM[match][:, :, 0], LIM[match][:, :, 1]; Mm = len(two)
        names = ["Z", "OE", "r", "g", "E1", "E2", "H"]
        null = {k: np.array([t[0][k] for t in two]) for k in names}; alt = {k: np.array([t[1][k] for t in two]) for k in names}
        short = {"Z": 1, "OE": -1, "r": -1, "g": 1, "E1": -1, "E2": -1, "H": -1}   # the direction of a shortfall
        rec = {"candidates": [len(c) for c in cand], "null": {k: q5(v) for k, v in null.items()},
               "p": {k: round(float((np.sum(short[k] * null[k] >= short[k] * real[k]) + 1) / (Mm + 1)), 4) for k in names},
               "power": {k: round(float(np.mean(short[k] * alt[k] > np.nanquantile(short[k] * null[k], 0.95))), 4) for k in names}}
        def limit(b):
            b = np.minimum.accumulate(b)                                     # monotone in g0
            return float(np.interp(0.05, b[::-1], G0[::-1])) if b[-1] < 0.05 else None
        as_small = (OEg >= real["OE"]).mean(0)                               # replicas with a shortfall no larger than the real one
        lims = [limit((OEg >= v).mean(0)) for v in OEg[:, 0][:2000]]
        zl = [limit((Zg <= v).mean(0)) for v in Zg[:, 0][:2000]]
        rec.update(limit95=limit(as_small), as_small_by_g0=dict(zip([str(g) for g in G0], np.round(as_small, 4).tolist())),
                   expected_limit=q5([x for x in lims if x is not None]), no_limit_share=round(float(np.mean([x is None for x in lims])), 4),
                   z_limit95=limit((Zg <= real["Z"]).mean(0)), z_expected_limit=q5([x for x in zl if x is not None]))
        out["matches"][match] = rec
        print(match, json.dumps({k: v for k, v in rec.items() if k != "as_small_by_g0"}), flush=True)
        if match == "count_dec":
            gl = rec["limit95"] or G0[-1]; k = int(np.argmin(np.abs(G0 - gl)))
            draws = {"null_oe": null["OE"], "null_z": null["Z"], "limit_oe": OEg[:, k], "limit_g0": G0[k], "oe_by_g0": OEg, "z_by_g0": Zg, "g0": G0}
    gl = out["matches"]["count_dec"]["limit95"]
    if gl is not None:                                                       # the replicas at the limit itself, for the chart
        with mp.get_context("fork").Pool(int(os.environ.get("POSS_REPLICA_WORKERS", "10"))) as pool:
            draws.update(limit_oe=np.array(sum(pool.map(_at, [("count_dec", c, gl) for c in range(nl)]), [])), limit_g0=gl)
    np.savez(D / "replica_draws.npz", **draws)
    (D / "replicas.json").write_text(json.dumps(out, indent=1))
