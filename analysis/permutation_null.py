"""The shadow test's cross-check with shuffled plate times: the plate times permuted among the plates, so that real
shadows fall on other plates and in other locations while every plate keeps its own transients. Under no effect the
transients carry no memory of when their plate was taken, so a shuffle is as likely as the real assignment.
For each plate with transients and each of the times, the shadow's reach is tested as crossed.py does (the widest
shadow, anywhere on the central square) and, where it reaches, the plate enters as shadow_fit.py makes it, by the same
rules (the plates it leaves out, wholly in shadow or barely touched, kept for the pooled counts only). Each shuffle
gives the statistics of shadow_stats.py (Z and O / E, and the pooled counts E1 and E2). The shuffles
mix geometries unlike the real one (other numbers of plates, other shares in shadow), so the p-values are also given
over the shuffles whose counting variance of Z lies within 20% of the real one's; replicas.py conditions on the real
geometry exactly and gives the test's null. POSS_PERM_SEED fixes the permutations.
python analysis/permutation_null.py [PERMUTATIONS] -> $POSS_SAMPLE/permutation_null.json, permutation_draws.npz"""
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402
import shadow_stats as ST                                                    # noqa: E402
from shadow import HALF                                                      # noqa: E402
from shadow_fit import acceptance, plate, usable                             # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
B = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 10000
MAIN, WIDE = SH.MAIN, SH.wide("topo")
cat, pl = sample.load()
cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)]
G = {q: g for q, g in cat.groupby("plate")}
TIMES = pl[["ut", "exp"]].reset_index(drop=True)
PLATES = [q for q in pl.index if q in G]
A_MAP, _ = acceptance(pl)


def plate_pairs(q):
    """Every time at which the shadow reaches plate q: (time index, the plate's record, flagged "aside" where the test
    leaves it out)."""
    p0 = pl.loc[q]; out = []
    for j, (ut, ex) in enumerate(zip(TIMES.ut.values, TIMES.exp.values)):
        p = p0.copy(); p["ut"] = ut; p["exp"] = ex
        if SH.near_miss(p, WIDE) or not SH.dark(SH.GX, SH.GY, p, WIDE).any():     # reach: tested on cell centres
            continue
        sub = SH.dark((SH.GX[:, None] + SH.SUB[:, 0]).ravel(), (SH.GY[:, None] + SH.SUB[:, 1]).ravel(), p, MAIN).reshape(-1, len(SH.SUB))
        rule = usable(sub)
        if rule == "none":
            continue
        d = plate(G[q], p, sub, A_MAP); d.pop("a"); d["aside"] = rule != "ok"
        out.append((j, d))
    return q, out


def load_pairs():
    import hashlib
    import pickle
    cache = D / "permutation_pairs.pkl"                                      # the pairs depend only on the sample and its times
    key = (sorted(PLATES), [str(x) for x in TIMES.ut.values], (D / "plates.parquet").stat().st_mtime, (D / "transients.parquet").stat().st_mtime,
           hashlib.sha256(b"".join(Path(f).read_bytes() for f in (__file__, SH.__file__, ST.__file__, sys.modules["shadow_fit"].__file__))).hexdigest())
    if cache.exists():
        k0, T0 = pickle.loads(cache.read_bytes())
        if k0 == key:
            return T0
    with Pool(min(int(os.environ.get("POSS_PERM_WORKERS", "8")), os.cpu_count())) as pool:
        TAB = dict(pool.map(plate_pairs, PLATES))
    cache.write_bytes(pickle.dumps((key, TAB)))
    return TAB


def var_z(P):
    return float(sum(float((p["u"] ** 2).sum()) * (p["Es2"] - p["F"] ** 2) for p in P if not p["aside"]))


def both(P):
    """Z and O / E over the plates the test uses, the pooled counts over all."""
    st = ST.stats([p for p in P if not p["aside"]])
    st.update({k: v for k, v in ST.stats(P, totals=totals).items() if k in ("E1", "E2")})
    return st


if __name__ == "__main__":
    TAB = load_pairs()
    mu, y, alpha = ST.pointing(pl, cat); totals = (float(y.sum()), len(pl))
    for q in TAB:
        for (_, d) in TAB[q]:
            d["mu"] = float(mu[q])
    idx = {q: dict(TAB[q]) for q in PLATES}; pos = {q: i for i, q in enumerate(pl.index)}
    assign = lambda perm: [idx[q][perm[pos[q]]] for q in PLATES if perm[pos[q]] in idx[q]]
    obs_P = assign(np.arange(len(pl))); obs = both(obs_P)
    g = np.random.default_rng(int(os.environ.get("POSS_PERM_SEED", "20261024")))
    names = ["Z", "OE", "E1", "E2"]; short = {"Z": 1, "OE": -1, "E1": -1, "E2": -1}
    rows, info = [], []
    for _ in range(B):
        P = assign(g.permutation(len(pl)))
        st = both(P) if P else {}
        rows.append([st.get(k, np.nan) for k in names]); info.append(var_z(P) if P else np.nan)
    X = np.array(rows); info = np.array(info); ok = np.isfinite(X[:, 0])
    near = ok & (np.abs(info - var_z(obs_P)) <= 0.2 * var_z(obs_P))
    pv = lambda m, i: round(float((np.sum(short[names[i]] * X[m, i] >= short[names[i]] * obs[names[i]]) + 1) / (m.sum() + 1)), 4)
    q5 = lambda v: {k: round(float(x), 4) for k, x in zip(("q05", "q16", "median", "q84", "q95"), np.nanquantile(v, [0.05, 0.16, 0.5, 0.84, 0.95]))}
    out = {"permutations": int(ok.sum()), "near": int(near.sum()), "plates_observed": sum(not p["aside"] for p in obs_P),
           "observed": {k: round(float(obs[k]), 4) for k in names},
           "null": {k: q5(X[ok, i]) for i, k in enumerate(names)}, "null_near": {k: q5(X[near, i]) for i, k in enumerate(names)},
           "p": {k: pv(ok, i) for i, k in enumerate(names)}, "p_near": {k: pv(near, i) for i, k in enumerate(names)}}
    np.savez(D / "permutation_draws.npz", **{k: X[:, i] for i, k in enumerate(names)}, near=near)
    (D / "permutation_null.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
