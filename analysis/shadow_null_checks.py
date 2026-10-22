"""Checks on the pooled shadow test's null (shadow_pooled.py): random placements of the real shadow outlines on plates the
shadow never crossed, each count scaled to the shadow plate's pointing (plate_properties.pointing, fitted without the
plates either shadow crosses). Every transient counts by its rating, as on the page (shadow_pooled.pooled).

1. Calibration: ten plates the shadow never crossed act as pretend shadow plates, with the real outlines on them, and are
   tested against the other uncrossed plates in the same way, the pointing curve refitted without them. With no signal,
   the share of tests below a level should be that level.
2. Leak: no host plate is crossed by either shadow.
3. Injection: a share of the summed ratings inside the real shadow removed, and the test run as on the page.
4. The pointing curve built term by term: how much each part moves the test.

python analysis/shadow_null_checks.py -> data/shadow_null_checks.json (not part of the default build: make shadow-checks)"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm

import plate_properties as PP
import sample
from shadow_pooled import counts, pooled

ROOT = Path(__file__).resolve().parents[1]
NCAL, NTRIAL, LEVELS = 500, 4000, (0.05, 0.11, 0.25, 0.5)

cat, plates = sample.ours()
(res, t, c, masks), topo = pooled(cat, plates, "geo"), pooled(cat, plates, "topo")
crossed = set(masks) | set(topo[3])
qi = pd.Series(np.arange(len(plates)), plates.index)
H = qi[[q for q in plates.index if q not in crossed]].values                   # host plates, as positions
K = list(masks); C = counts(c, plates, masks)                                   # outline x orientation x plate
P = PP.properties(cat, plates); S = PP.standardised(P); use = ~P.index.isin(crossed)
E = pd.Series(PP.pointing(P, use=use, y="nw")[0], P.index).reindex(plates.index).values
rng = np.random.default_rng(3)


def null(Ek, pool, n, E_=E):
    """n trials: each outline on a random plate from pool, turned and flipped at random, scaled to Ek over that plate's E."""
    q = pool[rng.integers(0, len(pool), (n, len(K)))]; d = rng.integers(0, 8, (n, len(K)))
    return (C[np.arange(len(K))[None, :], d, q] * Ek[None, :] / E_[q]).sum(1)


out = {"observed": res["observed"]}
Ek = E[qi[K].values]; s = null(Ek, H, 20000)
out["page"] = {"null_mean": round(float(s.mean()), 1), "p": round(float(np.mean(s <= res["observed"])), 3)}
# 1. calibration
ps = []
for _ in range(NCAL):
    fake = rng.choice(H, len(K), replace=False)
    o = C[np.arange(len(K)), 0, fake].sum()
    # the pointing curve refitted without the pretend shadow plates, as the page's is without the real ones
    E_f = pd.Series(PP.pointing(P, use=use & ~P.index.isin(plates.index[fake]), y="nw")[0], P.index).reindex(plates.index).values
    ps.append(np.mean(null(E_f[fake], np.setdiff1d(H, fake), NTRIAL, E_f) <= o))
out["calibration"] = {"runs": NCAL, "share_below": {str(v): round(float(np.mean(np.array(ps) <= v)), 3) for v in LEVELS}}
# 2. leak
out["hosts_crossed"] = len(set(plates.index[H]) & crossed)
# 3. injection
out["injected"] = {}
for f in (0.2, 0.4, 0.6):
    o = round(res["observed"] * (1 - f), 1)
    out["injected"][str(f)] = {"observed": o, "p": round(float(np.mean(s <= o)), 3)}
# 4. the pointing curve built term by term
y = P.nw.values.astype(float); out["terms"] = {}
for name, f in [("milky_way", "bs(b, df=4)"), ("milky_way_dust", "bs(b, df=4) + lebv"),
                ("milky_way_dust_ecliptic", "bs(b, df=4) + lebv + ecl_lat"), ("ecliptic", "ecl_lat"),
                ("centre_dust_with_equatorial_row", "bs(b, df=4) + lebv + ecl_lat + eq_row"), ("page_curve", PP.POINTING)]:
    X = np.asarray(patsy.dmatrix(f, S)); a = 0.5
    for _ in range(4):
        e_ = sm.GLM(y[use], X[use], family=sm.families.NegativeBinomial(alpha=a)).fit().predict(X)
        a = float(np.sum(((y[use] - e_[use]) ** 2 - e_[use]) / e_[use] ** 2) / (use.sum() - X.shape[1]))
    E_ = pd.Series(e_, P.index).reindex(plates.index).values
    q = H[rng.integers(0, len(H), (20000, len(K)))]; d = rng.integers(0, 8, (20000, len(K)))
    s_ = (C[np.arange(len(K))[None, :], d, q] * E_[qi[K].values][None, :] / E_[q]).sum(1)
    out["terms"][name] = {"expected": round(float(sum(t.f_area[k] * E_[qi[k]] for k in K)), 1), "p": round(float(np.mean(s_ <= res["observed"])), 3)}

(ROOT / "data/shadow_null_checks.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
