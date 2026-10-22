"""Every set of ratings from ml_retrain.py, and the released ratings, put through the page's tests (analysis/ML_PLAN.md,
sections 2a and 5): VASCO's 2026 nuclear test on their construction (per night, their dates and list, against random
nights); ours per photograph and against a pointing curve fitted to the summed ratings (against the shifted calendars);
the top tenth of each set's ratings (edges from VASCO's construction), each transient counted once, against its own
pointing curve; and the shadow's plate-by-plate ratio with every transient counted by its rating.

python analysis/ml_retrain_eval.py -> data/ml_retrain.json"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import plate_properties as PP
import sample
import shadow_fit as SF
import shadow as SH
from shadow import HALF

ROOT = Path(__file__).resolve().parents[1]
RAT = ROOT / "data/external/vasco_ml/ratings"
rng = np.random.default_rng(7)
r3 = lambda x: round(float(x), 3)

cat, plates = sample.ours()
P = PP.properties(cat, plates)
vcat, _ = sample.vasco(); vn = sample.vasco_nights()

# the calendars, as nuclear.py builds them
DAY0 = pd.Timestamp("1949-01-01")
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
shots["day"] = np.floor((pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601") - DAY0).dt.total_seconds() / 86400).astype(int)
TESTS = np.unique(shots.day.values)
nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
TESTS_V = np.unique(((pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date) - DAY0).dt.days).values)
near = lambda days, off=(0, 1, 2), tests=TESTS: np.isin(days, (tests[:, None] + np.array(off)[None, :]).ravel())
day = ((P.date - DAY0).dt.days).values; m0 = near(day)
SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-60, 61)]
clear = ~near(day, tuple(range(-14, 15)))
slid = [near(day, tuple(np.array((0, 1, 2)) + s)) & clear for s in SHIFTS]
keep = [i for i, m in enumerate(slid) if abs(m.sum() - m0.sum()) <= 0.3 * m0.sum()]
vwin = near(((vn - DAY0).days).values, (-1, 0, 1), TESTS_V)
perm = [rng.choice(len(vn), int(vwin.sum()), replace=False) for _ in range(5000)]

# the shadow plate by plate, as shadow_fit.measure fits it (the fitted ratio only), the shadow seen from Palomar
c6 = (cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)
hosts = ~cat.plate.isin(SF.crossed_plates()) & c6
geo = {}
for q in sorted(SF.crossed_plates() & set(plates.index)):
    p = plates.loc[q]; sub = SH.sub_dark(p, SH.MAIN)
    if not sub.any() or (sub > 1 - 1e-9).all():
        continue
    c = cat[c6 & (cat.plate == q)]
    geo[q] = dict(sub=sub.mean(1), idx=c.index.values, s=SH.dark(c.xi.values, c.eta.values, p, SH.MAIN))
hidx = cat.index[hosts].values; hcode = pd.factorize(cat.plate[hosts].values)[0]
hcell = np.floor((cat.eta[hosts].values + HALF) / 0.1).astype(int) * 60 + np.floor((cat.xi[hosts].values + HALF) / 0.1).astype(int)


def shadow_r(ws):
    maps = np.zeros((hcode.max() + 1, 3600)); np.add.at(maps, (hcode, hcell), ws.reindex(hidx).values)
    acc = (maps / maps.sum(1, keepdims=True)).mean(0); acc /= acc.mean()
    L = 0.0
    for o in geo.values():
        u = ws[o["idx"]].values
        L = L + u.sum() / (u ** 2).sum() * SF.curve(o["s"], u, float(o["sub"] @ acc / acc.sum()))
    return float(SF.best(L))


def evaluate(prob):
    w_all = pd.Series(prob); ws = w_all.reindex(cat.index)
    nw = ws.groupby(cat.plate).sum().reindex(P.index).fillna(0).values
    vals = w_all.reindex(vcat.index).groupby(vcat.date).sum().reindex(vn, fill_value=0).values
    st = lambda m: vals[m].mean() / vals[~m].mean(); r_v = st(vwin)
    E = PP.pointing(P.assign(nw=nw), y="nw")[0]
    ratio = lambda y, e, m: (y[m].sum() / e[m].sum()) / (y[~m].sum() / e[~m].sum())
    cal = lambda y, e, r: (np.sum([ratio(y, e, slid[i]) >= r for i in keep]) + 1) / (len(keep) + 1)
    r_e = ratio(nw, E, m0)
    cut = np.quantile(w_all.reindex(vcat.index).values, 0.9)
    yt = (ws >= cut).groupby(cat.plate).sum().reindex(P.index).fillna(0).values.astype(float)
    Et = PP.pointing(P.assign(n=yt))[0]; r_t = ratio(yt, Et, m0)
    return {"vasco_ratio": r3(r_v), "vasco_p": round(float((np.sum([st(np.isin(np.arange(len(vn)), i)) >= r_v for i in perm]) + 1) / (len(perm) + 1)), 4),
            "per_photograph": r3(ratio(nw, np.ones(len(nw)), m0)), "pointing": r3(r_e), "pointing_p": r3(cal(nw, E, r_e)),
            "top_cut": r3(cut), "top_nuclear": r3(r_t), "top_nuclear_p": r3(cal(yt, Et, r_t)), "shadow_r": r3(shadow_r(ws))}


released = sample.ratings()
out = {"released": evaluate(released), "sets": {}}
runs = {}
for f in sorted(RAT.glob("*.npy")):
    labels, night, seed = f.stem.rsplit("_", 2)
    pr = np.load(f).astype(float); assert len(pr) == len(released)
    runs.setdefault(f"{labels}|{night}", []).append({"seed": int(seed), "rank_corr_released": r3(pd.Series(pr).rank().corr(pd.Series(released).rank())), **evaluate(pr)})
for k, rs in runs.items():
    labels, night = k.split("|")
    out["sets"][k] = {"labels": labels, "night_inputs": night, "runs": len(rs),
                      **{m: {"median": r3(np.median([r[m] for r in rs])), "min": r3(min(r[m] for r in rs)), "max": r3(max(r[m] for r in rs))}
                         for m in rs[0] if m != "seed"}}
(ROOT / "data/ml_retrain.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: {m: v[m]["median"] for m in ("vasco_p", "pointing", "pointing_p", "top_nuclear", "shadow_r", "rank_corr_released")} for k, v in out["sets"].items()}, indent=1))
