"""VASCO-ML retrained (vasco/ml_retrain.py, with VASCO's own training code): every saved set of ratings
(data/external/vasco_ml/ratings/LABELS_NIGHT_RUN.npy, in the catalogue's order; four sets of labels, with and without the
model's per-night inputs, several runs each) put in place of the released ratings, and the rated tests run again: the
nuclear test's ratio of observed to expected per plate near tests (VASCO's window, the pointing model refitted to the
set's counts) with its chance from the shifted calendars, as analysis/nuclear.py measures them, and the shadow's
ratio in shadow, observed over expected (analysis/shadow_fit.py). Each set is summarised by its runs' median, lowest, and highest.
POSS_SAMPLE=DIR (the rated sample) python diagnostics/vasco_ratings.py -> DIR/vasco_ratings.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "analysis"), str(ROOT / "vasco")]
import expected as EX                                                        # noqa: E402
import sample                                                                # noqa: E402
import shadow_fit as SF                                                      # noqa: E402
import vasco_sample                                                          # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
RAT = ROOT / "data/external/vasco_ml/ratings"
r3 = lambda x: round(float(x), 3)
cat, pl = sample.load()
pl = pl.copy()
X, off = EX.design(pl)
# the calendars, as nuclear.py builds them (VASCO's window as applied: the test day to two days after, UT)
DAY0 = pd.Timestamp("1949-01-01")
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
TESTS = np.unique(np.floor((pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601") - DAY0).dt.total_seconds() / 86400).astype(int))
near = lambda days, off_: np.isin(days, (TESTS[:, None] + np.array(off_)[None, :]).ravel())
day = ((pl.date - DAY0).dt.days).values; m0 = near(day, (0, 1, 2))
SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-60, 61)]
clear = ~near(day, tuple(range(-14, 15)))
slid = [near(day, tuple(np.array((0, 1, 2)) + s)) & clear for s in SHIFTS]
keep = [i for i, m in enumerate(slid) if abs(m.sum() - m0.sum()) <= 0.3 * m0.sum()]
ratio = lambda y, e, m: (y[m].sum() / e[m].sum()) / (y[~m].sum() / e[~m].sum())


def evaluate(prob):
    c = cat.assign(w=prob[cat.row.values])
    y = c.groupby("plate").w.sum().reindex(pl.index).fillna(0).to_numpy(float)
    E, _ = EX.fit(y, X, off, np.ones(len(y), bool)); r = ratio(y, E, m0)
    C = np.array([ratio(y, E, slid[i]) for i in keep])
    return {"nuclear": r, "nuclear_p": (np.sum(C >= r) + 1) / (len(C) + 1), "shadow_r": SF.measure(c, pl, extra=False)["fit"]["r"]}


released = vasco_sample.ratings()
runs = {}
for f in sorted(RAT.glob("*.npy")):
    labels, night, _ = f.stem.rsplit("_", 2)
    pr = np.load(f).astype(float); assert len(pr) == len(released)
    runs.setdefault((labels, night), []).append(evaluate(pr))
    print(f.stem, {k: r3(v) for k, v in runs[(labels, night)][-1].items()}, flush=True)
out = {"released": {k: r3(v) for k, v in evaluate(released).items()}, "sets": {}}
for (labels, night), rs in runs.items():
    out["sets"][f"{labels}|{night}"] = {"labels": labels, "night_inputs": night, "runs": len(rs),
                                        **{m: {"median": r3(np.median([r[m] for r in rs])), "min": r3(min(r[m] for r in rs)), "max": r3(max(r[m] for r in rs))}
                                           for m in rs[0]}}
(D / "vasco_ratings.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: {m: v[m]["median"] for m in ("nuclear", "nuclear_p", "shadow_r")} for k, v in out["sets"].items()}, indent=1))
