"""The nuclear-date test (POSS_SAMPLE: the sample).

Unit: a plate, counted over its central square; each transient counted by its weight (its VASCO-ML rating, as VASCO's
2026 analysis counts, or once). Each plate's expected count comes from the pointing model (expected.py): a negative
binomial with offset log(area) and terms in where the plate points; the per-plate scatter is in the dispersion, fitted
by maximum likelihood. Measure: observed over expected on the plates near a US test
(VASCO's window: the test day to two days after, UT) over the same on all other plates. Nulls: random
nights from Palomar's log, and the real test calendar slid by 1 to 5 years and up to two months (calendars holding
within 30% of the real one's plates), leaving out plates within two weeks of a real test; the curve refitted without the
plates near tests as a check; an excess put into the plates near tests (x1.5, x2), measured against the calendars; and
each plate left out in turn (the curve fixed).
python analysis/nuclear.py -> $POSS_SAMPLE/nuclear.json (POSS_NUC_WINDOW=-1,0,1: nuclear_window-1_0_1.json)"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import expected as EX                                                        # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
D = Path(os.environ["POSS_SAMPLE"])
NRAND, MATCH, CLEAR = 20000, 0.3, 14
DAY0 = pd.Timestamp("1949-01-01")
VWIN = tuple(int(v) for v in os.environ.get("POSS_NUC_WINDOW", "0,1,2").split(","))   # UT days from each test: VASCO's as applied,
SUFFIX = "" if VWIN == (0, 1, 2) else "_window" + "_".join(str(v) for v in VWIN)       # or -1,0,1, a day either side in UT
r3 = lambda x: round(float(x), 3)
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
TESTS = np.unique(np.floor((pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601") - DAY0).dt.total_seconds() / 86400).astype(int))
if os.environ.get("POSS_NUC_LIST") == "vasco":                               # VASCO's list of US tests, as their code has it
    nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")      # (it leaves out Operation Upshot-Knothole)
    TESTS = np.unique(((pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date) - DAY0).dt.days).values)
    SUFFIX += "_vascolist"


def near(days, offsets=VWIN):
    return np.isin(days, (TESTS[:, None] + np.array(offsets)[None, :]).ravel())


cat, pl = sample.load()
pl = pl.copy()
pl["n"] = cat.groupby("plate").w.sum().reindex(pl.index).fillna(0)
UNIT = bool((cat.w == 1).all())
X, off = EX.design(pl); y = pl.n.to_numpy(float); TERMS = EX.TERMS(pl)          # the pointing model (expected.py)
curve = lambda use, yy=y: EX.fit(yy, X, off, use)
E, alpha = curve(np.ones(len(y), bool))
day = ((pl.date - DAY0).dt.days).values
nights = sample.nights(); ndays = ((nights - DAY0).days).values
m0 = near(day)
ratio = lambda m, e, yy=y: (yy[m].sum() / e[m].sum()) / (yy[~m].sum() / e[~m].sum())
SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-60, 61)]
clear = ~near(day, tuple(range(-CLEAR, CLEAR + 1)))
slid = [near(day, tuple(np.array(VWIN) + s)) & clear for s in SHIFTS]
n0 = int(m0.sum()); keep = [i for i, m in enumerate(slid) if abs(m.sum() - n0) <= MATCH * n0]
assert int(np.abs(day[np.any([slid[i] for i in keep], 0)][:, None] - TESTS[None, :]).min()) > CLEAR
rng = np.random.default_rng(43)
nn = int(near(ndays).sum())
rand = [np.isin(day, ndays[rng.choice(len(ndays), nn, replace=False)]) for _ in range(NRAND)]
summ = lambda x, r: {"p": r3((np.sum(x >= r) + 1) / (len(x) + 1)), "median": r3(np.median(x)), "q16": r3(np.quantile(x, 0.16)),
                     "q84": r3(np.quantile(x, 0.84)), "q025": r3(np.quantile(x, 0.025)), "q975": r3(np.quantile(x, 0.975))}
r = ratio(m0, E)
C = np.array([ratio(slid[i], E) for i in keep]); V = np.array([ratio(m, E) for m in rand])
# without the pointing model: the raw excess per plate against the same shifted calendars (which keep the tests' bunching
# but not their season); and across the calendars, how closely the pointing model's expected excess per plate follows
# the raw one, how often a calendar's plates point at sky as rich as the real one's, and each calendar's share of spring
# plates (February to June) against the real one's
one = np.ones(len(y)); Craw = np.array([ratio(slid[i], one) for i in keep]); Cexp = np.array([ratio(slid[i], one, E) for i in keep])
spring = lambda m: float(np.mean(pl.date[m].dt.month.between(2, 6)))
out = {"plates": len(pl), "plates_near": n0, "transients": int(y.sum()) if UNIT else round(float(y.sum()), 1), "dispersion": r3(alpha), "ratio": r3(r),
       "calendar_null": summ(C[np.isfinite(C)], r), "random_nights_null": summ(V[np.isfinite(V)], r),
       "sigma_calendar": round(float(norm.isf((np.sum(C >= r) + 1) / (len(C) + 1))), 2), "calendars": len(keep),
       "held_out": r3(ratio(m0, curve(~m0)[0])), "per_photograph": r3(ratio(m0, np.ones(len(y)))), "terms": TERMS,
       "per_photograph_calendar_null": summ(Craw[np.isfinite(Craw)], ratio(m0, one)),
       "calendar_pointing": {"expected_real": r3(ratio(m0, one, E)), "expected_median": r3(np.median(Cexp)),
                             "corr_log_raw_expected": r3(np.corrcoef(np.log(Craw), np.log(Cexp))[0, 1]),
                             "share_expected_as_high": r3(np.mean(Cexp >= ratio(m0, one, E))),
                             "spring_real": r3(spring(m0)), "spring_median": r3(np.median([spring(slid[i]) for i in keep]))},
       "calendar_draws": np.round(C[np.isfinite(C)], 4).tolist(), "random_draws": np.round(V[np.isfinite(V)], 4).tolist()}   # for the page's figure
# why counting per night and ignoring the sky would show an excess:
# photographs per observing night near tests and on other nights, the share of nights with no plate, the share of plates
# near tests taken from February to June, how much more the sky terms alone expect per unit of exposure near tests, and
# the excess counted per night (all transients over all nights) beside the one per photograph
nm = near(ndays); has = np.isin(ndays, day); pn = np.array([np.sum(day == d_) for d_ in ndays])
rate = E / pl.expo.to_numpy()
out["context"] = {"nights_near": int(nm.sum()), "nights_other": int((~nm).sum()),
                  "plates_per_night_near": round(float(pn[nm].mean()), 2), "plates_per_night_other": round(float(pn[~nm].mean()), 2),
                  "empty_nights_near": r3(1 - has[nm].mean()), "empty_nights_other": r3(1 - has[~nm].mean()),
                  "spring_share_near": r3(pl.date[m0].dt.month.between(2, 6).mean()), "spring_share_other": r3(pl.date[~m0].dt.month.between(2, 6).mean()),
                  "sky_expects_near": r3(rate[m0].mean() / rate[~m0].mean()),
                  "per_night": r3((y[np.isin(day, ndays[nm])].sum() / nm.sum()) / (y[np.isin(day, ndays[~nm])].sum() / (~nm).sum()))}
# each plate left out in turn, the curve fixed: the largest change, and which plate
loo = np.array([ratio(np.delete(m0, j), np.delete(E, j), np.delete(y, j)) for j in range(len(y))])
j = int(np.argmax(np.abs(loo - r)))
out["leave_one_out"] = {"max_change": r3(loo[j] - r), "plate": pl.index[j], "near": bool(m0[j]), "range": [r3(loo.min()), r3(loo.max())]}
# an excess put into the plates near tests: their counts scaled, the curve refitted, measured against the calendars
for s in (1.5, 2.0):
    ys = y.copy(); ys[m0] = np.round(ys[m0] * s) if UNIT else ys[m0] * s
    Es, _ = curve(np.ones(len(y), bool), ys)
    Cs = np.array([ratio(slid[i], Es, ys) for i in keep]); rs = ratio(m0, Es, ys)
    out[f"injected_x{s}"] = {"ratio": r3(rs), "p_calendar": r3((np.sum(Cs >= rs) + 1) / (len(Cs) + 1))}
print(json.dumps(out, indent=1))
out["window_ut_days"] = list(VWIN)
(D / f"nuclear{SUFFIX}.json").write_text(json.dumps(out, indent=1))
