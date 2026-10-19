"""Do transients cluster around US nuclear tests? VASCO's test as they ran it (Bruehl & Villarroel 2025; Bruehl et al.
2026), reproduced on their construction and on the full catalogue. The measurement itself is nuclear_fit.py.

Unit: an observing night; its value is the number of transients on its plates (each transient counted once; the public
machine-learning scores are not usable). Nights near tests: within one day of a US test. VASCO's test: the ratio of
transients per night near tests to other nights, against the same number of nights picked at random. Also how many
photographs the nights near tests hold, against other nights.

VASCO's construction: archive dates, their list of nights. Ours: every catalogue plate on its corrected date, nights
from Palomar's log (sample.py).

Then two corrections with physical reasons, one at a time: transients per photograph, and against what each
photograph's pointing predicts; each against VASCO's null and the test calendar slid by years.

python analysis/nuclear.py -> data/nuclear.json"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

import sample

ROOT = Path(__file__).resolve().parents[1]
NRAND = 20000
DAY0 = pd.Timestamp("1949-01-01")
r3 = lambda x: round(float(x), 3)

nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
TESTS = np.unique(((pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date) - DAY0).dt.days).values)


def near(days, offsets=(-1, 0, 1)):
    return np.isin(days, (TESTS[:, None] + np.array(offsets)[None, :]).ravel())


def random_nights(values, lab, rng):
    """VASCO's test: transients per night on the nights near tests over the same for all other nights, against random
    sets of nights of the same size; nights with no plates count as zero on both sides."""
    stat = lambda m: values[m].mean() / values[~m].mean()
    obs, k = stat(lab), int(lab.sum())
    R = np.empty(NRAND)
    for i in range(NRAND):
        m = np.zeros(len(values), bool); m[rng.choice(len(values), k, replace=False)] = True; R[i] = stat(m)
    p = np.mean(R >= obs)
    return {"real": r3(obs), "nights": k, "p": r3(p), "sigma": round(float(norm.isf(p)), 1),
            "hist": np.histogram(R, np.arange(0, 3.01, 0.05))[0].tolist()}


out = {}
# VASCO's construction: archive dates, their list of nights
vcat, _ = sample.vasco()
vn = sample.vasco_nights()
vv = vcat.groupby("date").size().reindex(vn, fill_value=0).values.astype(float)
out["vasco"] = random_nights(vv, near(((vn - DAY0).days).values), np.random.default_rng(42))
out["vasco"].pop("hist")

# ours: corrected dates, Palomar's nights
cat, plates = sample.ours()
nights = sample.nights()
days = ((nights - DAY0).days).values
values = cat.groupby("date").size().reindex(nights, fill_value=0).values.astype(float)
w = near(days)
out["perm"] = {"raw": random_nights(values, w, np.random.default_rng(42))}

# how many photographs the nights near tests hold
npl = plates.groupby("date").size().reindex(nights).fillna(0).values
out["steps"] = {"nights": len(nights), "nights_with_plates": int((npl > 0).sum()), "max_plates": int(npl.max()),
                "plates_per_night": round(float(npl.mean()), 2), "nights_window": int(w.sum()),
                "plates_window": round(float(npl[w].mean()), 2), "plates_other": round(float(npl[~w].mean()), 2),
                "per_plate": r3((values[w].sum() / npl[w].sum()) / (values[~w].sum() / npl[~w].sum())),
                "empty_window": int((npl[w] == 0).sum()), "empty_other": int((npl[~w] == 0).sum()), "nights_other": int((~w).sum())}

# corrections with physical reasons, one at a time: transients per photograph, then against what each photograph's
# pointing predicts (VASCO keep only transients with no known star nearby, so crowded fields near the Milky Way hold
# few). Each against two nulls: VASCO's random nights from the whole survey, and the real test calendar slid by 1 to 5
# years and up to 20 days, which keeps each test's season and the bunching of tests into series.
import patsy                                                                  # noqa: E402
import statsmodels.api as sm                                                  # noqa: E402
import plate_properties as PP                                                 # noqa: E402
from nuclear_model import P, SHIFTS, day, y_raw                               # noqa: E402

X = np.asarray(patsy.dmatrix("bs(b, df=4) + lebv + ecl_lat", PP.standardised(P))); a = 0.5
for _ in range(4):
    E = np.asarray(sm.GLM(y_raw, X, family=sm.families.NegativeBinomial(alpha=a)).fit().fittedvalues)
    a = float(np.sum(((y_raw - E) ** 2 - E) / E ** 2) / (len(y_raw) - X.shape[1]))
m0 = near(day)
per_night = lambda lab: values[lab].mean() / values[~lab].mean()
ratio = lambda m, e: (y_raw[m].sum() / e[m].sum()) / (y_raw[~m].sum() / e[~m].sum())
slid_nights = [near(days, tuple(np.array([-1, 0, 1]) + s)) for s in SHIFTS]
slid_plates = [near(day, tuple(np.array([-1, 0, 1]) + s)) for s in SHIFTS]
keep = [i for i, m in enumerate(slid_plates) if m.sum() >= 5 and (m & m0).sum() < 0.2 * m.sum()]
rng = np.random.default_rng(43)
random_nights = [rng.choice(len(nights), int(w.sum()), replace=False) for _ in range(NRAND)]
on_nights = lambda idx: np.isin(day, days[idx])
summ = lambda x, r: {"p": r3((np.sum(x >= r) + 1) / (len(x) + 1)), "median": r3(np.median(x)), "q16": r3(np.quantile(x, 0.16)),
                     "q84": r3(np.quantile(x, 0.84)), "q025": r3(np.quantile(x, 0.025)), "q975": r3(np.quantile(x, 0.975))}
out["corrections"] = []
for lab, f_real, f_vasco, f_cal in (
        ("Per night, as VASCO", lambda: per_night(w), lambda idx: per_night(np.isin(np.arange(len(nights)), idx)), lambda i: per_night(slid_nights[i])),
        ("Per photograph", lambda: ratio(m0, np.ones(len(y_raw))), lambda idx: ratio(on_nights(idx), np.ones(len(y_raw))), lambda i: ratio(slid_plates[i], np.ones(len(y_raw)))),
        ("Per photograph, allowing for pointing", lambda: ratio(m0, E), lambda idx: ratio(on_nights(idx), E), lambda i: ratio(slid_plates[i], E))):
    r = f_real(); V = np.array([f_vasco(idx) for idx in random_nights]); C = np.array([f_cal(i) for i in keep])
    out["corrections"].append({"measure": lab, "r": r3(r), "vasco_null": summ(V[np.isfinite(V)], r), "calendar_null": summ(C[np.isfinite(C)], r)})
b_ = np.abs(P.b.values)
out["pointing"] = {"median_low": float(np.median(y_raw[b_ < 10])), "median_high": float(np.median(y_raw[b_ > 40])),
                   "expected_near_over_other": r3(E[m0].mean() / E[~m0].mean()), "slides": len(keep),
                   "near_spring": r3(np.isin(P.date[m0].dt.month, [2, 3, 4, 5, 6]).mean())}

(ROOT / "data/nuclear.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "hist"} if k == "perm" else v) for k, v in out.items()}, default=str, indent=1))
