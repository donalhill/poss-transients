"""Do transients cluster around US nuclear tests? VASCO's test (Bruehl et al. 2026) on the full catalogue, and what it
measures once the number of photographs, where the telescope pointed and the other properties of each plate are allowed
for.

Unit: an observing night; its value is the number of transients on its plates (each transient counted once; the public
machine-learning scores are not usable). Nights near tests: within one day of a US test. VASCO's test: the ratio of
transients per night near tests to other nights, against the same number of nights picked at random.

Ours: every catalogue plate on its corrected date, nights from Palomar's log (sample.py). VASCO's construction (archive
dates, their list of nights) is run once, to reproduce their number.

python analysis/nuclear.py -> data/nuclear.json"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
from scipy.stats import norm

import plate_properties as PP
import sample

ROOT = Path(__file__).resolve().parents[1]
NRAND = 20000
DAY0 = pd.Timestamp("1949-01-01")
r3 = lambda x: round(float(x), 3)

nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
us = nt[nt.Country.str.upper() == "USA"]
TESTS = np.unique(((pd.to_datetime(us.Date) - DAY0).dt.days).values)


def near(days, offsets=(-1, 0, 1), shift=0):
    return np.isin(days, (TESTS[:, None] + shift + np.array(offsets)[None, :]).ravel())


def random_nights(values, e, lab, rng):
    """VASCO's test: transients on the nights near tests over their expected number, divided by the same for all other
    nights, against random sets of nights of the same size. With e = 1 for every night this is VASCO's ratio of
    transients per night; nights with no plates count as zero on both sides."""
    V, E, k = values.sum(), e.sum(), int(lab.sum())
    stat = lambda m: (values[m].sum() / e[m].sum()) / ((V - values[m].sum()) / (E - e[m].sum()))
    obs = stat(lab)
    R = np.empty(NRAND)
    for i in range(NRAND):
        m = np.zeros(len(values), bool); m[rng.choice(len(values), k, replace=False)] = True; R[i] = stat(m)
    p = np.mean(R >= obs)
    return {"real": r3(obs), "nights": k, "p": r3(p), "sigma": round(float(norm.isf(p)), 1),
            "hist": np.histogram(R, np.arange(0, 3.01, 0.05))[0].tolist()}


out = {}

# ---- VASCO's construction: archive dates, their list of nights
vcat, _ = sample.vasco()
vn = sample.vasco_nights()
vv = vcat.groupby("date").size().reindex(vn, fill_value=0).values.astype(float)
out["vasco"] = random_nights(vv, np.ones(len(vv)), near(((vn - DAY0).days).values), np.random.default_rng(42))
out["vasco"].pop("hist")

# ---- ours: corrected dates, Palomar's nights
cat, plates = sample.ours()
nights = sample.nights()
days = ((nights - DAY0).days).values
values = cat.groupby("date").size().reindex(nights, fill_value=0).values.astype(float)
lab = near(days)

# each plate's expected number of transients from every property of a plate (plate_properties.py), one negative
# binomial model fitted to all plates; a night's expected number is the sum over its plates
plates = PP.properties(cat, plates)
Z = PP.standardised(plates)
model = smf.negativebinomial(PP.FORMULA, Z).fit(disp=0, maxiter=500)
assert model.mle_retvals["converged"]
plates["e"] = model.predict(Z)
expected = plates.groupby("date").e.sum().reindex(nights).fillna(0).values
out["model"] = {k: {"coef": r3(model.params[k]), "z": round(float(model.tvalues[k]), 1)} for k in PP.TERMS}
# the distance to the Milky Way alone (one smooth curve rising to a plateau), for how much pointing on its own explains
sat = lambda x, A, b0, C: C + A * (1 - np.exp(-x / b0))
psat, _ = curve_fit(sat, plates.b, plates.n, p0=[200, 20, 30], maxfev=20000)
expected_b = plates.assign(e=sat(plates.b, *psat)).groupby("date").e.sum().reindex(nights).fillna(0).values
npl = plates.groupby("date").size().reindex(nights).fillna(0).values

rng = np.random.default_rng(42)
out["perm"] = {"raw": random_nights(values, np.ones(len(values)), lab, rng), "adj": random_nights(values, expected, lab, rng)}

# the nights near tests by single day: the test day, and the days before and after
out["lags"] = {str(k): r3(values[near(days, (k,))].mean() / values[~near(days, (k,))].mean()) for k in (-1, 0, 1)}

# VASCO's excess, step by step: per night = plates per night x transients per plate; per plate = expected per plate
# (where the plates pointed) x transients per plate over expected
w, o = lab, ~lab
ym = nights.to_period("M"); same = o & np.isin(ym, ym[w].unique())
out["steps"] = {"nights": len(nights), "nights_with_plates": int((npl > 0).sum()), "max_plates": int(npl.max()),
                "plates_per_night": round(float(npl.mean()), 2),
                "per_night": r3(values[w].mean() / values[o].mean()),
                "per_plate": r3((values[w].sum() / npl[w].sum()) / (values[o].sum() / npl[o].sum())),
                "per_expected": r3((values[w].sum() / expected[w].sum()) / (values[o].sum() / expected[o].sum())),
                "pointing": r3((expected_b[w].sum() / npl[w].sum()) / (expected_b[o].sum() / npl[o].sum())),
                "plates_window": round(float(npl[w].mean()), 2), "plates_other": round(float(npl[o].mean()), 2),
                "plates_same_months": round(float(npl[same].mean()), 2), "nights_window": int(w.sum()),
                "empty_window": int((npl[w] == 0).sum()), "empty_other": int((npl[o] == 0).sum()), "nights_other": int(o.sum())}

# what is left once photographs and pointing are allowed for: the plates near tests holding the most over expected, and
# the excess with each of them left out
pw = plates[plates.date.isin(nights[w])].assign(ratio=lambda d: d.n / d.e).sort_values("ratio", ascending=False)
base_n, base_e = values[o].sum(), expected[o].sum()
def excess_without(q):
    k = pw.drop(q); return (k.n.sum() / k.e.sum()) / (base_n / base_e)
out["top_plates"] = [{"plate": q, "n": int(r.n), "expected": round(float(r.e), 1), "ratio": round(float(r.ratio), 2),
                      "excess_without": r3(excess_without(q))} for q, r in pw.head(3).iterrows()]

# by month of the year: transients per plate, plates per night, distance from the Milky Way, nights near tests
mt = pd.DataFrame({"m": nights.month, "n": values, "p": npl, "w": w}).groupby("m").sum()
mb = plates.groupby(plates.date.dt.month).b.mean()
out["months"] = [{"m": "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()[i - 1], "n": round(float(mt.n[i] / (nights.month == i).sum()), 1), "pp": round(float(mt.n[i] / mt.p[i]), 1),
                  "ppn": round(float(mt.p[i] / (nights.month == i).sum()), 2), "b": round(float(mb[i]), 1), "win": int(mt.w[i])}
                 for i in range(1, 13)]
out["window_nights_by_month"] = pd.Series(nights[w].strftime("%Y-%m")).value_counts().sort_index().to_dict()

(ROOT / "data/nuclear.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "hist"} if k == "perm" else v) for k, v in out.items() if k != "months"}, default=str, indent=1))
