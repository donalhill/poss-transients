"""Do transients cluster around US nuclear tests? VASCO's test as they ran it (Bruehl & Villarroel 2025; Bruehl et al.
2026), reproduced on their construction and on the full catalogue. The measurement itself is nuclear_fit.py.

Unit: an observing night; its value is the number of transients on its plates (each transient counted once; the public
machine-learning scores are not usable). Nights near tests: within one day of a US test. VASCO's test: the ratio of
transients per night near tests to other nights, against the same number of nights picked at random. Also how many
photographs the nights near tests hold, against other nights.

VASCO's construction: archive dates, their list of nights. Ours: every catalogue plate on its corrected date, nights
from Palomar's log (sample.py).

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

(ROOT / "data/nuclear.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "hist"} if k == "perm" else v) for k, v in out.items()}, default=str, indent=1))
