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
# years and up to two months, which keeps each test's time of year and the bunching of tests into series; slid
# calendars leave out plates taken within two weeks of a real test.
import patsy                                                                  # noqa: E402
import statsmodels.api as sm                                                  # noqa: E402
import plate_properties as PP                                                 # noqa: E402
from nuclear_model import P, day, y_raw                                       # noqa: E402

SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-60, 61)]

X = np.asarray(patsy.dmatrix("bs(b, df=4) + lebv + ecl_lat", PP.standardised(P))); a = 0.5
for _ in range(4):
    E = np.asarray(sm.GLM(y_raw, X, family=sm.families.NegativeBinomial(alpha=a)).fit().fittedvalues)
    a = float(np.sum(((y_raw - E) ** 2 - E) / E ** 2) / (len(y_raw) - X.shape[1]))
m0 = near(day)
per_night = lambda lab: values[lab].mean() / values[~lab].mean()
ratio = lambda m, e: (y_raw[m].sum() / e[m].sum()) / (y_raw[~m].sum() / e[~m].sum())
# a slid calendar never counts a plate (or night) taken within two weeks of a real test, so the null holds none of them
CLEAR = 14
clear_plates, clear_nights = ~near(day, tuple(range(-CLEAR, CLEAR + 1))), ~near(days, tuple(range(-CLEAR, CLEAR + 1)))
slid_nights = [near(days, tuple(np.array([-1, 0, 1]) + s)) & clear_nights for s in SHIFTS]
slid_plates = [near(day, tuple(np.array([-1, 0, 1]) + s)) & clear_plates for s in SHIFTS]
keep = [i for i, m in enumerate(slid_plates) if m.sum() >= 5]
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
# the corrected measure against the slid calendars, in full, for the chart
r_c = ratio(m0, E); C = np.array([ratio(slid_plates[i], E) for i in keep])
out["corrected"] = {"real": r3(r_c), "p": r3((np.sum(C >= r_c) + 1) / (len(C) + 1)), "sigma": round(float(norm.isf((np.sum(C >= r_c) + 1) / (len(C) + 1))), 1),
                    "median": r3(np.median(C)), "hist": np.histogram(C, np.arange(0, 3.01, 0.05))[0].tolist()}


# everything else checked, each on top of the corrected measure: every other plate property, other definitions of a
# transient, and each test series left out; each against its own slid calendars
def expected(y, terms=()):
    X_ = np.asarray(patsy.dmatrix(" + ".join(["bs(b, df=4)", "lebv", "ecl_lat", *terms]), PP.standardised(P))); a_ = 0.5
    for _ in range(4):
        mu = np.asarray(sm.GLM(y, X_, family=sm.families.NegativeBinomial(alpha=a_)).fit().fittedvalues)
        a_ = float(np.sum(((y - mu) ** 2 - mu) / mu ** 2) / (len(y) - X_.shape[1]))
    return mu


def check(group, label, y=y_raw, e=None, tests=TESTS):
    e = expected(y) if e is None else e
    w_ = lambda s: np.isin(day, (tests[:, None] + s + np.array([-1, 0, 1])[None, :]).ravel()) & (clear_plates if s else True)
    gone = m0 & ~w_(0)                                                   # plates near a left-out series: in neither group
    rt = lambda m: (y[m & ~gone].sum() / e[m & ~gone].sum()) / (y[~m & ~gone].sum() / e[~m & ~gone].sum())
    r = rt(w_(0)); Cs = np.array([rt(w_(SHIFTS[i])) for i in keep])
    return {"group": group, "label": label, "r": r3(r), "p": r3((np.sum(Cs >= r) + 1) / (len(Cs) + 1))}


from scipy.spatial import cKDTree                                             # noqa: E402
from nuclear_model import cat as pcat, shots, y_cells, y_inner               # noqa: E402
iso = np.zeros(len(pcat), bool)
for _, idx in pcat.groupby("plate").indices.items():
    xy = np.c_[pcat.xi.values[idx], pcat.eta.values[idx]] * 60
    iso[idx] = cKDTree(xy).query(xy, k=2)[0][:, 1] >= 5 if len(xy) > 1 else True
central = (np.maximum(pcat.xi.abs(), pcat.eta.abs()) < 2.5).values
count = lambda k: pcat[k].groupby("plate").size().reindex(P.index).fillna(0).values.astype(float)
props = {"depth": "How deep the plate goes", "fwhm": "How sharp its stars are", "exp": "Exposure time", "airmass": "Height above the horizon",
         "moon_up": "Moon above the horizon", "moon_sep": "Distance from the Moon", "quality_bad": "VASCO’s plate-quality flag"}
out["checks"] = [check("Plate properties", lab, e=expected(y_raw, (t,))) for t, lab in props.items()]
out["checks"].append(check("Plate properties", "All of them together", e=expected(y_raw, tuple(props))))
out["checks"] += [check("What counts as a transient", "2′ cells holding a transient", y=y_cells), check("What counts as a transient", "Central 5 × 5° of each plate", y=y_inner),
                  check("What counts as a transient", "No other transient within 5′", y=count(iso)), check("What counts as a transient", "Central, with no other within 5′", y=count(iso & central))]
for sname, g_ in shots.groupby("series"):
    t_ = np.setdiff1d(TESTS, g_.day.values)
    if (m0 & ~np.isin(day, (t_[:, None] + np.array([-1, 0, 1])[None, :]).ravel())).any():
        out["checks"].append(check("Test series left out", f"Without Operation {sname}", e=E, tests=t_))

# every window of consecutive days [a, b] within 14 days of a test, each scored as above against the same slid
# calendars: where VASCO's window ranks among them, and how often the best window does as well by chance
R = 14
OFF = np.arange(-R, R + 1)


def offsets(shift):
    """Plates x offsets: whether the plate was taken that many days from a (shifted) test; cumulated over offsets."""
    H = np.isin(day[:, None] - OFF[None, :], TESTS + shift) & (clear_plates[:, None] if shift else True)
    return np.c_[np.zeros(len(day), int), np.cumsum(H, 1)]


WIN = [(a_, b_) for a_ in OFF for b_ in OFF if b_ >= a_]
ia, ib = np.array([a_ + R for a_, _ in WIN]), np.array([b_ + R + 1 for _, b_ in WIN])
Y, Et = y_raw.sum(), E.sum()


def scan(shift):
    C = offsets(shift); M = (C[:, ib] - C[:, ia] > 0).T.astype(float)                  # windows x plates
    n1, e1, k = M @ y_raw, M @ E, M.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(k >= 5, (n1 / e1) / ((Y - n1) / (Et - e1)), np.nan)


real_w = scan(0)
null_w = np.array([scan(SHIFTS[i]) for i in keep])
use = np.isfinite(real_w)
loc = np.array([(np.sum(null_w[:, j] >= real_w[j]) + 1) / (np.isfinite(null_w[:, j]).sum() + 1) if use[j] else np.nan for j in range(len(WIN))])
# each slide's own local p for every window, by ranking it among the slides
nw = np.where(np.isfinite(null_w), null_w, -np.inf)
rank = (np.sum(nw[None, :, :] >= nw[:, None, :], 1)) / (len(keep) + 1)                  # slides x windows
iV = WIN.index((-1, 1))
span = np.array([max(abs(a_), abs(b_)) for a_, b_ in WIN])
out["scan"] = {"windows": int(use.sum()), "vasco_r": r3(real_w[iV]), "vasco_p_local": r3(loc[iV]),
               "vasco_rank": int(np.sum(loc[use] <= loc[iV])),
               "share_slides_rank_as_high": r3(np.mean([np.sum(rank[t, use] <= rank[t, iV]) <= np.sum(loc[use] <= loc[iV]) for t in range(len(keep))])),
               "grid": [{"a": int(a_), "b": int(b_), "r": r3(real_w[j]), "p": r3(loc[j]), "z": round(float(np.clip(norm.isf(loc[j]), -4, 4)), 2)}
                        for j, (a_, b_) in enumerate(WIN) if use[j]], "by_size": []}
for d_ in (3, 7, 14):
    c_ = np.flatnonzero(use & (span <= d_)); best = c_[np.nanargmin(loc[c_])]
    pg = (np.sum(rank[:, c_].min(1) <= loc[best]) + 1) / (len(keep) + 1)
    out["scan"]["by_size"].append({"days": d_, "windows": len(c_), "best": [int(WIN[best][0]), int(WIN[best][1])], "best_r": r3(real_w[best]),
                                   "best_p_local": r3(loc[best]), "p_global": r3(pg)})

# local to global: each slid calendar's best window over the whole scan, as a local significance
best_slide = norm.isf(rank[:, use].min(1)); jb = np.flatnonzero(use)[np.nanargmin(loc[use])]
pg = out["scan"]["by_size"][2]["p_global"]
out["global"] = {"best_local_sigma": r3(norm.isf(loc[jb])), "vasco_local_sigma": r3(norm.isf(loc[iV])), "p_global": pg, "global_sigma": r3(norm.isf(pg)),
                 "median_best_sigma": r3(np.median(best_slide)), "best_local_p": r3(loc[jb]), "vasco_local_p": r3(loc[iV]),
                 "curve": [{"sigma": r3(z_), "p_best": r3((np.sum(best_slide >= z_) + 1) / (len(best_slide) + 1))}
                           for z_ in np.arange(0, norm.isf(5 / (len(keep) + 1)), 0.05)]}            # while the slides resolve it
bw = out["scan"]["by_size"][2]["best"]; inb = np.isin(day, (TESTS[:, None] + np.arange(bw[0], bw[1] + 1)[None, :]).ravel())
top = np.argmax(np.where(inb, y_raw - E * (Y - y_raw[inb].sum()) / (Et - E[inb].sum()), -np.inf))
out["scan"]["best_plate"] = {"plate": P.index[top], "date": str(P.date.iloc[top].date()), "n": int(y_raw[top]), "expected": round(float(E[top])),
                             "share_of_window": r3(y_raw[top] / y_raw[inb].sum()), "share_of_catalogue": r3(y_raw[top] / Y)}

b_ = np.abs(P.b.values)
out["pointing"] = {"median_low": float(np.median(y_raw[b_ < 10])), "median_high": float(np.median(y_raw[b_ > 40])),
                   "expected_near_over_other": r3(E[m0].mean() / E[~m0].mean()), "slides": len(keep),
                   "near_spring": r3(np.isin(P.date[m0].dt.month, [2, 3, 4, 5, 6]).mean())}

(ROOT / "data/nuclear.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "hist"} if k == "perm" else v) for k, v in out.items()}, default=str, indent=1))
