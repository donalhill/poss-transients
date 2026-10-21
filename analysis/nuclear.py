"""Do transients cluster around US nuclear tests? Every number in the page's nuclear section.

1. VASCO's test as they ran it (Bruehl & Villarroel 2025; Bruehl et al. 2026), reproduced on their construction and on
   the full catalogue.

Unit: an observing night; its value is the number of transients on its plates (each transient counted once; the public
machine-learning scores are not usable). Nights near tests: within one day of a US test. VASCO's test: the ratio of
transients per night near tests to other nights, against the same number of nights picked at random. Also how many
photographs the nights near tests hold, against other nights.

VASCO's construction: archive dates, their list of nights. Ours: every catalogue plate on its corrected date, nights
from Palomar's log (sample.py).

2. Two corrections with physical reasons, one at a time: transients per photograph, and against what each photograph's
   pointing predicts; each against VASCO's null and against the real test calendar shifted by 1 to 5 years and up to two
   months, leaving out plates taken within two weeks of a real test. The pointing curve refitted without the plates
   near tests, as a check that they do not shape their own expected counts. Checks on the null: how close any plate a
   shifted calendar counts comes to a real test, and an excess put into the plates near tests (scaled by 1.5 and 2),
   measured against the shifted calendars.
3. Everything else measurable about a plate, other definitions of a transient, and each test series left out.
4. Every window that includes the test day, up to a week either side; the rich plates that move the result between
   windows, and the windows without them.

python analysis/nuclear.py -> data/nuclear.json"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.spatial import cKDTree
from scipy.stats import norm

import plate_properties as PP
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
P = PP.properties(cat, plates)                                                # every plate and its properties
day = ((P.date - DAY0).dt.days).values
y_raw = P.n.values.astype(float)

SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-60, 61)]

E, X, a = PP.pointing(P)                                                      # the pointing curve (plate_properties.py)
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


def curve(y, use):
    """The pointing curve fitted to the plates in use, predicting every plate."""
    a_ = 0.5
    for _ in range(4):
        fit = sm.GLM(y[use], X[use], family=sm.families.NegativeBinomial(alpha=a_)).fit()
        E_ = fit.predict(X[use]); a_ = float(np.sum(((y[use] - E_) ** 2 - E_) / E_ ** 2) / (use.sum() - X.shape[1]))
    return fit.predict(X)


# the same curve fitted without the plates near tests, so they cannot shape their own expected counts
out["corrected"]["held_out"] = r3(ratio(m0, curve(y_raw, ~m0)))
# the null holds no real effect: every plate a shifted calendar counts is more than two weeks from a real test
counted = np.any([slid_plates[i] for i in keep], 0)
gap = int(np.abs(day[counted][:, None] - TESTS[None, :]).min())
assert gap > CLEAR
# and it does not absorb one put in: transients on the plates near tests scaled up, the curve refitted, measured as above
inj = []
for f in (1.5, 2.0):
    y_ = np.where(m0, np.round(y_raw * f), y_raw); E_ = curve(y_, np.ones(len(y_), bool))
    rat = lambda m: (y_[m].sum() / E_[m].sum()) / (y_[~m].sum() / E_[~m].sum())
    r_ = rat(m0); C_ = np.array([rat(slid_plates[i]) for i in keep])
    inj.append({"factor": f, "measured": r3(r_), "measured_over_real": r3(r_ / r_c), "null_median": r3(np.median(C_)),
                "p": r3((np.sum(C_ >= r_) + 1) / (len(C_) + 1))})
out["null_checks"] = {"closest_test_days": gap, "injected": inj}


# everything else checked, each on top of the corrected measure: every other plate property, other definitions of a
# transient, and each test series left out; each against its own slid calendars
def expected(y, terms=()):
    X_ = np.asarray(patsy.dmatrix(" + ".join([PP.POINTING, *terms]), PP.standardised(P))); a_ = 0.5
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


pcat = cat
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
shots["day"] = np.floor((pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601") - DAY0).dt.total_seconds() / 86400).astype(int)
shots["series"] = shots.source_url.str.extract(r"Operation_([A-Za-z]+)")[0].fillna(shots.Shot_Name)
y_cells = cat.assign(i=np.floor(cat.xi * 30).astype(int), j=np.floor(cat.eta * 30).astype(int)).groupby("plate") \
             .apply(lambda g: len(set(zip(g.i, g.j)))).reindex(P.index).fillna(0).values.astype(float)
y_inner = cat[(cat.xi.abs() < 2.5) & (cat.eta.abs() < 2.5)].groupby("plate").size().reindex(P.index).fillna(0).values.astype(float)
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

# every window that includes the test day, set by how many days before and after it reaches (up to a week), scored as
# above against the same shifted calendars: how often the strongest of them does at least as well as VASCO's by chance
R = 7
OFF = np.arange(-R, R + 1)
WIN = [(-b_, a_) for b_ in range(R + 1) for a_ in range(R + 1)]                # (first day, last day) from the test
ia, ib = np.array([a_ + R for a_, _ in WIN]), np.array([b_ + R + 1 for _, b_ in WIN])
Y, Et = y_raw.sum(), E.sum()


def scan(shift):
    """The ratio for every window, for the test calendar moved by shift days (plates near real tests left out)."""
    H = np.isin(day[:, None] - OFF[None, :], TESTS + shift) & (clear_plates[:, None] if shift else True)
    C = np.c_[np.zeros(len(day), int), np.cumsum(H, 1)]
    M = (C[:, ib] - C[:, ia] > 0).T.astype(float)                                   # windows x plates
    n1, e1, k = M @ y_raw, M @ E, M.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):                         # fewer than 5 plates: undetermined
        return np.where(k >= 5, (n1 / e1) / ((Y - n1) / (Et - e1)), np.nan)


real_w = scan(0)
null_w = np.array([scan(SHIFTS[i]) for i in keep])
iV = WIN.index((-1, 1))
p_vasco = (np.sum(null_w[:, iV] >= real_w[iV]) + 1) / (np.isfinite(null_w[:, iV]).sum() + 1)
# for each shifted calendar, its strongest window's local chance, by ranking every window among the shifted calendars;
# an undetermined window never counts as striking
nw = np.where(np.isfinite(null_w), null_w, -np.inf)
rank = np.sum(nw[None, :, :] >= nw[:, None, :], 1) / (len(keep) + 1)                  # shifted calendars x windows
fam = {"windows": len(WIN), "chance_vs_vasco": r3(np.mean(rank.min(1) <= p_vasco))}

# the rich plates that split the windows into two groups: the plate within a week of a test that most exceeds what its
# pointing predicts, and the two that do so most within a day of a test
lvl = (Y - y_raw[m0].sum()) / (Et - E[m0].sum())                             # the level of plates away from tests


def rich(offsets_, k):
    near_ = near(day, tuple(offsets_)); ex = np.where(near_, y_raw - lvl * E, -np.inf)
    out_ = []
    for j in np.argsort(ex)[::-1][:k]:
        off = day[j] - TESTS; off = off[np.abs(off) <= R]
        out_.append({"plate": P.index[j], "date": str(P.date.iloc[j].date()), "n": int(y_raw[j]), "expected": int(round(lvl * E[j])),
                     "days_from_test": int(off[np.argmin(np.abs(off))])})
    return out_


fam["rich_week"], fam["rich_day"] = rich(range(-R, R + 1), 1)[0], rich((-1, 0, 1), 2)
two = P.index.isin([x["plate"] for x in fam["rich_day"]]); test_day = near(day, (0,))
fam["rich_day_share"] = {"plates_test_day": int(test_day.sum()), "test_day": r3(y_raw[two & test_day].sum() / y_raw[test_day].sum())}
# every window again without those three plates (left out of both sides): a description, since they were picked by looking
drop3 = P.index.isin([fam["rich_week"]["plate"]] + [x["plate"] for x in fam["rich_day"]])


def ratio_without(a_, b_):
    m_ = near(day, tuple(range(a_, b_ + 1))); w_, o_ = m_ & ~drop3, ~m_ & ~drop3
    return (y_raw[w_].sum() / E[w_].sum()) / (y_raw[o_].sum() / E[o_].sum())


fam["cells"] = [{"before": -a_, "after": b_, "r": r3(real_w[j]), "r_without": r3(ratio_without(a_, b_))} for j, (a_, b_) in enumerate(WIN)]
fam["without"] = {"plates": list(P.index[drop3]), "median": r3(np.median([c_["r_without"] for c_ in fam["cells"]])),
                  "vasco": r3(ratio_without(-1, 1))}
out["windows"] = fam

b_ = np.abs(P.b.values)
out["pointing"] = {"median_low": float(np.median(y_raw[b_ < 10])), "median_high": float(np.median(y_raw[b_ > 40])),
                   "expected_near_over_other": r3(E[m0].mean() / E[~m0].mean()), "slides": len(keep),
                   "near_spring": r3(np.isin(P.date[m0].dt.month, [2, 3, 4, 5, 6]).mean())}

(ROOT / "data/nuclear.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "hist"} if k == "perm" else v) for k, v in out.items()}, default=str, indent=1))
