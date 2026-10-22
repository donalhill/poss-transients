"""Do transients cluster around US nuclear tests? Every number in the page's nuclear section.

1. VASCO's test as they ran it (Bruehl & Villarroel 2025; Bruehl et al. 2026), reproduced on their construction and on
   the full catalogue.

Unit: an observing night; its value is the summed machine-learning ratings of the transients on its plates (VASCO-ML's
rating of how likely each is to be real, sample.py), as Bruehl et al. (2026) count theirs; how far the result depends on
the labels that trained the ratings is in ml_retrain_eval.py, and the same test with each transient counted once in
section 6. Nights near tests: VASCO's window, within one day of a US test by local evening
date (in UT, the test day to two days after; VWIN). VASCO's test: the ratio of
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
5. A 95% upper limit on an excess near tests, with the limit expected with none; and the ratings in bands (VASCO-ML's
   decile edges): the excess in each band, every band at once, its trend with rating.
6. The same tests with each transient counted once: VASCO's construction, and ours at each step.

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

# US tests by UT date: DOE/NV-209 (data/us_shot_times.csv). VASCO's list (their repository) lacks all eleven shots of
# Operation Upshot-Knothole (March-June 1953); it is used only to reproduce their result
nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
TESTS_V = np.unique(((pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date) - DAY0).dt.days).values)
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
shots["day"] = np.floor((pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601") - DAY0).dt.total_seconds() / 86400).astype(int)
shots["series"] = shots.source_url.str.extract(r"Operation_([A-Za-z]+)")[0].fillna(shots.Shot_Name)
TESTS = np.unique(shots.day.values)
# VASCO's window, a night within one day of a test, on their nights' local evening dates: a night's UT date is its
# evening date + 1, and the tests fell in the Nevada (or Pacific) morning, on the same UT date, so in UT the window runs
# from the test day to two days after
VWIN = (0, 1, 2)
MATCH = 0.3                          # shifted calendars kept: those holding within 30% of the real calendar's plates


def near(days, offsets=VWIN, tests=None):
    tests = TESTS if tests is None else tests
    return np.isin(days, (tests[:, None] + np.array(offsets)[None, :]).ravel())


def random_nights(values, lab, rng):
    """VASCO's test: transients per night on the nights near tests over the same for all other nights, against random
    sets of nights of the same size; nights with no plates count as zero on both sides."""
    stat = lambda m: values[m].mean() / values[~m].mean()
    obs, k = stat(lab), int(lab.sum())
    R = np.empty(NRAND)
    for i in range(NRAND):
        m = np.zeros(len(values), bool); m[rng.choice(len(values), k, replace=False)] = True; R[i] = stat(m)
    p = (np.sum(R >= obs) + 1) / (NRAND + 1)
    return {"real": r3(obs), "nights": k, "trials": NRAND, "p": r3(p), "sigma": round(float(norm.isf(p)), 1),
            "hist": np.histogram(R, np.arange(0, 3.01, 0.05))[0].tolist()}


out = {}
# VASCO's construction: archive dates, their list of nights
vcat, _ = sample.vasco()
vn = sample.vasco_nights()
vv = vcat.groupby("date").w.sum().reindex(vn, fill_value=0).values.astype(float)
out["vasco"] = random_nights(vv, near(((vn - DAY0).days).values, (-1, 0, 1), tests=TESTS_V), np.random.default_rng(42))   # local evening dates
out["vasco"].pop("hist")

# ours: corrected dates, Palomar's nights
cat, plates = sample.ours()
nights = sample.nights()
days = ((nights - DAY0).days).values
values = cat.groupby("date").w.sum().reindex(nights, fill_value=0).values.astype(float)
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
y_raw = P.nw.values.astype(float)                                            # each plate's summed ratings

SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-60, 61)]

E, X, a = PP.pointing(P, y="nw")                                              # the pointing curve (plate_properties.py), on the ratings
m0 = near(day)
per_night = lambda lab: values[lab].mean() / values[~lab].mean()
ratio = lambda m, e: (y_raw[m].sum() / e[m].sum()) / (y_raw[~m].sum() / e[~m].sum())
# a slid calendar never counts a plate (or night) taken within two weeks of a real test, so the null holds none of them
CLEAR = 14
clear_plates, clear_nights = ~near(day, tuple(range(-CLEAR, CLEAR + 1))), ~near(days, tuple(range(-CLEAR, CLEAR + 1)))
slid_nights = [near(days, tuple(np.array(VWIN) + s)) & clear_nights for s in SHIFTS]
slid_plates = [near(day, tuple(np.array(VWIN) + s)) & clear_plates for s in SHIFTS]
# calendars that hold about as many plates as the real one: a small calendar's ratio swings more, so unmatched sizes
# would widen the null
n0 = int(m0.sum())
keep = [i for i, m in enumerate(slid_plates) if abs(m.sum() - n0) <= MATCH * n0]
rng = np.random.default_rng(43)
rand_idx = [rng.choice(len(nights), int(w.sum()), replace=False) for _ in range(NRAND)]
on_nights = lambda idx: np.isin(day, days[idx])
summ = lambda x, r: {"p": r3((np.sum(x >= r) + 1) / (len(x) + 1)), "median": r3(np.median(x)), "q16": r3(np.quantile(x, 0.16)),
                     "q84": r3(np.quantile(x, 0.84)), "q025": r3(np.quantile(x, 0.025)), "q975": r3(np.quantile(x, 0.975))}
out["corrections"] = []
for lab, f_real, f_vasco, f_cal in (
        ("Per night, as VASCO", lambda: per_night(w), lambda idx: per_night(np.isin(np.arange(len(nights)), idx)), lambda i: per_night(slid_nights[i])),
        ("Per photograph", lambda: ratio(m0, np.ones(len(y_raw))), lambda idx: ratio(on_nights(idx), np.ones(len(y_raw))), lambda i: ratio(slid_plates[i], np.ones(len(y_raw)))),
        ("Per photograph, allowing for pointing", lambda: ratio(m0, E), lambda idx: ratio(on_nights(idx), E), lambda i: ratio(slid_plates[i], E))):
    r = f_real(); V = np.array([f_vasco(idx) for idx in rand_idx]); C = np.array([f_cal(i) for i in keep])
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
out["calendars"] = {"made": len(SHIFTS), "kept": len(keep), "plates_real": n0, "match": MATCH}
# the same measure with VASCO's test list, against its own shifted calendars (matched in the same way)
mV = near(day, tests=TESTS_V); nV = int(mV.sum())
slid_V = [near(day, tuple(np.array(VWIN) + s), tests=TESTS_V) & clear_plates for s in SHIFTS]
CV = np.array([ratio(m, E) for m in slid_V if abs(m.sum() - nV) <= MATCH * nV]); rV = ratio(mV, E)
out["vasco_list"] = {"plates": nV, "r": r3(rV), "p": r3((np.sum(CV >= rV) + 1) / (len(CV) + 1)), "calendars": len(CV),
                     "missing_shots": int(len(np.setdiff1d(TESTS, TESTS_V))), "extra_plates": n0 - nV}
# the null holds no real effect: every plate a shifted calendar counts is more than two weeks from a real test
counted = np.any([slid_plates[i] for i in keep], 0)
gap = int(np.abs(day[counted][:, None] - TESTS[None, :]).min())
assert gap > CLEAR
# and it does not absorb one put in: transients on the plates near tests scaled up, the curve refitted, measured as above
inj = []
for f in (1.5, 2.0):
    y_ = np.where(m0, y_raw * f, y_raw); E_ = curve(y_, np.ones(len(y_), bool))
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
    w_ = lambda s: np.isin(day, (tests[:, None] + s + np.array(VWIN)[None, :]).ravel()) & (clear_plates if s else True)
    gone = m0 & ~w_(0)                                                   # plates near a left-out series: in neither group
    rt = lambda m: (y[m & ~gone].sum() / e[m & ~gone].sum()) / (y[~m & ~gone].sum() / e[~m & ~gone].sum())
    r = rt(w_(0)); Cs = np.array([rt(w_(SHIFTS[i])) for i in keep])
    return {"group": group, "label": label, "r": r3(r), "p": r3((np.sum(Cs >= r) + 1) / (len(Cs) + 1))}


pcat = cat
y_cells = cat.assign(i=np.floor(cat.xi * 30).astype(int), j=np.floor(cat.eta * 30).astype(int)).groupby("plate") \
             .apply(lambda g: len(set(zip(g.i, g.j)))).reindex(P.index).fillna(0).values.astype(float)
y_inner = cat[(cat.xi.abs() < 2.5) & (cat.eta.abs() < 2.5)].groupby("plate").w.sum().reindex(P.index).fillna(0).values.astype(float)
iso = np.zeros(len(pcat), bool)
for _, idx in pcat.groupby("plate").indices.items():
    xy = np.c_[pcat.xi.values[idx], pcat.eta.values[idx]] * 60
    iso[idx] = cKDTree(xy).query(xy, k=2)[0][:, 1] >= 5 if len(xy) > 1 else True
central = (np.maximum(pcat.xi.abs(), pcat.eta.abs()) < 2.5).values
count = lambda k: pcat[k].groupby("plate").w.sum().reindex(P.index).fillna(0).values.astype(float)
props = {"depth": "How deep the plate goes", "fwhm": "How sharp its stars are", "exp": "Exposure time", "airmass": "Height above the horizon",
         "moon_up": "Moon above the horizon", "moon_sep": "Distance from the Moon", "quality_bad": "VASCO’s plate-quality flag"}
out["checks"] = [check("Plate properties", lab, e=expected(y_raw, (t,))) for t, lab in props.items()]
out["checks"].append(check("Plate properties", "All of them together", e=expected(y_raw, tuple(props))))
out["checks"] += [check("What counts as a transient", "2′ cells holding a transient", y=y_cells), check("What counts as a transient", "Central 5 × 5° of each plate", y=y_inner),
                  check("What counts as a transient", "No other transient within 5′", y=count(iso)), check("What counts as a transient", "Central, with no other within 5′", y=count(iso & central))]
for sname, g_ in shots.groupby("series"):
    t_ = np.setdiff1d(TESTS, g_.day.values)
    if (m0 & ~np.isin(day, (t_[:, None] + np.array(VWIN)[None, :]).ravel())).any():
        out["checks"].append(check("Test series left out", f"Without Operation {sname}", e=E, tests=t_))

# every window that includes the test day, set by how many days before and after it reaches (up to a week), scored as
# above against the same shifted calendars. Windows are counted as VASCO count theirs, by local evening date, so theirs
# is (-1, 1); in UT each runs a day later (VWIN)
R = 7
UT = 1                                                                        # UT offset = local evening offset + 1
OFF = np.arange(-R + UT, R + UT + 1)                                          # UT offsets the windows can reach
WIN = [(-b_, a_) for b_ in range(R + 1) for a_ in range(R + 1)]                # (first, last) local day from the test
ia, ib = np.array([a_ + UT - OFF[0] for a_, _ in WIN]), np.array([b_ + UT - OFF[0] + 1 for _, b_ in WIN])
assert WIN.index((-1, 1)) == [w_ for w_, (a_, b_) in enumerate(WIN) if (a_ + UT, b_ + UT) == (VWIN[0], VWIN[-1])][0]
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
# and the real calendar's strongest window, by the same local chance: how often chance gives a strongest window as striking
rw = np.where(np.isfinite(real_w), real_w, -np.inf)
p_real = (np.sum(nw >= rw[None, :], 0) + 1) / (len(keep) + 1)                          # each real window's local chance
fam = {"windows": len(WIN), "p_vasco": r3(p_vasco), "chance_vs_vasco": r3(np.mean(rank.min(1) <= p_vasco)),
       "best": {"window": list(WIN[int(np.argmin(p_real))]), "r": r3(real_w[int(np.argmin(p_real))]), "p_local": r3(p_real.min()),
                "p_global": r3(np.mean(rank.min(1) <= p_real.min()))}}

# the rich plates that split the windows into two groups: the plate within a week of a test that most exceeds what its
# pointing predicts, and the two that do so most in VASCO's window
lvl = (Y - y_raw[m0].sum()) / (Et - E[m0].sum())                             # the level of plates away from tests


def rich(offsets_, k):
    near_ = near(day, tuple(offsets_)); ex = np.where(near_, y_raw - lvl * E, -np.inf)
    out_ = []
    for j in np.argsort(ex)[::-1][:k]:
        off = day[j] - TESTS - UT; off = off[np.abs(off) <= R]                     # local days from the test
        out_.append({"plate": P.index[j], "date": str(P.date.iloc[j].date()), "n": round(float(y_raw[j]), 1), "expected": int(round(lvl * E[j])),
                     "days_from_test": int(off[np.argmin(np.abs(off))])})
    return out_


fam["rich_week"], fam["rich_day"] = rich(OFF, 1)[0], rich(VWIN, 2)
two = P.index.isin([x["plate"] for x in fam["rich_day"]])
fam["rich_day_share"] = {"plates_window": int(m0.sum()), "window": r3(y_raw[two & m0].sum() / y_raw[m0].sum())}
# every window again without those three plates (left out of both sides): a description, since they were picked by looking
drop3 = P.index.isin([fam["rich_week"]["plate"]] + [x["plate"] for x in fam["rich_day"]])


def ratio_without(a_, b_):
    m_ = near(day, tuple(range(a_ + UT, b_ + UT + 1))); w_, o_ = m_ & ~drop3, ~m_ & ~drop3
    return (y_raw[w_].sum() / E[w_].sum()) / (y_raw[o_].sum() / E[o_].sum())


fam["cells"] = [{"before": -a_, "after": b_, "r": r3(real_w[j]), "r_without": r3(ratio_without(a_, b_))} for j, (a_, b_) in enumerate(WIN)]
fam["without"] = {"plates": list(P.index[drop3]), "median": r3(np.median([c_["r_without"] for c_ in fam["cells"]])),
                  "vasco": r3(ratio_without(-1, 1))}
out["windows"] = fam

# the ratings in bands (VASCO-ML's decile edges), each transient counted once within its band, as their decile analysis
# counts them: in each band the excess near tests against that band's own pointing curve, scored against the shifted
# calendars; every band at once (Fisher's combination of the bands' chances, against the same for each calendar with
# every band shifted together, each calendar's chances ranked among the others); the trend of the excess with each
# band's mean rating; and the most extreme band, allowing for its being one of ten
band = sample.decile(cat.w)
YB = np.stack([cat[band == k].groupby("plate").size().reindex(P.index).fillna(0).values.astype(float) for k in range(10)])
EB = np.stack([PP.pointing(P.assign(n=YB[k]))[0] for k in range(10)])
rat = lambda y, e, m: (y[m].sum() / e[m].sum()) / (y[~m].sum() / e[~m].sum())
bands = lambda m: np.array([rat(YB[k], EB[k], m) for k in range(10)])
R0 = bands(m0); RC = np.array([bands(slid_plates[i]) for i in keep])                       # calendars x bands
p0 = (np.sum(RC >= R0, 0) + 1) / (len(RC) + 1)
pC = np.array([(np.sum(np.delete(RC, i, 0) >= RC[i], 0) + 1) / len(RC) for i in range(len(RC))])
F0, FC = -2 * np.log(p0).sum(), -2 * np.log(pC).sum(1)
mr = np.array([cat.w[band == k].mean() for k in range(10)])
slope = lambda r_: np.polyfit(mr, np.log(r_), 1)[0]
s0, sC = slope(R0), np.array([slope(r_) for r_ in RC])
out["bands"] = {"r": [r3(v) for v in R0], "p": [r3(v) for v in p0], "mean_rating": [r3(v) for v in mr],
                "joint_p": r3((np.sum(FC >= F0) + 1) / (len(FC) + 1)),
                "trend": {"slope": r3(s0), "p": r3((np.sum(sC >= s0) + 1) / (len(sC) + 1))},
                "lowest_p": r3(p0.min()), "lowest_band": int(p0.argmin()), "lowest_p_any": r3(np.mean(pC.min(1) <= p0.min()))}
# a 95% upper limit on an excess near tests (every plate near a test scaled by f, which the measure follows, as the
# injection check shows): f such that the shifted calendars, scaled by f, fall at or below the real ratio 5% of the time;
# and the limit expected with no excess, from each calendar in turn as the real one
q05 = np.quantile(C, 0.05)
out["limit"] = {"f95": r3(r_c / q05), "expected": [r3(np.quantile(C / q05, q)) for q in (0.025, 0.16, 0.5, 0.84, 0.975)]}

b_ = np.abs(P.b.values)
out["pointing"] = {"median_low": float(np.median(y_raw[b_ < 10])), "median_high": float(np.median(y_raw[b_ > 40])),
                   "expected_near_over_other": r3(E[m0].mean() / E[~m0].mean()), "slides": len(keep),
                   "near_spring": r3(np.isin(P.date[m0].dt.month, [2, 3, 4, 5, 6]).mean())}

# 6. the same with each transient counted once: VASCO's construction, and ours per night, per photograph, and against a
# pointing curve fitted to the counts, each against the shifted calendars. The ratings rise with distance from the Milky
# Way, so the curve fitted to counts expects fewer near tests than the one fitted to ratings
vc = vcat.groupby("date").size().reindex(vn, fill_value=0).values.astype(float)
values_c = cat.groupby("date").size().reindex(nights, fill_value=0).values.astype(float)
yc = P.n.values.astype(float); Ec = PP.pointing(P)[0]
ratio_c = lambda m, e: (yc[m].sum() / e[m].sum()) / (yc[~m].sum() / e[~m].sum())
cal_p = lambda r, f: r3((np.sum(np.array([f(slid_plates[i]) for i in keep]) >= r) + 1) / (len(keep) + 1))
nohist = lambda d: {k: v for k, v in d.items() if k != "hist"}
r_ph, r_pt = ratio_c(m0, np.ones(len(yc))), ratio_c(m0, Ec)
out["counted"] = {"vasco": nohist(random_nights(vc, near(((vn - DAY0).days).values, (-1, 0, 1), tests=TESTS_V), np.random.default_rng(42))),
                  "per_night": nohist(random_nights(values_c, w, np.random.default_rng(42))),
                  "per_photograph": {"r": r3(r_ph), "p": cal_p(r_ph, lambda m: ratio_c(m, np.ones(len(yc))))},
                  "pointing": {"r": r3(r_pt), "p": cal_p(r_pt, lambda m: ratio_c(m, Ec))},
                  "expected_near_over_other": r3(Ec[m0].mean() / Ec[~m0].mean())}

(ROOT / "data/nuclear.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "hist"} if k == "perm" else v) for k, v in out.items()}, default=str, indent=1))
