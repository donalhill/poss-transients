"""Checks on VASCO-ML's ratings that need no retraining (analysis/ML_PLAN.md, section 2).

2a. The labelled transients: which match the catalogue, which lie on nights near US tests (VASCO's window and list, and
    the full list), and how they were labelled.
2b. The nuclear test sees only plate-level ratings: the weighted excess as the plain excess times the near-test plates'
    mean rating, at each step; whether that mean rating is predicted by plate and night properties alone (a model of
    plate mean rating fitted without the plates near tests), against shifted calendars; and how ratings follow a night's
    transient count.
2c. The shadow with ratings shuffled among transients on each plate (every plate, so the density maps are rebuilt): the
    pooled weighted count and the plate-by-plate weighted ratio, real against shuffled.
2d. Labelled transients on the shadow plates and the near-test plates (their ratings are in-sample).
2e. Sculpting: on plates the shadow never crosses, ratings inside the real shadow outlines against outside, real against
    shuffled within each plate.

python analysis/ml_checks.py -> data/ml_checks.json"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm

import plate_properties as PP
import sample
import shadow_fit as SF
from shadow import HALF, SHADOW_DEG, antisun_axis, cell_counts, deproject, shadow_share, unit, vasco_plate
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
NSHUF = 500
rng = np.random.default_rng(23)
r3 = lambda x: round(float(x), 3)
out = {}

cat, plates = sample.ours()
P = PP.properties(cat, plates)
P["m"] = P.nw / P.n.where(P.n > 0)                                            # each plate's mean rating

# ------------------------------------------------------------------ the tests' calendars, as nuclear.py builds them
DAY0 = pd.Timestamp("1949-01-01")
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
shots["day"] = np.floor((pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601") - DAY0).dt.total_seconds() / 86400).astype(int)
TESTS = np.unique(shots.day.values)
nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
TESTS_V = np.unique(((pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date) - DAY0).dt.days).values)
VWIN = (0, 1, 2)
near = lambda days, off=VWIN, tests=TESTS: np.isin(days, (tests[:, None] + np.array(off)[None, :]).ravel())
day = ((P.date - DAY0).dt.days).values
m0 = near(day)
SHIFTS = [365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-60, 61)]
clear = ~near(day, tuple(range(-14, 15)))
slid = [near(day, tuple(np.array(VWIN) + s)) & clear for s in SHIFTS]
keep = [i for i, m in enumerate(slid) if abs(m.sum() - m0.sum()) <= 0.3 * m0.sum()]

# ------------------------------------------------------------------ 2a, 2d: the labelled transients
lab = sample.labelled()
full = pd.read_csv(ROOT / "data/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv", usecols=["obs_date"])
lab["plate"] = lab.row.map(lambda r: cat.plate.get(r, None) if r >= 0 else None)
arch = pd.to_datetime(full.obs_date.astype(str).str[:10])
lab["near_vasco"] = (lab.row >= 0) & np.isin(((arch.reindex(lab.row.clip(0)).values - np.datetime64(DAY0)) // np.timedelta64(1, "D")),
                                               (TESTS_V[:, None] + np.array([-1, 0, 1])[None, :]).ravel())
lab["near_ours"] = lab.plate.map(pd.Series(m0, P.index)).fillna(False).astype(bool)
grp = lambda m: {"real": int((lab.label[m] == 1).sum()), "flaw": int((lab.label[m] == 0).sum())}
SPL = pd.read_parquet(ROOT / "data/shadow_plates.parquet")
shadow_pl = set(SPL.index[SPL.f_area > 0])
out["labels"] = {"total": grp(lab.index == lab.index), "matched": grp(lab.row >= 0),
                 "by_file": {f: grp(lab.label_file == f) for f in lab.label_file.unique()},
                 "near_tests_vasco": grp(lab.near_vasco), "near_tests_ours": grp(lab.near_ours),
                 "near_tests_by_file": {f: grp(lab.near_vasco & (lab.label_file == f)) for f in lab.label_file.unique()},
                 "on_shadow_plates": grp(lab.plate.isin(shadow_pl)),
                 "share_of_window_transients_labelled": r3(lab.near_vasco.sum() / max(1, near(((arch - DAY0).dt.days).values, (-1, 0, 1), TESTS_V).sum())),
                 "share_of_all_transients_labelled": r3((lab.row >= 0).sum() / len(full)),
                 "mean_rating": {"real": r3(lab.prob[lab.label == 1].mean()), "flaw": r3(lab.prob[lab.label == 0].mean())}}

# what drives the ratings: each transient's leading input by SHAP value, from VASCO-ML's scored file; the inputs that
# describe a whole plate (its quality flag) or a whole night (pooled over that night's plates, common.compute_plate_features)
import zipfile
with zipfile.ZipFile(ROOT / "data/external/scored_catalog.zip") as z, z.open("scored_catalog.csv") as f_:
    top1 = pd.read_csv(f_, usecols=["shap_top1_name"]).shap_top1_name
WHOLE = ["plate_quality_GOOD", "plate_quality_MODERATE", "plate_low_snr_frac", "plate_snr_std", "plate_elongation_mean", "plate_n_high_snr"]
out["drivers"] = {"leading": {k: r3(v) for k, v in top1.value_counts(normalize=True).head(8).items()},
                  "plate_or_night_share": r3(top1.isin(WHOLE).mean()),
                  "rating_variance_between_plates": r3(cat.groupby("plate").w.transform("mean").var() / cat.w.var())}

# ------------------------------------------------------------------ 2b: the nuclear test sees only plate-level ratings
E_n, _, _ = PP.pointing(P)                                                     # expected counts
E_w, _, _ = PP.pointing(P, y="nw")                                             # expected summed ratings
ratio = lambda y, e, m: (y[m].sum() / e[m].sum()) / (y[~m].sum() / e[~m].sum())
one = np.ones(len(P))
steps = {}
for name, e_n, e_w in (("per_photograph", one, one), ("against_pointing", E_n, E_w)):
    r_n, r_w = ratio(P.n.values, e_n, m0), ratio(P.nw.values, e_w, m0)
    steps[name] = {"counts": r3(r_n), "ratings": r3(r_w), "rating_factor": r3(r_w / r_n)}
# the near-test plates' mean rating, against a model of plate mean rating from plate and night properties alone, fitted
# without the plates near tests; scored like the nuclear test, against the shifted calendars
nights = P.groupby("date").n.transform("sum")                                  # the night's transient count
X_ = PP.standardised(P).assign(lnight=np.log(nights.values), ln=np.log(P.n.clip(1).values))
form = "bs(b, df=4) + lebv_plate + ecl_lat + eq_row + quality_bad + depth + fwhm + exp + airmass + lnight + ln"
Xm = np.asarray(patsy.dmatrix(form, X_)); ok = P.n.values > 0


def excess(test, held):
    """The test plates' mean rating over the model's prediction, against the same for the other plates; the model is
    fitted without the plates in held, and neither group holds them except the test plates themselves."""
    u = ok & ~held
    fitm = sm.GLM(P.m.values[u], Xm[u], family=sm.families.Binomial(), freq_weights=P.n.values[u]).fit()
    pred = fitm.predict(Xm)
    rated = lambda m: P.nw.values[m].sum() / (P.n.values * pred)[m].sum()
    return rated(test & ok) / rated(~held & ok), 1 - fitm.deviance / fitm.null_deviance


r_m, dev = excess(m0, m0)
C_m = np.array([excess(slid[i], slid[i] | m0)[0] for i in keep])
steps["mean_rating_near_tests"] = {"actual_over_others": r3((P.nw[m0].sum() / P.n[m0].sum()) / (P.nw[~m0].sum() / P.n[~m0].sum())),
                                   "after_plate_and_night_properties": r3(r_m), "calendars_median": r3(np.median(C_m)),
                                   "calendars_90": [r3(np.quantile(C_m, 0.05)), r3(np.quantile(C_m, 0.95))],
                                   "p_as_high": r3((np.sum(C_m >= r_m) + 1) / (len(C_m) + 1)),
                                   "model_deviance_explained": r3(dev)}
# how ratings follow a night's transient count, raw and after each plate's pointing
nt_ = P.assign(e=E_n).groupby("date")[["n", "nw", "e"]].sum()
rank = lambda a, b: r3(pd.Series(a).rank().corr(pd.Series(b).rank()))
steps["night_count_vs_mean_rating"] = {"raw": rank(nt_.n.values, (nt_.nw / nt_.n).values),
                                       "count_over_pointing": rank((nt_.n / nt_.e).values, (nt_.nw / nt_.n).values)}
out["nuclear"] = steps

# ------------------------------------------------------------------ 2c, 2e: the shadow, ratings shuffled within each plate
codes = pd.factorize(cat.plate)[0]
order = np.argsort(codes, kind="stable")


def shuffled():
    """Ratings shuffled among transients on the same plate."""
    w = np.empty(len(cat)); w[order] = cat.w.values[np.argsort(codes + rng.random(len(codes)))]
    return w


c6 = (cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)
hosts_mask = ~cat.plate.isin(SPL.query("f_area > 0 or f_area_topo > 0").index) & c6
partial = sorted(SPL.query("0 < f_area < 1").index)
geo = {}
for q in partial:
    p = plates.loc[q]; V = vasco_plate(q, p, cat[c6]); s = shadow_share(V["sub"], p.ut)
    ax = np.asarray(antisun_axis(Time(pd.Timestamp(p.ut).to_datetime64()))).reshape(3)
    c = cat[c6 & (cat.plate == q)]
    geo[q] = dict(s=s, idx=c.index.values, inside=(unit(*deproject(c.xi.values, c.eta.values, p.ra_c, p.dec_c)) @ ax) > np.cos(np.radians(SHADOW_DEG)))
whole = sorted(SPL.query("f_area >= 1").index)
hidx = cat.index[hosts_mask].values; hpl = cat.plate[hosts_mask].values; hx, hy = cat.xi[hosts_mask].values, cat.eta[hosts_mask].values
hcode = pd.factorize(hpl)[0]


def shadow_stats(w):
    """The pooled weighted count in shadow (every crossed plate) and the plate-by-plate weighted ratio, for ratings w."""
    ws = pd.Series(w, cat.index)
    g = int(round(2 * HALF / 0.1)); maps = np.zeros((hcode.max() + 1, g * g))
    ix = np.floor((hx + HALF) / 0.1).astype(int); iy = np.floor((hy + HALF) / 0.1).astype(int)
    np.add.at(maps, (hcode, iy * g + ix), ws.reindex(hidx).values)
    acc = (maps / maps.sum(1, keepdims=True)).mean(0); acc /= acc.mean()
    f = np.array([np.sum(geo[q]["s"] * acc) / acc.sum() for q in partial])
    O = np.array([ws[geo[q]["idx"]][geo[q]["inside"]].sum() for q in partial]); n = np.array([ws[geo[q]["idx"]].sum() for q in partial])
    k = np.array([n[j] / (ws[geo[q]["idx"]] ** 2).sum() for j, q in enumerate(partial)])
    r = float(SF.make_fit(f)((O * k)[None], (n * k)[None])[0])
    pooled = O.sum() + sum(ws[cat.index[c6 & (cat.plate == q)]].sum() for q in whole)
    return pooled, r


real = shadow_stats(cat.w.values)
sh = np.array([shadow_stats(shuffled()) for _ in range(NSHUF)])
w_nolab = cat.w.copy(); w_nolab[lab.row[(lab.row >= 0) & lab.plate.isin(shadow_pl)].values] = 0.0   # 2d: weight 0 removes them
nolab = shadow_stats(w_nolab.values)
out["shadow"] = {"without_labelled": {"pooled": r3(nolab[0]), "plate_by_plate_r": r3(nolab[1]), "removed": int(lab.plate.isin(shadow_pl).sum())},
                 "pooled": {"real": r3(real[0]), "shuffled_median": r3(np.median(sh[:, 0])),
                            "shuffled_90": [r3(np.quantile(sh[:, 0], 0.05)), r3(np.quantile(sh[:, 0], 0.95))],
                            "p_as_low": r3((np.sum(sh[:, 0] <= real[0]) + 1) / (NSHUF + 1))},
                 "plate_by_plate_r": {"real": r3(real[1]), "shuffled_median": r3(np.median(sh[:, 1])),
                                      "shuffled_90": [r3(np.quantile(sh[:, 1], 0.05)), r3(np.quantile(sh[:, 1], 0.95))],
                                      "p_as_low": r3((np.sum(sh[:, 1] <= real[1]) + 1) / (NSHUF + 1))},
                 "shuffles": NSHUF}

# 2e: on plates the shadow never crosses, each real outline laid at its place on the plate face: the mean rating of host
# transients inside it against outside, pooled over every host plate; shuffling within plates keeps each plate's ratings
# and removes any tie to position
cats = cat[hosts_mask]


def sculpt(w):
    ws = pd.Series(w, cat.index).reindex(cats.index).values
    res = []
    for q in partial:
        s = geo[q]["s"]
        g = int(round(2 * HALF / 0.1)); ix = np.floor((cats.xi.values + HALF) / 0.1).astype(int); iy = np.floor((cats.eta.values + HALF) / 0.1).astype(int)
        inside = s[iy * g + ix] > 0.5
        res.append([ws[inside].mean(), ws[~inside].mean()])
    a = np.array(res)
    return float(np.mean(a[:, 0] / a[:, 1]))


s_real = sculpt(cat.w.values)
s_null = np.array([sculpt(shuffled()) for _ in range(200)])
out["sculpting"] = {"inside_over_outside_mean_rating": r3(s_real), "shuffled_median": r3(np.median(s_null)),
                    "shuffled_90": [r3(np.quantile(s_null, 0.05)), r3(np.quantile(s_null, 0.95))],
                    "p_two_sided": r3(2 * min(np.mean(s_null >= s_real), np.mean(s_null <= s_real)))}

(ROOT / "data/ml_checks.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
