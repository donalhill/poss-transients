"""Build docs/data/site.json and docs/img/* for the walkthrough page from the analysis outputs.
Run from the repo root after the analysis scripts: python analysis/build_site.py"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2, norm
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import astropy.units as u
from astropy.coordinates import SkyCoord, get_sun
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from geometry import deproject, tangent                                      # noqa: E402
OUT = ROOT / "docs"
(OUT / "data").mkdir(parents=True, exist_ok=True)
(OUT / "img").mkdir(parents=True, exist_ok=True)
site = {}


r1 = lambda x: round(float(x), 1)
r3 = lambda x: round(float(x), 3)

# ------------------------------------------------------------------ sky map + shadow plates
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402
cat, pl = sample.ours()
pl["N"] = cat.groupby("plate").size().reindex(pl.index).fillna(0).astype(int)
vasco_plates = set(sample.vasco()[1].index)
SP = json.loads((ROOT / "data/shadow_pooled.json").read_text())
pt = pd.read_parquet(ROOT / "data/shadow_plates.parquet")
avg = SP["ours_geo"]["avg_per_plate"]
edge = np.linspace(-SH.HALF, SH.HALF, 5)
sq = np.r_[np.c_[edge, np.full(5, -SH.HALF)], np.c_[np.full(5, SH.HALF), edge][1:], np.c_[edge[::-1], np.full(5, SH.HALF)][1:],
           np.c_[np.full(5, -SH.HALF), edge[::-1]][1:-1]]
plates = []
for name, p in pl.iterrows():
    ra, dec = deproject(sq[:, 0], sq[:, 1], p.ra_c, p.dec_c)
    rec = {"id": name, "N": int(p.N), "ra": r1(p.ra_c), "dec": r1(p.dec_c), "date": p.ut.strftime("%Y-%m-%d"),
           "poly": [[r1(a), r1(b)] for a, b in zip(ra, dec)], "vasco": name in vasco_plates}
    if pt.loc[name, "f_area"] > 0:
        s = get_sun(Time(p.ut))
        ara, adec = (s.ra.deg + 180) % 360, -s.dec.deg
        # small circle of radius 8.5 deg about the anti-sun
        circ = SkyCoord(ara * u.deg, adec * u.deg).directional_offset_by(np.radians(np.arange(0, 361, 6)) * u.rad, SH.SHADOW_DEG * u.deg)
        rec.update({"frac": r3(pt.loc[name, "f_area"]), "obs": int(pt.loc[name, "n_in"]), "exp": r1(pt.loc[name, "f_area"] * avg),
                    "circle": [[r1(a), r1(b)] for a, b in zip(circ.ra.deg, circ.dec.deg)]})
    plates.append(rec)
site["plates"] = plates

# ------------------------------------------------------------------ the plate-level measurement (shadow_fit.py): per-plate
# counts and ratios, the common fit, its calibration on toys from real plates, and the cross-checks; and, for each plate
# partly in shadow, every transient in its central 6 x 6 deg (xi east, eta north, deg), which side of the shadow's edge
# it lies, and the edge itself
F_ = json.loads((ROOT / "data/shadow_fit.json").read_text())
site["shadow"] = {k: F_[k] for k in ("fit", "closure", "whole")}
site["shadow"]["toys"] = {k: F_["toys"][k] for k in ("null", "p_as_high", "hist_null", "interval68", "interval90")}
site["shadow"]["injected"] = F_["injected"]                              # a shortfall as large as VASCO's, put into the same trials
site["shadow"]["checks"] = F_["checks"]
plc = {p_["plate"]: p_ for p_ in F_["plates"]}
up_close = []
for q in sorted(plc, key=lambda q: -plc[q]["shadow"]):
    p_ = pl.loc[q]; t_ = Time(pd.Timestamp(p_.ut).to_datetime64())
    c_ = cat[(cat.plate == q) & (cat.xi.abs() < SH.HALF) & (cat.eta.abs() < SH.HALF)]
    # the shadow's edge: the 8.5 deg circle about the anti-sun point at mid-exposure, on the plate's tangent plane
    ax = np.asarray(SH.antisun_axis(t_)).reshape(3)
    c0 = SkyCoord(np.degrees(np.arctan2(ax[1], ax[0])) * u.deg, np.degrees(np.arcsin(ax[2])) * u.deg)
    ring_ = c0.directional_offset_by(np.radians(np.arange(0, 360.5, 1.5)) * u.rad, SH.SHADOW_DEG * u.deg)
    ex, ey = tangent(ring_.ra.deg, ring_.dec.deg, p_.ra_c, p_.dec_c)
    dark = SH.antisun_geocentric(c_.ra.values, c_.dec.values, Time(np.full(len(c_), pd.Timestamp(p_.ut).to_datetime64()))) < SH.SHADOW_DEG
    up_close.append({**plc[q], "sunlit": plc[q]["n"] - plc[q]["O"], "xi": [round(float(v_), 3) for v_ in c_.xi],
                     "eta": [round(float(v_), 3) for v_ in c_.eta], "dark": "".join("1" if d_ else "0" for d_ in dark),
                     "edge": [[round(float(x_), 3), round(float(y_), 3)] for x_, y_ in zip(ex, ey)]})
site["shadow"]["plates"] = up_close
site["shadow"]["half"] = SH.HALF

# ------------------------------------------------------------------ two ideas of chance (shadow_pooled.py)
sums = np.load(ROOT / "data/shadow_paste_scaled.npy")                       # scaled to each shadow plate's pointing
bins = np.arange(0, 2001, 25)
site["paste"] = {"bins": bins.tolist(), "counts": np.histogram(sums, bins)[0].tolist(), "n": int(len(sums)),
                 "observed": SP["ours_geo"]["observed"], "area_expected": SP["ours_geo"]["expected"],
                 "p": SP["ours_geo"]["pointing"]["p"], "binom_sd": SP["ours_geo"]["binom_sd"], "binom_z": abs(SP["ours_geo"]["binom_z"]),
                 "sig": round(float(norm.isf(SP["ours_geo"]["pointing"]["p"])), 1)}

# ------------------------------------------------------------------ nuclear tests (nuclear.py)
NU = json.loads((ROOT / "data/nuclear.json").read_text())
site["perm"] = NU["perm"]
site["perm"]["corrected"] = NU["corrected"]

# ------------------------------------------------------------------ alignments
rb = lambda a: a.reshape(*a.shape[:-1], -1, 2).sum(-1)                 # 0.2-dex bins
# VASCO's own score on the released catalogue (alignments.py): every track candidate scored with their formula,
# arbitrated within each size, real catalogue against copies that keep the catalogue's clumps
AZ = np.load(ROOT / "data/alignments.npz")
rng = np.random.default_rng(3)
pear = []
for th in range(0, 91):
    k = 0
    for _ in range(800):
        s = np.sort(rng.uniform(0, 300, 4))
        x = s * np.cos(np.radians(th)) + rng.normal(0, 1, 4); y = s * np.sin(np.radians(th)) + rng.normal(0, 1, 4)
        k += abs(np.corrcoef(x, y)[0, 1]) > 0.99
    pear.append(r3(k / 800))
site["pearson"] = pear
# the same perfect lines through a straight-line (track) fit: kept if the hits sit within 2" (RMS per degree of freedom)
# of the fitted line, which does not depend on the line's direction
rng = np.random.default_rng(3)
fit_keep = []
for th in range(0, 91):
    k = 0
    for _ in range(800):
        s_ = np.sort(rng.uniform(0, 300, 4))
        x = s_ * np.cos(np.radians(th)) + rng.normal(0, 1, 4); y = s_ * np.sin(np.radians(th)) + rng.normal(0, 1, 4)
        ev = np.linalg.eigvalsh(np.cov(np.c_[x, y].T, bias=True))
        k += 4 * ev[0] / 2 <= 4.0
    fit_keep.append(r3(k / 800))
site["trackfit"] = fit_keep
# alignments as track reconstruction (analysis/alignments.py): distribution of the track score (the expected number of
# chance tracks as long and as straight, from the density around each candidate; straight to 1" precision; any
# direction) for track candidates of 3-6 transients after arbitration, real catalogue against 1,000 copies with every
# real line broken and the clumps kept; and the tail, with the chance of so large a difference in any of its 8 cells
site["trackdist"] = {"lo": [r1(v_) for v_ in AZ["bins"][:-1].reshape(-1, 2)[:, 0]], "all": []}
for k, r in enumerate((3, 4, 5, 6)):
    B = rb(AZ["null_l"][:, k])
    n = rb(AZ["real_l"][k])
    # the real count's place among the copies, as standard deviations: mid-rank, so a count equal to every copy's
    # (0 against 0 included) sits at 0, and a count beyond every copy at the most 1,000 copies can resolve
    z = norm.ppf(((B < n).sum(0) + 0.5 * (B == n).sum(0) + 0.5) / (len(B) + 1))
    site["trackdist"]["all"].append({"r": r, "real": n.tolist(), "bmean": B.mean(0).round(2).tolist(), "bsd": B.std(0).round(2).tolist(),
                                     "z": np.round(z, 2).tolist()})
cells = [(k, int(round((lvl - AZ["bins"][0]) / 0.1))) for k in range(4) for lvl in (-3, -4)]
Rc = np.array([AZ["real_l"][k, :c].sum() for k, c in cells]); Nc = np.stack([AZ["null_l"][:, k, :c].sum(1) for k, c in cells], 1)
ploc = lambda x: min(np.mean(Nc[:, j] >= x[j]) for j in range(len(cells)))
p_real = ploc(Rc)
site["trackdist"]["tail"] = {"cells": [{"r": (3, 4, 5, 6)[k], "lvl": lvl, "real": int(Rc[i]), "copies": r1(Nc[:, i].mean()), "p": r3(np.mean(Nc[:, i] >= Rc[i]))}
                                       for i, ((k, _), lvl) in enumerate(zip(cells, [-3, -4] * 4))],
                             "p_local": r3(p_real), "p_global": r3(np.mean([ploc(Nc[i]) <= p_real for i in range(len(Nc))]))}

# the track score against VASCO's score (alignments.py): candidate by candidate where both score, in 0.2-dex bins; and the
# share of injected tracks each recovers, by direction, at 3 false alarms per survey (4-transient candidates)
AJ = json.loads((ROOT / "data/alignments.json").read_text())
H2d = AZ["ml_hist"].reshape(40, 2, 40, 2).sum((1, 3))
site["compare"] = {"lo": [r1(v_) for v_ in AZ["bins"][:-1].reshape(-1, 2)[:, 0]], "hist": H2d.tolist(), "rho": r3(float(AZ["ml_rho"])),
                   "top50": int(AZ["ml_top50"]), "counts": AZ["ml_counts"].tolist(), "vasco_strong": AZ["ml_vasco_strong"].tolist(),
                   "vasco_test": AJ["counts"]["candidates_1e-3"]}
SIGS = [s_.split(":") for s_ in AZ["signals"]]
# one real track candidate scoring 1 in 1,000 or better, to show how the track score works: its transients and the
# others around it (arcmin from the box centre, east to the left), each transient's distance from the fitted line, and
# what goes into the score; the score is recomputed from these and must match the search's
EX = next(c_ for c_ in AJ["track"]["candidates"] if c_["plate"] == "XE011")
P_ = np.c_[EX["x"], EX["y"]]; ctr_ = P_.mean(0); Vt_ = np.linalg.svd(P_ - ctr_)[2]
al_, pp_ = (P_ - ctr_) @ Vt_[0], (P_ - ctr_) @ Vt_[1]; L_ = float(np.hypot(*(P_[:, None] - P_[None]).transpose(2, 0, 1)).max()); side_ = max(600.0, 1.5 * L_)
mid_ = (P_[al_.argmin()] + P_[al_.argmax()]) / 2
g_ = cat[cat.plate == EX["plate"]]; Q_ = np.c_[g_.xi * 3600, g_.eta * 3600] - mid_
N_ = int(((np.abs(Q_) <= side_ / 2).all(1)).sum()); n_ = N_ / (side_ / 60) ** 2; w_ = float(np.abs(pp_).max()); r_ = EX["r"]
lam_ = np.pi * N_ * n_ * (2 * n_ * max(w_, 1.0) / 60) ** (r_ - 2) * (L_ / 60) ** r_ / (r_ * math.factorial(r_ - 2))
assert abs(lam_ / EX["score"] - 1) < 0.01, (lam_, EX["score"])
chi_ = float((pp_ ** 2).sum()); assert abs(chi_ / (r_ - 2) / EX["chi2ndof"] - 1) < 0.01   # sigma = 1"
view_ = (np.abs(Q_) <= 0.8 * side_).all(1)
import alignments as AL                                                      # noqa: E402  the search's settings, for the figure
site["method"] = {"plate": EX["plate"], "side": round(side_ / 60, 3), "r": r_, "N": N_, "L": round(L_ / 60, 2), "w": round(w_, 2), "score": EX["score"], "chi2": round(chi_, 2), "p_chi2": round(float(chi2.sf(chi_, r_ - 2)), 3),
                  "sigma": 1.0, "hmax": AL.HMAX, "occupied": AL.BOX, "lmax": round(AL.LMAX / 60), "p_min": round(float(chi2.sf(AL.CHI2MAX[r_ - 2], r_ - 2)), 3),
                  "members": np.round(np.c_[-(P_[:, 0] - mid_[0]), P_[:, 1] - mid_[1]] / 60, 3).tolist(),
                  "others": np.round(np.c_[-Q_[view_, 0], Q_[view_, 1]] / 60, 3).tolist(),
                  "res": np.round(np.c_[(al_ - al_.min()) / 60, pp_], 3).tolist()}


def recovered(key, sigs, k=1, Bt=3):
    n = np.cumsum(AZ[f"null_{key}"][:, k].mean(0)); i_ = np.flatnonzero(n <= Bt)[-1]
    return r1(np.mean([np.cumsum(AZ[f"signal_{key}"][s_, k] - AZ[f"real_{key}"][k])[i_] for s_ in sigs]))


angles = sorted({float(s_[3]) for s_ in SIGS if s_[3] != "any"})
site["recovery"] = {"angles": angles, "false_alarms": 3,
                    "track": [recovered("l", [j for j, s_ in enumerate(SIGS) if s_[3] != "any" and float(s_[3]) == a_]) for a_ in angles],
                    "vasco": [recovered("v", [j for j, s_ in enumerate(SIGS) if s_[3] != "any" and float(s_[3]) == a_]) for a_ in angles],
                    "any": {"track": recovered("l", [j for j, s_ in enumerate(SIGS) if s_[1] == "1.0" and s_[3] == "any"]),
                            "vasco": recovered("v", [j for j, s_ in enumerate(SIGS) if s_[1] == "1.0" and s_[3] == "any"])}}
# ------------------------------------------------------------------ dates
t = pd.read_csv(ROOT / "data/plate_date_corrections.csv", index_col=0, parse_dates=["header_date", "corrected_date"])
used = t.index.isin(pl.index) | t.index.isin(pd.read_parquet(ROOT / "data/plates.parquet").index)
an = t[t.index.isin(pd.read_parquet(ROOT / "data/plates.parquet").index)]
# the plates whose date breaks the rule of their run of serial numbers
exc = [{"plate": p, "e": int(r.label), "header": r.header_date.strftime("%Y-%m-%d"), "correct": r.corrected_date.strftime("%Y-%m-%d"),
        "days": int(r.offset_days)} for p, r in t[t.offset_days.notna() & t.rule_offset.notna() & (t.offset_days != t.rule_offset)].sort_values("label").iterrows()]
# Palomar's plate log against the evidence that does not depend on it: asteroid trails and labels read off the scans
direct = t.label_offset.fillna(t.asteroid_offset)
chk = t.log_offset.notna() & direct.notna()
site["date_summary"] = {"log_matched": int(t.log_offset.notna().sum()), "log_checked": int(chk.sum()),
                        "log_agree": int((t.log_offset[chk] == direct[chk]).sum()),
                        "analysis_plates": int(len(an)), "early": int((an.offset_days == 1).sum()), "correct": int((an.offset_days == 0).sum()),
                        "gross": int((an.offset_days.abs() > 1).sum()), "exceptions": exc}

# ------------------------------------------------------------------ catalogue make-up
ctrl = pd.read_csv(ROOT / "data/external/XE_final_random_sources.csv")
p_all = pd.read_parquet(ROOT / "data/plates.parquet")
ctrl = ctrl.join(p_all[["ra_c", "dec_c"]], on="Name").dropna()
cx, cy = tangent(ctrl.RA.values, ctrl.Dec.values, ctrl.ra_c.values, ctrl.dec_c.values)
c_r = np.maximum(np.abs(cx), np.abs(cy))
t_r = np.maximum(cat.xi.abs(), cat.eta.abs()).values
ann = np.arange(0, 3.01, 0.25)
area = (2 * ann[1:]) ** 2 - (2 * ann[:-1]) ** 2
nt_, _ = np.histogram(t_r, ann); nc_, _ = np.histogram(c_r, ann)
dt, dc = nt_ / area, nc_ / area
ref = (ann[:-1] >= 0.5) & (ann[:-1] < 1.0)
per_plate = np.array([np.histogram(t_r[cat.plate.values == q], ann)[0] for q in cat.plate.unique() if q != ""]) / area
rngb = np.random.default_rng(1)
boot = []
for _ in range(2000):
    sm = per_plate[rngb.integers(0, len(per_plate), len(per_plate))].sum(0); boot.append(sm / sm[ref].mean())
# dust towards each plate (plate_dust.py): away from the Milky Way (|b| > 20 deg), transients on the dustiest plates
# (E(B-V) >= 0.15) against the least dusty (E(B-V) < 0.05)
du = pd.read_csv(ROOT / "data/plate_dust.csv").set_index("plate")
du["b"] = np.abs(SkyCoord(pl.loc[du.index].ra_c.values * u.deg, pl.loc[du.index].dec_c.values * u.deg).galactic.b.deg)
du = du[du.b > 20]; dusty, clear = du[du.ebv >= 0.15], du[du.ebv < 0.05]
dust = {"n_dusty": len(dusty), "median_dusty": int(dusty.n.median()), "n_clear": len(clear), "median_clear": int(clear.n.median())}
# share of plates whose band 2.0-2.75 deg from the centre is denser than their inner 1.5 deg
cs_ = cat; rs_ = np.maximum(cs_.xi.abs(), cs_.eta.abs())
share_edge = float(np.mean([((g < 2.75) & (g >= 2.0)).sum() / (5.5 ** 2 - 4.0 ** 2) > (g < 1.5).sum() / 3.0 ** 2 for _, g in rs_.groupby(cs_.plate)]))
site["radial"] = {"share_edge": r3(share_edge), "r": (ann[:-1] + 0.125).round(3).tolist(), "transients": (dt / dt[ref].mean()).round(3).tolist(),
                  "lo": np.percentile(boot, 5, 0).round(3).tolist(), "hi": np.percentile(boot, 95, 0).round(3).tolist(),
                  "control": (dc / dc[ref].mean()).round(3).tolist()}
# the same on the plate itself: density in 0.2 deg cells (north up, east left), each plate weighted equally, over the
# centre density (within 1 deg); VASCO's random positions alongside
G2 = np.arange(-3.4, 3.41, 0.2)
def plate_map(x, y, w):
    h = np.histogram2d(y, x, [G2, G2], weights=w)[0]
    ctr = (np.abs(G2[:-1] + 0.1)[:, None] < 1) & (np.abs(G2[:-1] + 0.1)[None, :] < 1)
    return (h / h[ctr].mean()).round(3)
w_ = 1 / cs_.groupby("plate").plate.transform("size").values
cc_ = np.isin(ctrl.Name.values, cs_.plate.unique())
site["radial"]["map"] = {"edges": G2.round(2).tolist(), "transients": plate_map(cs_.xi.values, cs_.eta.values, w_).tolist(),
                         "control": plate_map(cx[cc_], cy[cc_], None).tolist()}
cat_all = cat[cat.plate != ""]
pa = p_all.copy()
pa["N"] = cat_all.groupby("plate").size().reindex(pa.index).fillna(0).astype(int)
gal = SkyCoord(pa.ra_c.values * u.deg, pa.dec_c.values * u.deg).galactic
pa["absb"], pa["l"] = np.abs(gal.b.deg), gal.l.deg
pa = pa[pa.N > 0]
import plate_properties as PP                                                # noqa: E402
# the pointing curve (plate_properties.pointing) fitted to every plate: each plate's count expected from its pointing,
# and the middle half of plates about it (actual over expected)
P_pc = PP.properties(*sample.ours()); E_pc = pd.Series(PP.pointing(P_pc)[0], P_pc.index)
lr_ = np.log(P_pc.n[P_pc.n > 0] / E_pc[P_pc.n > 0])
# the equatorial row: each plate's count over what pointing predicts without the row term, by the declination of its
# centre, with each row of plates' median
E_norow = pd.Series(PP.pointing(P_pc, formula=PP.POINTING.replace(" + eq_row", ""))[0], P_pc.index)
_, pl_s = sample.ours(); dec_ = pl_s.dec_c.reindex(P_pc.index)
rr_ = (P_pc.n / E_norow); row_ = (np.round(dec_ / 6) * 6).astype(int)
site["rows"] = {"plates": [{"id": q, "dec": r1(dec_[q]), "r": r3(max(float(rr_[q]), 0.01))} for q in P_pc.index],
                "rows": []}
rng_rows = np.random.default_rng(9)
for k in sorted(row_.unique()):                                               # each row's median, and its spread from resampling the row's plates
    v_ = rr_[row_ == k].values; bm_ = np.median(v_[rng_rows.integers(0, len(v_), (4000, len(v_)))], 1)
    site["rows"]["rows"].append({"dec": int(k), "median": r3(float(np.median(v_))), "lo": r3(float(np.percentile(bm_, 16))), "hi": r3(float(np.percentile(bm_, 84))), "n": int(len(v_))})
eqm_ = rr_[row_ == 0].median(); othm_ = rr_[(row_ > 0)].median()
# the same sky seen on both rows: where the equatorial plates overlap the row to their north (dec +2.25 to +3), each row's
# transients per random point (VASCO's random positions measure the area each plate searched)
from geometry import deproject as dep_
c_s, _ = sample.ours(); pc_ = pl_s.reindex(c_s.plate); dec_t = dep_(c_s.xi.values, c_s.eta.values, pc_.ra_c.values, pc_.dec_c.values)[1]
rnd_ = pd.read_csv(ROOT / "data/external/XE_final_random_sources.csv", usecols=["Dec", "Name"])
def per_random(S):
    nt_ = ((dec_t > 2.25) & (dec_t < 3.0) & c_s.plate.isin(S).values).sum(); nr_ = ((rnd_.Dec > 2.25) & (rnd_.Dec < 3.0) & rnd_.Name.isin(S)).sum()
    return nt_ / nr_
same_sky = per_random(set(row_.index[row_ == 0])) / per_random(set(row_.index[row_ == 6]))
site["crowding"] = [{"id": i, "b": r1(r.absb), "l": r1(r.l), "dec": r1(r.dec_c), "N": int(r.N), "E": r1(E_pc.get(i, np.nan))} for i, r in pa.iterrows()]
# calibration: plates in eight groups by predicted count; each group's average prediction and average actual count, with
# the spread of the latter from resampling the group's plates
ord_ = E_pc.sort_values().index; grp_ = np.array_split(ord_, 8); rng_cal = np.random.default_rng(5)
cal_ = []
for g_ in grp_:
    n_ = P_pc.n[g_].values; bs_ = n_[rng_cal.integers(0, len(n_), (4000, len(n_)))].mean(1)
    cal_.append({"E": r1(E_pc[g_].mean()), "N": r1(n_.mean()), "lo": r1(np.percentile(bs_, 16)), "hi": r1(np.percentile(bs_, 84)), "n": int(len(n_))})
site["pointing_fit"] = {"q25": r3(float(np.exp(np.percentile(lr_, 25)))), "q75": r3(float(np.exp(np.percentile(lr_, 75)))), "calibration": cal_,
                       "share": r3(float(np.corrcoef(np.log(P_pc.n[P_pc.n > 0]), np.log(E_pc[P_pc.n > 0]))[0, 1] ** 2))}   # of the spread in log counts
# transients per plate against distance from the Milky Way: a straight-line fit of log(count) on a smooth curve of that
# distance, so the curve is the typical (median) plate and the band the middle half of plates it predicts
import patsy, statsmodels.api as sm
bnd_ = (0.0, float(np.ceil(pa.absb.max())))
Xc_ = patsy.dmatrix(f"bs(absb, df=4, lower_bound={bnd_[0]}, upper_bound={bnd_[1]})", pa, return_type="dataframe")
fc_ = sm.OLS(np.log(pa.N.values), Xc_).fit(); res_ = np.log(pa.N.values) - fc_.fittedvalues.values
sd_ = 1.4826 * float(np.median(np.abs(res_ - np.median(res_))))         # robust to a few extreme plates
bg_ = np.arange(0.0, bnd_[1] + 0.01, 1.0)
lg_ = np.asarray(fc_.predict(patsy.build_design_matrices([Xc_.design_info], {"absb": bg_})[0]))
site["crowd_fit"] = {"b": bg_.round(1).tolist(), "mid": np.exp(lg_).round(1).tolist(), "q25": np.exp(lg_ - 0.6745 * sd_).round(1).tolist(),
                     "q75": np.exp(lg_ + 0.6745 * sd_).round(1).tolist(), "inside": r3(float(np.mean(np.abs(np.log(pa.N.values) - fc_.fittedvalues.values) <= 0.6745 * sd_)))}

# ------------------------------------------------------------------ numbers quoted in the text ({n:key} in site/copy.md)
WORDS = "zero one two three four five six seven eight nine ten eleven twelve".split()
word = lambda n: WORDS[int(n)] if 0 <= int(n) < len(WORDS) else f"{int(n):,}"
odds = lambda p_: word(round(1 / p_)) if round(1 / p_) <= 10 else f"{round(1 / p_):,}"    # "one in {odds}"
num = lambda v, d=0: f"{v:,.{d}f}"
def frac(p_):
    """A share in words: "one in seven" when small, otherwise the nearest simple fraction ("two in five")."""
    if p_ < 0.25:
        return f"one in {odds(p_)}"
    k, n = min(((k, n) for n in range(2, 6) for k in range(1, n)), key=lambda kn: abs(kn[0] / kn[1] - p_))
    return f"{word(k)} in {word(n)}"
from scipy.stats import norm                                                 # noqa: E402
sig = lambda p_: f"{norm.isf(p_):.1f}"
N_ = pl.N
avg_all = N_.mean(); sd_all = np.sqrt(avg_all); lo_b, hi_b = avg_all - 2 * sd_all, avg_all + 2 * sd_all
pa_b = pa.set_index(pa.index).absb
du_all = pd.read_csv(ROOT / "data/plate_dust.csv").set_index("plate")
G, T_, Fp = SP["ours_geo"], SP["ours_topo"], {p_["plate"]: p_ for p_ in F_["plates"]}
V_ = SP["vasco_geo"]
ST = NU["steps"]
numbers = {
    # the catalogue
    "cat.transients": num(len(cat)), "cat.plates": num(len(pl)),
    "cat.left_out": num(len(pl) - len(vasco_plates)), "cat.min": num(N_.min()), "cat.max": num(N_.max()),
    "cat.avg": num(avg_all), "cat.noise": num(sd_all), "cat.band_lo": num(lo_b), "cat.band_hi": num(hi_b),
    "cat.in_band": num(((N_ >= lo_b) & (N_ <= hi_b)).sum()),
    "cat.ml_default": num(100 * pd.read_csv(ROOT / "data/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv", usecols=["ml_probability"]).ml_probability.value_counts(normalize=True).iloc[0], 1),
    "nc2.photo_pc": num(100 * NU["corrections"][1]["calendar_null"]["p"]),
    "cat.over1000": word((N_ > 1000).sum()), "cat.edge_share": num(100 * share_edge),
    "crowd.inner": num(pa[pa.absb < 5].N.median()), "crowd.outer": num(round(pa[pa.absb > 40].N.median(), -1)),
    "crowd.ratio": word(round(pa[pa.absb > 40].N.median() / pa[pa.absb < 5].N.median())),
    "dust.n_dusty": num(dust["n_dusty"]), "dust.median_dusty": num(dust["median_dusty"]),
    "dust.n_clear": num(dust["n_clear"]), "cat.eq_plates": num(int((row_ == 0).sum())), "cat.eq_share": num(100 * eqm_ / othm_), "cat.eq_same_sky": num(100 * same_sky), "sh.eq_plates": word(sum(1 for q in pt.index[pt.f_area > 0] if abs(pl.loc[q].dec_c) < 3)), "cat.pointing_share": num(100 * site["pointing_fit"]["share"]), "dust.median_clear": num(dust["median_clear"]),
    # Earth's shadow, pooled (shadow_pooled.py)
    "sh.plates": word(G["plates"]), "sh.plates_n": num(G["plates"]), "sh.partial": word(len(F_["plates"])),
    "sh.never": num(G["hosts"]["plates"]), "sh.observed": num(G["observed"]), "sh.worth": num(G["plates_worth"], 2),
    "sh.avg": num(G["avg_per_plate"], 1), "sh.expected": num(G["expected"]), "sh.binom_sd": num(G["binom_sd"]),
    "sh.binom_z": num(abs(G["binom_z"]), 1), "sh.p_random": num(100 * G["p_random"]),
    "sh.p_point": num(100 * G["pointing"]["p"]), "sh.trials": num(len(sums)), "sh.sig_point": sig(G["pointing"]["p"]),
    "sh.odds_point": odds(G["pointing"]["p"]), "sh.odds_near": odds(G["p_near"]), "sh.sig_near": sig(G["p_near"]), "sh.near_max": num(G["near"]["max_deg"]), "sh.near_hosts": num(G["near"]["hosts"]),
    "sh.shadow_lo": num(G["near"]["shadow_deg"][0]), "sh.shadow_hi": num(round(G["near"]["shadow_deg"][1])),
    "sh.centre_obs": num(G["centre"]["observed"]), "sh.centre_exp": num(G["centre"]["expected"]), "sh.odds_centre": odds(G["centre"]["p_random"]),
    "sh.topo_plates": word(T_["plates"]), "sh.topo_obs": num(T_["observed"]), "sh.topo_exp": num(T_["expected"]),
    "sh.odds_topo": odds(T_["p_random"]), "sh.sig_topo": sig(T_["p_random"]),
    "sh.vasco_obs": num(V_["observed"]), "sh.vasco_exp": num(V_["expected"]), "sh.vasco_plates": word(V_["plates"]), "sh.vasco_plates_cap": word(V_["plates"]).capitalize(),
    # Earth's shadow, plate by plate (shadow_fit.py)
    "fit.r": num(F_["fit"]["r"], 2), "fit.count_lo": num(F_["fit"]["count_lo"], 2), "fit.count_hi": num(F_["fit"]["count_hi"], 2),
    "fit.null_median": num(F_["toys"]["null"]["median"], 2), "fit.null_lo": num(F_["toys"]["null"]["q16"], 2),
    "fit.null_hi": num(F_["toys"]["null"]["q84"], 2), "fit.ci90_lo": num(F_["toys"]["interval90"][0], 2),
    "fit.ci90_hi": num(F_["toys"]["interval90"][1], 2), "fit.inj_r": num(F_["injected"]["r_true"], 2),
    "fit.inj_pct": num(100 * F_["injected"]["p_as_high"]), "fit.inj_cl": num(100 * (1 - F_["injected"]["p_as_high"])),
    "fit.max_plate": max(Fp, key=lambda q: Fp[q]["r"] or 0), "fit.max_r": num(max(p_["r"] or 0 for p_ in F_["plates"]), 1), "fit.trials": num(F_["toys"]["trials"]),
    "fit.ex_plate": "XE358", "fit.ex_area": num(100 * Fp["XE358"]["shadow"], 1), "fit.ex_weighted": num(100 * Fp["XE358"]["weighted"], 1),
    "fit.ex_share": num(100 * Fp["XE358"]["O"] / Fp["XE358"]["n"]), "fit.ex_O": num(Fp["XE358"]["O"]), "fit.ex_n": num(Fp["XE358"]["n"]),
    "pl.XE415_n": num(Fp["XE415"]["n"]), "pl.XE355_n": num(Fp["XE355"]["n"]), "pl.XE617_area": num(100 * Fp["XE617"]["shadow"]),
    "pl.XE617_O": num(Fp["XE617"]["O"]), "pl.XE617_n": num(Fp["XE617"]["n"]), "pl.XE640_area": num(100 * Fp["XE640"]["shadow"]),
    "pl.XE640_n": num(Fp["XE640"]["n"]),
    # nuclear tests (nuclear.py)
    "nuc.vasco_nights": num(NU["vasco"]["nights"]), "nuc.vasco_ratio": num(NU["vasco"]["real"], 2),
    "nuc.vasco_odds": odds(NU["vasco"]["p"]), "nuc.vasco_sig": num(NU["vasco"]["sigma"], 1),
}
south = pd.read_csv(ROOT / "data/shadow_plates_log.csv")
# shadow plates holding under a tenth of what their pointing predicts (pointing curve fitted to every plate)
empty_sh = [q for q in pt.index[pt.f_area > 0] if P_pc.n[q] / E_pc[q] < 0.1]
numbers.update({"sh.empty_plates": " and ".join(empty_sh), "sh.empty_n": word(len(empty_sh)), "sh.taurus_lc": word(sum(1 for q in pt.index[pt.f_area > 0] if 40 < pl.loc[q].ra_c < 75)), "sh.south": num((~south.in_catalogue).sum()),
                "sh.south_ratio": word(round(south[~south.in_catalogue].shadow_fraction.sum() / south[south.in_catalogue].shadow_fraction.sum())),
                })
# the dates (date_corrections.py and its evidence)
DS, DC = site["date_summary"], pd.read_csv(ROOT / "data/date_check_summary.csv")
src = t.source.value_counts()
off_ = t.offset_days.reindex(pl.index); left_ = ~pl.index.isin(vasco_plates)
numbers.update({
    "dates.early_pct": num(5 * round(20 * (off_ == 1).mean())), "dates.ast_late_trails": num(DC[DC.best_day_offset == 1].n_at_best.sum()),
    "dates.ast_late_plates": num((DC.best_day_offset == 1).sum()), "dates.ast_ok_trails": num(DC[DC.best_day_offset == 0].n_at_best.sum()),
    "dates.ast_ok_plates": num((DC.best_day_offset == 0).sum()), "dates.log_agree": num(DS["log_agree"]), "dates.log_checked": num(DS["log_checked"]),
    "dates.dated": num(len(t)), "dates.by_log": num(src["Palomar log"]), "dates.by_asteroids": num(src["asteroids"]),
    "dates.by_label": num(src["label"]), "dates.by_runs": num(src["block rule"] + src.get("serial neighbours", 0)), "dates.runs_cat": (lambda k_: "none of them in the catalogue, so every catalogue plate is dated directly" if k_ == 0 else f"{word(k_)} of them in the catalogue")(int(t[t.source.isin(["block rule", "serial neighbours"]) & t.index.isin(pl.index)].shape[0])), "dates.unresolved": word(src.get("unresolved", 0)),
    "dates.exceptions": word(len(DS["exceptions"])).capitalize(), "dates.gross": word(sum(abs(e_["days"]) > 1 for e_ in DS["exceptions"])).capitalize(),
    "dates.vasco_nights": num(len(sample.vasco_nights())), "dates.dropped": num((left_ & (off_ == 0)).sum()),
    "dates.shifted": num((~left_ & (off_ == 0)).sum())})
# the alignments (alignments.py, alignment_fit.py, vasco_candidates.py)
VC, TK, CM, RC = json.loads((ROOT / "data/vasco_candidates.json").read_text()), AJ["track"], site["compare"], site["recovery"]
i3, i4 = TK["levels"].index(-3.0), TK["levels"].index(-4.0)
# the same null made other ways (alignment_null_checks.py): other move sizes and a Gaussian move, with transients moved
# together as on the page; every transient moved on its own by up to 2' (none of the bunching kept); and the copies of
# the catalogue with artificial tracks added
NCK = json.loads((ROOT / "data/alignment_null_checks.json").read_text())
moves = [r_ for r_ in NCK["results"] if r_["link_arcsec"] == AJ["link_arcsec"]]
alone = next(r_ for r_ in NCK["results"] if r_["alone"])
p3 = [100 * r_["track_share_ge"][i3] for r_ in moves + [alone]]   # every other way of making the random catalogues
ab = NCK["absorption"]; added = [NCK["injected_track"][i] - NCK["real_track"][i] for i in (i3, i4)]
pr_ = NCK["pairs"][0]["copies_over_real"]; i20 = NCK["separation_arcsec"].index(20.0)            # the page's copies, 20-30"
# the closest pair of points in each of VASCO's five shortlisted alignments (from their tables): those closer than the
# plates' two-hit resolution cannot be in the catalogue
VS = pd.read_csv(ROOT / "data/vasco_2025_candidates.csv"); close = []
for _, g_ in VS.groupby("candidate"):
    c_ = SkyCoord(g_.ra_hms.values, g_.dec_dms.values, unit=(u.hourangle, u.deg))
    d_ = min(c_[i_].separation(c_[j_]).arcsec for i_ in range(len(c_)) for j_ in range(i_ + 1, len(c_)))
    if d_ < 8: close.append(d_)
# the fit of the four track score distributions: the background alone, and the interval on the number of glint tracks
AF = json.loads((ROOT / "data/alignment_fit.json").read_text())
site["limit"] = {f_: {"p": v_["p_post"], "s95": v_["s95"]} for f_, v_ in AF["flashes"].items()}
c4 = next(c_ for c_ in site["trackdist"]["tail"]["cells"] if c_["r"] == 4 and c_["lvl"] == -3)
BX = {b_["box_arcsec"]: b_ for b_ in NCK["boxes"]}                    # the track score with smaller boxes
# VASCO's score counted without arbitration: overlapping candidates grouped into families sharing transients
vk_ = [a_ for a_ in AL.score(AL.cat)[5] if a_["mu"] <= 1e-3]; fam_ = list(range(len(vk_))); own_ = {}
root_ = lambda i: i if fam_[i] == i else root_(fam_[i])
for i_, a_ in enumerate(vk_):
    for t_ in a_["ids"]:
        if t_ in own_: fam_[root_(i_)] = root_(own_[t_])
        else: own_[t_] = i_
fams_ = pd.Series([vk_[i_]["plate"] for i_ in range(len(vk_))]).groupby([root_(i_) for i_ in range(len(vk_))]).agg(["first", "size"]).sort_values("size", ascending=False)
AC = AJ["counts"]
numbers.update({
    "align.alt_p3_lo": num(min(p3)), "align.alt_p3_hi": num(max(p3)),
    "align.alone_move": num(alone["shift_arcsec"] / 60), "align.alone_exp3": num(alone["track_mean"][i3], 1),
    "align.vasco_real": num(NCK["real_vasco_1e-3"]), "align.vasco_exp": num(moves[0]["vasco_mean"], 1), "align.vasco_p": num(100 * moves[0]["vasco_share_ge"]),
    "align.vasco_alone": num(alone["vasco_mean"], 1),
    "align.vasco_all": num(AC["all_candidates_1e-3"]["real"]), "align.vasco_all_exp": num(AC["all_candidates_1e-3"]["copies_mean"], 1),
    "align.vasco_fam1": fams_.iloc[0]["first"], "align.vasco_fam1_n": num(fams_.iloc[0]["size"]),
    "align.vasco_fam2": fams_.iloc[1]["first"], "align.vasco_fam2_n": num(fams_.iloc[1]["size"]),
    "align.vasco_straight": num(AC["straight_candidates_1e-3"]["real"]), "align.vasco_straight_exp": num(AC["straight_candidates_1e-3"]["copies_mean"], 1),
    "align.vasco_straight_p": num(100 * AC["straight_candidates_1e-3"]["share_copies_ge_real"]),
    "align.vasco_alone_phrase": (lambda k_: "more than any of them holds" if k_ == 0 else f"more than all but {num(k_)} of them hold")(
        round(alone["vasco_share_ge"] * AJ["ncopies"])),
    "align.eff10": num(BX[10.0]["found_diagonal"]),
    "align.absorb3": num(100 * (ab["copies_of_injected_mean"][i3] - ab["copies_of_plain_mean"][i3]) / added[0]),
    "align.absorb4": num(100 * (ab["copies_of_injected_mean"][i4] - ab["copies_of_plain_mean"][i4]) / added[1]),
    "align.keep3": num(100 - 100 * (ab["copies_of_injected_mean"][i3] - ab["copies_of_plain_mean"][i3]) / added[0]),
    "align.keep4": num(100 - 100 * (ab["copies_of_injected_mean"][i4] - ab["copies_of_plain_mean"][i4]) / added[1]),
    "align.pair2030": num(pr_[i20], 1),
    "align.vasco_close_n": word(len(close)).capitalize(), "align.vasco_close": " and ".join(num(d_) for d_ in sorted(close)),
    "align.fit_bg_p": num(100 * AF["background"]["p"]), **{f"align.lim{f_}": word(v_["s95"]) for f_, v_ in AF["flashes"].items()}, **{f"align.lim{f_}_n": num(v_["s95"]) for f_, v_ in AF["flashes"].items()}, **{f"align.lim{f_}_2": num(v_["s95_scatter2"]) for f_, v_ in AF["flashes"].items()},
    "align.fit_none": "zero" if all(v_["s_hat"] == 0 for v_ in AF["flashes"].values()) else "above zero, though no more than chance allows",
    "align.fit_cls": "the same limits" if all(v_["s95"] == v_["s95_cls"] for v_ in AF["flashes"].values()) else
        "limits of " + ", ".join(word(AF["flashes"][f_]["s95_cls"]) for f_ in ("6", "5", "4")) + f', and {word(AF["flashes"]["3"]["s95_cls"])}'
        + " respectively",
    "align.c4_real": num(c4["real"]), "align.c4_exp": num(c4["copies"], 1), "align.c4_p": num(100 * c4["p"]),
    "align.p_global": num(100 * site["trackdist"]["tail"]["p_global"]),
    "align.ex_plate": site["method"]["plate"]})
numbers.update({
    "align.vasco_transients": num(VC["n_transients"]), "align.vasco_in_cat": word(VC["n_in_catalogue"]),
    "align.copies": num(AJ["ncopies"]), "align.shift": num(AJ["shift_arcsec"]), "align.link": num(AJ["link_arcsec"]),
    "align.n3": word(TK["real"][i3]).capitalize(), "align.n3_lc": word(TK["real"][i3]), "align.exp3": num(TK["copies_mean"][i3], 1), "align.p3": num(100 * TK["share_copies_ge_real"][i3]),
    "align.exp4": num(TK["copies_mean"][i4], 1),
    "align.rho": num(CM["rho"], 2), "align.top50": num(CM["top50"]), "align.vasco_strong": num(CM["vasco_strong"][0]),
    "align.vasco_bent": num(CM["vasco_strong"][1]), "align.both": num(CM["counts"][2]),
    "align.false_alarms": word(RC["false_alarms"]), "align.vasco_ew": num(RC["vasco"][RC["angles"].index(0.0)]),
    "align.track_ew": num(RC["track"][RC["angles"].index(0.0)]), "align.track_diag": num(RC["track"][RC["angles"].index(45.0)]),
    "align.vasco_diag": num(RC["vasco"][RC["angles"].index(45.0)]), "align.vasco_2": num(RC["vasco"][RC["angles"].index(2.0)]),
})
# the nuclear tests (nuclear.py)
pc = lambda p_: num(100 * p_)
ordinal = lambda n: f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"
longdate = lambda d: pd.Timestamp(d).strftime("%-d %B %Y")
numbers.update({"align.pear_diag": num(100 * site["pearson"][45], 1), "align.pear_2": pc(site["pearson"][2]), "align.pear_ew": pc(site["pearson"][0])})
numbers.update({
    "nuc.nights": num(ST["nights_window"]), "nuc.ratio": num(NU["perm"]["raw"]["real"], 2), "nuc.odds": odds(NU["perm"]["raw"]["p"]),
    "nuc.sig": num(NU["perm"]["raw"]["sigma"], 1), "nuc.trials": num(NU["perm"]["raw"]["trials"]), "nuc.plates_window": num(ST["plates_window"], 1), "nuc.plates_other": num(ST["plates_other"], 1),
    **{f"nc2.{k}_r": num(x["r"], 2) for k, x in zip(("night", "photo", "point"), NU["corrections"])},
    "nc2.point_pv": pc(NU["corrections"][2]["vasco_null"]["p"]), "nc2.point_pc": pc(NU["corrections"][2]["calendar_null"]["p"]),
    "nc2.point_pv_sig": num(norm.isf(NU["corrections"][2]["vasco_null"]["p"]), 1),
    "nc2.empty_near": num(100 * ST["empty_window"] / ST["nights_window"]), "nc2.empty_other": num(100 * ST["empty_other"] / ST["nights_other"]),
    "nc2.med_low": num(NU["pointing"]["median_low"]), "nc2.med_high": num(NU["pointing"]["median_high"]),
    "nc2.exp_pct": num(round(100 * (NU["pointing"]["expected_near_over_other"] - 1))), "nc2.spring": pc(NU["pointing"]["near_spring"]),
    "nc2.slides": num(NU["pointing"]["slides"]), "nc2.sig": num(NU["corrected"]["sigma"], 1), "nc2.held_out": num(NU["corrected"]["held_out"], 2), "nc2.inj_f": num(NU["null_checks"]["injected"][0]["factor"], 1), "nc2.inj_r": num(NU["null_checks"]["injected"][0]["measured"], 2),
    "nc2.inj_p": num(100 * NU["null_checks"]["injected"][0]["p"], 1), "nc2.odds": odds(NU["corrected"]["p"]), "nc2.frac": frac(NU["corrected"]["p"]),
    "nc4.windows": num(NU["windows"]["windows"]), "nc4.chance": pc(NU["windows"]["chance_vs_vasco"]),
    "nc4.rw_plate": NU["windows"]["rich_week"]["plate"], "nc4.rw_date": longdate(NU["windows"]["rich_week"]["date"]), "nc4.rw_n": num(NU["windows"]["rich_week"]["n"]),
    "nc4.rw_e": num(round(NU["windows"]["rich_week"]["expected"], -1)), "nc4.rw_days": word(NU["windows"]["rich_week"]["days_from_test"]),
    "nc4.rd_plates": num(NU["windows"]["rich_day_share"]["plates_window"]), "nc4.rd_share": pc(NU["windows"]["rich_day_share"]["window"]),
    "nc4.best_before": word(-NU["windows"]["best"]["window"][0]), "nc4.best_after": word(NU["windows"]["best"]["window"][1]), "nc4.best_r": num(NU["windows"]["best"]["r"], 2),
    "nc4.best_span": (lambda b_, a_: f"from the test day to {word(a_)} days after" if b_ == 0 else f"from {word(b_)} days before a test to {word(a_)} after")(-NU["windows"]["best"]["window"][0], NU["windows"]["best"]["window"][1]),
    "nc4.rd_ratio": num(sum(x["n"] for x in NU["windows"]["rich_day"]) / sum(x["expected"] for x in NU["windows"]["rich_day"]), 1),
    "nc4.best_local": pc(NU["windows"]["best"]["p_local"]), "nc4.best_global": pc(NU["windows"]["best"]["p_global"]), "nc4.vasco_r": num(next(c["r"] for c in NU["windows"]["cells"] if c["before"] == 1 and c["after"] == 1), 2),
    "nuc.missing": word(NU["vasco_list"]["missing_shots"]), "nuc.vlist_r": num(NU["vasco_list"]["r"], 2), "nuc.vlist_pct": pc(NU["vasco_list"]["p"]),
    "nuc.vlist_extra": word(NU["vasco_list"]["extra_plates"]), "nc2.made": num(NU["calendars"]["made"]), "nc2.n0": num(NU["calendars"]["plates_real"]), "nc2.match": num(round(100 * NU["calendars"]["match"])),
    "nc4.wo_median": num(NU["windows"]["without"]["median"], 2), "nc4.wo_vasco": num(NU["windows"]["without"]["vasco"], 2),
    "nc4.rd1": NU["windows"]["rich_day"][0]["plate"], "nc4.rd2": NU["windows"]["rich_day"][1]["plate"], "nc4.r_lo": num(min(c["r"] for c in NU["windows"]["cells"]), 2), "nc4.r_hi": num(max(c["r"] for c in NU["windows"]["cells"]), 2), 
    **{f"nc3.{k}_{f}": v for k, g in (("defs", "What counts as a transient"), ("series", "Test series left out"))
       for f, v in (("lo", num(min(x["r"] for x in NU["checks"] if x["group"] == g), 2)), ("hi", num(max(x["r"] for x in NU["checks"] if x["group"] == g), 2)))},
    "nc3.all": num(next(x["r"] for x in NU["checks"] if x["label"] == "All of them together"), 2)})
site["numbers"] = numbers

(OUT / "data/site.json").write_text(json.dumps(site, separators=(",", ":"), allow_nan=False))
print("site.json", round((OUT / "data/site.json").stat().st_size / 1024), "kB")

# ------------------------------------------------------------------ images (DSS plate material, figure-sized crops)
IMGS = ("label_xe356.png", "label_xe617.png", "asteroid_isabella_wide.png", "asteroid_erigone_wide.png")
if "--images" not in sys.argv and all((OUT / "img" / n).exists() for n in IMGS):
    print("images kept (redraw with --images)"); raise SystemExit
def save(fig, name):
    fig.savefig(OUT / "img" / name, dpi=150, bbox_inches="tight", pad_inches=0.02, facecolor="white"); plt.close(fig)

from fetch_labels import OUT as LABELS, fetch as fetch_label
for plate in ("XE356", "XE617"):                 # the label corner of each scan (fetch_labels.py), cropped to the label
    fetch_label(plate)
    a = np.load(LABELS / f"{plate}.npy").astype(float)[150:500, 100:900]
    fig = plt.figure(figsize=(8, 3.5)); ax = fig.add_axes([0, 0, 1, 1])
    lo, hi = np.percentile(a, [1, 99.7]); ax.imshow(a, origin="lower", cmap="gray_r", vmin=lo, vmax=hi); ax.axis("off")
    save(fig, f"label_{plate.lower()}.png")

from asteroid_clock import plate_header, ephemeris, CUT
# one wide view of each plate holding both predicted positions: the archive's date and one day later
for plate, num, nm, label in (("XE356", 210, "isabella", "(210) Isabella"), ("XE617", 163, "erigone", "(163) Erigone")):
    p = plate_header(plate)
    eph = [ephemeris(num, p["t_hdr"] + pd.Timedelta(days=dd), hours=4) for dd in (0, 1)]
    mid = [e.iloc[len(e) // 2] for e in eph]
    c0, c1 = (SkyCoord(m.ra * u.deg, m.dec * u.deg) for m in mid)
    cen = c0.directional_offset_by(c0.position_angle(c1), c0.separation(c1) / 2)
    wid, hgt = 23.0, 14.5   # the same frame for both plates; each holds both predicted positions with room to spare
    assert c0.separation(c1).arcmin + 4 < wid and abs(c0.dec - c1.dec).to(u.arcmin).value + 4 < hgt
    f = CUT / f"{plate}_frame_{num}.fits"
    if not f.exists():
        url = (f"https://archive.stsci.edu/cgi-bin/dss_search?v=poss1_red&r={cen.ra.deg:.6f}&d={cen.dec.deg:.6f}"
               f"&e=J2000&h={hgt:.1f}&w={wid:.1f}&f=fits&c=none")
        f.write_bytes(requests.get(url, timeout=180).content)
    hd = fits.getheader(f)
    assert hd["REGION"].strip() == plate, (plate, hd["REGION"])
    img, w = fits.getdata(f).astype(float), WCS(hd)
    fig = plt.figure(figsize=(9, 9 * img.shape[0] / img.shape[1])); ax = fig.add_axes([0, 0, 1, 1])
    lo, hi = np.percentile(img, [5, 99.8]); ax.imshow(img, origin="lower", cmap="gray_r", vmin=lo, vmax=hi); ax.axis("off")
    H_, W_ = img.shape
    pts = [w.world_to_pixel(SkyCoord(m.ra * u.deg, m.dec * u.deg)) for m in mid]
    corners = {(0.02, 0.05): (0.17 * W_, 0.26 * H_), (0.70, 0.05): (0.84 * W_, 0.26 * H_),
               (0.02, 0.53): (0.17 * W_, 0.74 * H_), (0.70, 0.53): (0.84 * W_, 0.74 * H_)}
    spot = max(corners, key=lambda k: min(np.hypot(corners[k][0] - px, corners[k][1] - py) for px, py in pts))
    for dd, e in enumerate(eph):
        x, y = w.world_to_pixel(SkyCoord(e.ra.values * u.deg, e.dec.values * u.deg))
        ax.plot(x, y, color="#c0392b", lw=1.4, alpha=0.95, ls=(0, (4, 3)) if dd == 0 else "-")
        xm, ym = pts[dd]
        if dd == 1:   # zoom on the trail
            ins = ax.inset_axes([spot[0], spot[1], 0.28, 0.42])
            ins.imshow(img, origin="lower", cmap="gray_r", vmin=lo, vmax=hi); ins.plot(x, y, color="#c0392b", lw=1.0, alpha=0.7)
            ins.set_xlim(xm - 28, xm + 28); ins.set_ylim(ym - 20, ym + 20); ins.set_xticks([]); ins.set_yticks([])
            for sp in ins.spines.values(): sp.set_color("#c0392b"); sp.set_linewidth(1.2)
            ax.indicate_inset_zoom(ins, edgecolor="#c0392b", alpha=0.9)
    pix = abs(hd["PLTSCALE"] * hd["XPIXELSZ"] / 1000)          # arcsec per pixel
    L = 120 / pix   # 2-arcminute scale bar, in the free corner furthest from the paths and the inset
    free = [k for k in corners if k != spot]
    bx, by = corners[max(free, key=lambda k: min(np.hypot(corners[k][0] - px, corners[k][1] - py) for px, py in pts))]
    bx = W_ * (0.04 if bx < W_ / 2 else 0.96) - (0 if bx < W_ / 2 else L); by = H_ * (0.06 if by < H_ / 2 else 0.94)
    ax.plot([bx, bx + L], [by, by], color="black", lw=2.5, solid_capstyle="butt")
    ax.text(bx + L / 2, by + 0.025 * H_, "2′", ha="center", va="bottom", fontsize=15, color="black",
            path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    save(fig, f"asteroid_{nm}_wide.png")
print("images written")
