"""Build docs/data/site.json and docs/img/* for the walkthrough page from the analysis outputs.
Run from the repo root after the analysis scripts: python analysis/build_site.py"""
import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np
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
OUT = ROOT / "docs"
(OUT / "data").mkdir(parents=True, exist_ok=True)
(OUT / "img").mkdir(parents=True, exist_ok=True)
site = {}


def quiet_exec(path, stop=None):
    g = {"__file__": str(path), "__name__": "build"}
    src = Path(path).read_text()
    if stop:
        src = src.split(stop)[0]
    with contextlib.redirect_stdout(io.StringIO()):
        exec(src, g)
    return g


def deproject(xi, eta, ra0, dec0):
    xi, eta, ra0, dec0 = map(np.radians, (xi, eta, ra0, dec0))
    rho = np.hypot(xi, eta); c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def tangent(ra, dec, ra0, dec0):
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    c = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    return (np.degrees(np.cos(dec) * np.sin(ra - ra0) / c),
            np.degrees((np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / c))


r1 = lambda x: round(float(x), 1)
r3 = lambda x: round(float(x), 3)

# ------------------------------------------------------------------ sky map + shadow plates
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
strict = set(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv").obs_date.astype(str).str[:10])
cs = cat[cat.obs_date.astype(str).str[:10].isin(strict)]
pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
pl = pl[pl.index.isin(cs.plate.unique())].copy()
pl["N"] = cs.groupby("plate").size().reindex(pl.index).fillna(0).astype(int)
pt = pd.read_parquet(ROOT / "data/plate_shadow_geo_corr.parquet")
avg = pt.N.sum() / len(pt)
edge = np.linspace(-3.0, 3.0, 5)
sq = np.r_[np.c_[edge, np.full(5, -3.0)], np.c_[np.full(5, 3.0), edge][1:], np.c_[edge[::-1], np.full(5, 3.0)][1:], np.c_[np.full(5, -3.0), edge[::-1]][1:-1]]
plates = []
for name, p in pl.iterrows():
    ra, dec = deproject(sq[:, 0], sq[:, 1], p.ra_c, p.dec_c)
    rec = {"id": name, "N": int(p.N), "ra": r1(p.ra_c), "dec": r1(p.dec_c), "date": p.ut.strftime("%Y-%m-%d"),
           "poly": [[r1(a), r1(b)] for a, b in zip(ra, dec)]}
    if name in pt.index and pt.loc[name, "f_area"] > 0:
        s = get_sun(Time(p.ut))
        ara, adec = (s.ra.deg + 180) % 360, -s.dec.deg
        th = np.radians(np.arange(0, 361, 6))
        # small circle of radius 8.5 deg about the anti-sun
        c0 = SkyCoord(ara * u.deg, adec * u.deg)
        circ = c0.directional_offset_by(th * u.rad, 8.5 * u.deg)
        rec.update({"frac": r3(pt.loc[name, "f_area"]), "obs": int(pt.loc[name, "n_in"]),
                    "exp": r1(pt.loc[name, "f_area"] * avg), "anti": [r1(ara), r1(adec)],
                    "circle": [[r1(a), r1(b)] for a, b in zip(circ.ra.deg, circ.dec.deg)]})
    plates.append(rec)
site["plates"] = plates
# plates in the catalogue that the published tests leave out (no transients on a night in the team's list)
allp = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
allp = allp[allp.index.isin(cat.plate.unique()) & ~allp.index.isin(pl.index)]
n_all = cat.groupby("plate").size()
site["plates_out"] = []
for name, p in allp.iterrows():
    ra, dec = deproject(sq[:, 0], sq[:, 1], p.ra_c, p.dec_c)
    site["plates_out"].append({"id": name, "N": int(n_all.get(name, 0)), "date": p.ut.strftime("%Y-%m-%d"),
                               "poly": [[r1(a), r1(b)] for a, b in zip(ra, dec)]})
site["avg_per_plate"] = r1(avg)

# ------------------------------------------------------------------ the plate-level measurement (shadow_fit.py): per-plate
# counts and ratios, the common fit, its calibration on toys from real plates, and the cross-checks; and, for each of the
# eight plates partly in shadow, every transient within 3 deg (xi east, eta north, deg), which side of the shadow's
# edge it lies, and the edge itself
import shadow as SH                                                          # noqa: E402
F_ = json.loads((ROOT / "data/shadow_fit.json").read_text())
site["shadow"] = {k: F_[k] for k in ("fit", "closure", "whole")}
site["shadow"]["toys"] = {k: F_["toys"][k] for k in ("null", "p_as_high", "hist_null", "interval68", "interval90")}
site["shadow"]["checks"] = F_["checks"]
plc = {p_["plate"]: p_ for p_ in F_["plates"]}
up_close = []


def project(ra, dec, ra0, dec0):
    """RA, Dec (deg) -> plate-frame offsets xi (east), eta (north) in deg, gnomonic about (ra0, dec0)."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(xi), np.degrees(eta)


for q in sorted(plc, key=lambda q: -plc[q]["shadow"]):
    p_ = pl.loc[q]; t_ = Time(pd.Timestamp(p_.ut).to_datetime64())
    c_ = cs[(cs.plate == q) & (cs.xi.abs() < SH.HALF) & (cs.eta.abs() < SH.HALF)]
    # the shadow's edge: the 8.5 deg circle about the anti-sun point at mid-exposure, on the plate's tangent plane
    ax = np.asarray(SH.antisun_axis(t_)).reshape(3)
    c0 = SkyCoord(np.degrees(np.arctan2(ax[1], ax[0])) * u.deg, np.degrees(np.arcsin(ax[2])) * u.deg)
    ring_ = c0.directional_offset_by(np.radians(np.arange(0, 360.5, 1.5)) * u.rad, SH.SHADOW_DEG * u.deg)
    ex, ey = project(ring_.ra.deg, ring_.dec.deg, p_.ra_c, p_.dec_c)
    dark = SH.antisun_geocentric(c_.ra.values, c_.dec.values, Time(np.full(len(c_), pd.Timestamp(p_.ut).to_datetime64()))) < SH.SHADOW_DEG
    up_close.append({**plc[q], "sunlit": plc[q]["n"] - plc[q]["O"], "xi": [round(float(v_), 3) for v_ in c_.xi],
                     "eta": [round(float(v_), 3) for v_ in c_.eta], "dark": "".join("1" if d_ else "0" for d_ in dark),
                     "edge": [[round(float(x_), 3), round(float(y_), 3)] for x_, y_ in zip(ex, ey)]})
site["shadow"]["plates"] = up_close
site["shadow"]["half"] = SH.HALF

# ------------------------------------------------------------------ two ideas of chance
pf = ROOT / "data/paste_sums_geo.npy"
sums = np.load(pf if pf.exists() else ROOT / "data/paste_sums_geo_all.npy")
# observed: transients in the shadowed parts, each tested against the shadow's edge (plate_shadow_geo_corr.parquet, from
# shadow_test.py); expected: the shadowed plates' worth of sky times the average plate
psh = pd.read_parquet(ROOT / "data/plate_shadow_geo_corr.parquet")
n_strict, area_exp, observed = int(((cs.xi.abs() < SH.HALF) & (cs.eta.abs() < SH.HALF)).sum()), round(float(psh.f_area.sum() * psh.N.mean()), 1), int(psh[psh.f_area > 0].n_in.sum())
pfrac = area_exp / n_strict
bins = np.arange(0, 2001, 25)
h, _ = np.histogram(sums, bins)
site["paste"] = {"bins": bins.tolist(), "counts": h.tolist(), "n": int(len(sums)), "observed": observed,
                 "area_expected": area_exp, "p": r3(np.mean(sums <= observed)),
                 "binom_sd": r1(np.sqrt(n_strict * pfrac * (1 - pfrac))), "n_transients": n_strict}

# ------------------------------------------------------------------ nuclear: nights, shots, lags, season
lee = quiet_exec(ROOT / "analysis/nuclear_lee.py", stop="# breadth statistic")
nights, values, sets, windows, members = lee["nights"], lee["values"], lee["sets"], lee["windows"], lee["members"]
ratios, local_p, shifts, ens = lee["ratios"], lee["local_p"], lee["shifts"], lee["ens"]
real = lee["real"]
site["nights"] = [{"d": d.strftime("%Y-%m-%d"), "n": int(v)} for d, v in zip(nights, values)]
shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
shots["t"] = pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601")
shots = shots[(shots.t >= "1949-11-01") & (shots.t <= "1957-05-31")]
shots["op"] = shots.source_url.str.extract(r"Operation_([A-Za-z–\-]+)")[0].str.replace("_", " ")
site["shots"] = [{"t": t.strftime("%Y-%m-%dT%H:%M"), "name": n, "site": "Nevada" if str(s).startswith("NTS") else "Pacific",
                  "op": (o if isinstance(o, str) else "")} for t, n, s, o in zip(shots.t, shots.Shot_Name, shots.Site, shots.op)]
wl = list(windows)
lagrows = []
for si, sname in enumerate(sets):
    for k in range(-7, 8):
        m = si * len(wl) + wl.index(f"lag {k:+d}")
        col = ens[:, m][np.isfinite(ens[:, m])]
        lagrows.append({"set": sname, "lag": k, "ratio": r3(real[m]), "p": r3(local_p(real[m], ens[:, m])),
                        "q05": r3(np.percentile(col, 5)), "q50": r3(np.percentile(col, 50)), "q95": r3(np.percentile(col, 95))})
site["lags"] = lagrows


def nsig_dist(shift_list):
    rr = np.concatenate([ratios(sets[s], 0) for s in sets])
    e = np.array([np.concatenate([ratios(sets[s], sh) for s in sets]) for sh in shift_list])
    p_real = np.array([local_p(rr[m], e[:, m]) for m in range(len(members))])
    rk = np.empty_like(e)
    for m in range(e.shape[1]):
        col = np.nan_to_num(e[:, m], nan=-np.inf)
        rk[:, m] = (np.argsort(np.argsort(-col)) + 1) / (len(col) + 1)
    return int(np.sum(p_real < 0.05)), np.sum(rk < 0.05, axis=1)


year_shift = np.array([365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-20, 21)])
n_any, d_any = nsig_dist(shifts)
n_yr, d_yr = nsig_dist(year_shift)
m1 = wl.index("+/-1")   # VASCO's window, all US tests
r_any = ens[:, m1][np.isfinite(ens[:, m1])]
r_yr = np.array([ratios(sets["all US"], sh)[m1] for sh in year_shift]); r_yr = r_yr[np.isfinite(r_yr)]
site["pm1"] = {"real": r3(real[m1]), "any": [r3(v) for v in r_any], "season": [r3(v) for v in r_yr],
               "p_any": r3(np.mean(r_any >= real[m1])), "p_season": r3(np.mean(r_yr >= real[m1]))}
site["lee"] = {"any": {"real": n_any, "dist": np.bincount(d_any, minlength=31)[:31].tolist(), "p": r3(np.mean(d_any >= n_any))},
               "season": {"real": n_yr, "dist": np.bincount(d_yr, minlength=31)[:31].tolist(), "p": r3(np.mean(d_yr >= n_yr))}}
pc = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
pc["b"] = np.abs(SkyCoord(pc.ra_c.values * u.deg, pc.dec_c.values * u.deg).galactic.b.deg)
nb = pc.groupby(pc.ut.dt.normalize()).b.mean().reindex(nights)
mm = pd.DataFrame({"m": nights.month, "n": values, "b": nb.values}).groupby("m").agg(n=("n", "mean"), b=("b", "mean"))
nv = shots[shots.Site.str.startswith("NTS")].t.dt.month.value_counts()
day0 = pd.Timestamp("1949-01-01")
in_win = np.isin(((nights - day0).days).values, (sets["all US"][:, None] + np.array([-1, 0, 1])).ravel())
win = pd.Series(nights[in_win].month).value_counts()
# expected transients for each plate from its distance to the Milky Way (one smooth curve fitted to all plates, rising to a
# plateau); a night's expected number is the sum over its plates. VASCO's window is then tested on count / expected.
from scipy.optimize import curve_fit
nuc = lee["g"]
pln = nuc["pl"][nuc["pl"].cor_date.isin(nights)].copy()
pln["n"] = nuc["cat"].groupby("plate").size().reindex(pln.index).fillna(0)
pln["b"] = np.abs(SkyCoord(pln.ra_c.values * u.deg, pln.dec_c.values * u.deg).galactic.b.deg)
sat = lambda x, A, b0, C: C + A * (1 - np.exp(-x / b0))
psat, _ = curve_fit(sat, pln.b, pln.n, p0=[200, 20, 30], maxfev=20000)
exp_night = pln.assign(e=sat(pln.b, *psat)).groupby("cor_date").e.sum().reindex(nights).values
adj = values / exp_night
day0 = pd.Timestamp("1949-01-01"); Nd = ((nights - day0).days).values


def win_ratio(v, shift):
    lab = np.isin(Nd, (sets["all US"][:, None] + shift + np.array([-1, 0, 1])).ravel()) & np.isfinite(v)
    rest = ~np.isin(Nd, (sets["all US"][:, None] + shift + np.array([-1, 0, 1])).ravel()) & np.isfinite(v)
    return v[lab].mean() / v[rest].mean() if 2 <= lab.sum() else np.nan


r_adj = win_ratio(adj, 0)
e_any = np.array([win_ratio(adj, sh) for sh in shifts]); e_yr = np.array([win_ratio(adj, sh) for sh in year_shift])
in_win0 = np.isin(Nd, (sets["all US"][:, None] + np.array([-1, 0, 1])).ravel())
win_dates = set(nights[in_win0])
site["skyadj"] = {"params": [r3(v_) for v_ in psat], "ratio": r3(r_adj),
                  "p_any": r3(np.nanmean(e_any[np.isfinite(e_any)] >= r_adj)), "p_season": r3(np.nanmean(e_yr[np.isfinite(e_yr)] >= r_adj)),
                  "plates": [{"b": r1(b_), "n": int(n_), "w": int(d_ in win_dates), "id": i_, "d": d_.strftime("%Y-%m-%d")}
                             for i_, b_, n_, d_ in zip(pln.index, pln.b, pln.n, pln.cor_date)]}
# VASCO's own test (window nights against the same number of observing nights picked at random), before and after
# allowing for the distance of each night's plates from the Milky Way
rngp = np.random.default_rng(42)


def random_nights(e, n=20000):
    """VASCO's test: transients on the nights near tests over their expected number, divided by the same for all other
    nights, against random sets of nights of the same size drawn from the whole survey. With e = 1 for every night this
    is VASCO's ratio of transients per night; nights with no plates count as zero on both sides."""
    lab = np.isin(Nd, (sets["all US"][:, None] + np.array([-1, 0, 1])).ravel()); k = int(lab.sum())
    V_, E_ = values.sum(), e.sum()
    stat = lambda m: (values[m].sum() / e[m].sum()) / ((V_ - values[m].sum()) / (E_ - e[m].sum()))
    obs = stat(lab); R = np.empty(n)
    for i in range(n):
        m = np.zeros(len(values), bool); m[rngp.choice(len(values), k, replace=False)] = True; R[i] = stat(m)
    from scipy.stats import norm
    pv = np.mean(R >= obs)
    return {"real": r3(obs), "nights": k, "p": r3(pv), "sigma": r1(norm.isf(pv)), "hist": np.histogram(R, np.arange(0, 3.01, 0.05))[0].tolist()}


site["perm"] = {"raw": random_nights(np.ones(len(values))), "adj": random_nights(np.nan_to_num(exp_night))}

# VASCO's unit is the night, and a night's transients count its photographs: transients per plate by month, and VASCO's
# excess split into its parts. Per night = (plates per night) x (transients per plate); per plate = (expected per plate,
# from where the plates pointed) x (transients per plate over expected). "Plates" are the catalogue's plates (those with
# at least one transient on VASCO's nights).
npl = pln.groupby("cor_date").size().reindex(nights).fillna(0).values
lab0 = np.isin(Nd, (sets["all US"][:, None] + np.array([-1, 0, 1])).ravel())
mt = pd.DataFrame({"m": nights.month, "n": values, "p": npl}).groupby("m").sum()
site["months"] = [{"m": "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()[i - 1], "n": r1(mm.n[i]), "b": r1(mm.b[i]),
                   "pp": r1(mt.n[i] / mt.p[i]), "ppn": round(float(mt.p[i] / (nights.month == i).sum()), 2),
                   "shots": int(nv.get(i, 0)), "win": int(win.get(i, 0))} for i in range(1, 13)]
en = np.nan_to_num(exp_night); w_, o_ = lab0, ~lab0
ym = nights.to_period("M"); same = o_ & np.isin(ym, ym[w_].unique())
site["nuc_steps"] = {"nights_with_plates": int((npl > 0).sum()), "nights_listed": len(npl), "max_plates": int(npl.max()),
                     "per_night": r3(values[w_].mean() / values[o_].mean()),
                     "per_plate": r3((values[w_].sum() / npl[w_].sum()) / (values[o_].sum() / npl[o_].sum())),
                     "per_expected": r3((values[w_].sum() / en[w_].sum()) / (values[o_].sum() / en[o_].sum())),
                     "pointing": r3((en[w_].sum() / npl[w_].sum()) / (en[o_].sum() / npl[o_].sum())),
                     "plates_window": round(float(npl[w_].mean()), 2), "plates_other": round(float(npl[o_].mean()), 2),
                     "plates_same_months": round(float(npl[same].mean()), 2), "nights_window": int(w_.sum()),
                     "empty_window": int((npl[w_] == 0).sum()), "empty_other": int((npl[o_] == 0).sum()), "nights_other": int(o_.sum())}

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
    site["trackdist"]["all"].append({"r": r, "real": rb(AZ["real_l"][k]).tolist(), "bmean": B.mean(0).round(2).tolist(), "bsd": B.std(0).round(2).tolist()})
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


def recovered(key, sigs, k=1, Bt=3):
    n = np.cumsum(AZ[f"null_{key}"][:, k].mean(0)); i_ = np.flatnonzero(n <= Bt)[-1]
    return r1(np.mean([np.cumsum(AZ[f"signal_{key}"][s_, k] - AZ[f"real_{key}"][k])[i_] for s_ in sigs]))


angles = sorted({float(s_[3]) for s_ in SIGS if s_[3] != "any"})
site["recovery"] = {"angles": angles, "false_alarms": 3,
                    "track": [recovered("l", [j for j, s_ in enumerate(SIGS) if s_[3] != "any" and float(s_[3]) == a_]) for a_ in angles],
                    "vasco": [recovered("v", [j for j, s_ in enumerate(SIGS) if s_[3] != "any" and float(s_[3]) == a_]) for a_ in angles],
                    "any": {"track": recovered("l", [j for j, s_ in enumerate(SIGS) if s_[1] == "1.0" and s_[3] == "any"]),
                            "vasco": recovered("v", [j for j, s_ in enumerate(SIGS) if s_[1] == "1.0" and s_[3] == "any"])}}
# the tracks that stand out (track score 1e-4 or better, each place once), for a blink comparison: the red plate the
# transients were found on, the blue plate of the same field taken the same night, and the second survey's red plate,
# all resampled onto the red plate's pixels (DSS cutouts sized to each track); and, where it shows what made a track,
# the red plate at high contrast (a scratch shows as a pale band) or a wide view (a trail runs far beyond the track)
from scipy.ndimage import map_coordinates, median_filter
from scipy.stats import chi2
EVIDENCE = {"XE017": "zoom", "XE243": "zoom", "XE409": "zoom", "XE493": "zoom", "XE492": "wide"}


def dss(survey, c, size, q):
    f = ROOT / f"data/cutouts/tracks/{q}_{survey}_{size:g}.fits"
    if not f.exists():
        f.parent.mkdir(exist_ok=True)
        f.write_bytes(requests.get(f"https://archive.stsci.edu/cgi-bin/dss_search?v={survey}&r={c.ra.deg:.5f}&d={c.dec.deg:.5f}"
                                   f"&e=J2000&h={size}&w={size}&f=fits&c=none", timeout=300).content)
    return fits.getheader(f), fits.getdata(f).astype(float)


def png(name, im, lo_, hi_):
    if "--images" in sys.argv or not (OUT / "img" / name).exists():
        fig = plt.figure(figsize=(3, 3)); ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(im, origin="lower", cmap="gray_r", vmin=lo_, vmax=hi_); ax.axis("off")
        fig.savefig(OUT / "img" / name, dpi=110, facecolor="white"); plt.close(fig)
    return f"img/{name}"


frac = lambda v_, n_: [round(float(x_), 4) for x_ in v_ / (n_ - 1)]
TT = AJ["track"]; tail = []
for a_ in [c_ for c_ in TT["candidates"] if c_["score"] <= 1e-4]:
    q = a_["plate"]; bx_, by_ = np.array(a_["x"]) / 3600, np.array(a_["y"]) / 3600
    ra_, dec_ = SH.deproject(bx_, by_, pl.loc[q].ra_c, pl.loc[q].dec_c); sk = SkyCoord(ra_ * u.deg, dec_ * u.deg)
    c0_ = SkyCoord(ra_.mean() * u.deg, dec_.mean() * u.deg); size = float(max(2, np.ceil(1.6 * a_["d_arcmin"])))
    hr_, R_ = dss("poss1_red", c0_, size, q); assert hr_["REGION"].strip() == q, (q, hr_["REGION"])
    wr_ = WCS(hr_); H_, W_ = R_.shape; yy_, xx_ = np.mgrid[0:H_, 0:W_]; world = wr_.pixel_to_world(xx_, yy_)
    blink, years = [png(f"tail_{q.lower()}_red.png", R_, *np.percentile(R_, [2, 99.7]))], [int(hr_["DATE-OBS"][:4])]
    for sv, tag in (("poss1_blue", "blue"), ("poss2ukstu_red", "red2")):
        h_, im_ = dss(sv, c0_, size, q); ox_, oy_ = WCS(h_).world_to_pixel(world)
        im_ = map_coordinates(im_, [oy_, ox_], order=1, mode="nearest")
        blink.append(png(f"tail_{q.lower()}_{tag}.png", im_, *np.percentile(im_, [2, 99.7]))); years.append(int(h_["DATE-OBS"][:4]))
    px_, py_ = wr_.world_to_pixel(sk)
    rec = {"plate": q, "r": a_["r"], "score": a_["score"], "chi2": a_["chi2ndof"], "length": a_["d_arcmin"], "size": size, "vasco": a_["vasco"],
           "xi": r3(float(bx_.mean())), "eta": r3(float(by_.mean())), "blink": blink, "years": years,
           "p_line": r3(float(chi2.sf(a_["chi2ndof"] * (a_["r"] - 2), a_["r"] - 2))), "mx": frac(px_, W_), "my": [1 - v_ for v_ in frac(py_, H_)]}
    kind = EVIDENCE.get(q)
    if kind == "zoom":                                                 # high contrast, sky removed, at most 6' across
        Rr = R_ - median_filter(R_, size=max(15, min(61, H_ // 3) | 1)); zn = 1.4826 * np.median(np.abs(Rr - np.median(Rr)))
        zs = min(size, 6.0); hz = int(round(W_ * zs / size / 2)); zc = (int(round(float(np.median(px_)))), int(round(float(np.median(py_)))))
        zx0, zx1, zy0, zy1 = max(zc[0] - hz, 0), min(zc[0] + hz, W_), max(zc[1] - hz, 0), min(zc[1] + hz, H_)
        rec["evidence"] = {"kind": kind, "size": zs, "img": png(f"tail_{q.lower()}_zoom.png", Rr[zy0:zy1, zx0:zx1], -3 * zn, 3 * zn),
                           "x": frac(px_ - zx0, zx1 - zx0), "y": [1 - v_ for v_ in frac(py_ - zy0, zy1 - zy0)]}
    elif kind == "wide":                                               # the red plate three times wider
        ws = 3 * size; h_, im_ = dss("poss1_red", c0_, ws, q); wx_, wy_ = WCS(h_).world_to_pixel(sk)
        rec["evidence"] = {"kind": kind, "size": ws, "img": png(f"tail_{q.lower()}_wide.png", im_, *np.percentile(im_, [2, 99.7])),
                           "x": frac(wx_, im_.shape[1]), "y": [1 - v_ for v_ in frac(wy_, im_.shape[0])]}
    tail.append(rec)
# how the tracks' image widths compare with the plate's stars (tail_shapes.py), and are they on a second, independent copy
# of the plate (second_copy.py)
SHP = json.loads((ROOT / "data/tail_shapes.json").read_text())
SEC = {r["plate"]: r for r in json.loads((ROOT / "data/second_copy.json").read_text())}
widths = {t["plate"]: float(np.median([r["ratio"] for r in t["transients"] if r])) for t in SHP["tracks"]}
# each track's transients against the star model of their plate (tail_shapes.py): the model (the average of the stars of
# similar peak brightness, at unit light), each transient on the same scale with its width as a multiple of the stars',
# how fast the track's images fade from the centre against the stars (1 for star-shaped images), and the track's call;
# then the same for transients picked at random. Images are scaled to the peak of the star model


def nan_none(a):
    """Rounded, with pixels left out (NaN) as nulls."""
    return [[None if not np.isfinite(v) else round(float(v), 3) for v in row] for row in a]


def star_model(rows):
    rows = [r for r in rows if r and r.get("img")]
    M = np.mean([r["model"] for r in rows], 0)
    return rows, (M / M.max()).round(3).tolist()


wt = []
for t in SHP["tracks"]:
    rows, M = star_model(t["transients"])
    wt.append({"plate": t["plate"], "call": t["call"], "p_wider": t["p_wider"], "p_narrower": t["p_narrower"], "model": M, "profile": t["profile"],
               "transients": [{"ratio": r3(r["ratio"]), "img": nan_none(np.array(r["img"], float) / np.max(r["model"]))} for r in rows]})
rows, M = star_model(SHP["control"]); CP = dict(SHP["control_profile"])
img_ = np.nanmean([np.array(r["img"], float) for r in rows], 0) / np.mean([r["model"] for r in rows], 0).max()
site["widths"] = {"tracks": wt, "level": SHP["level"], "r": SHP["r_arcsec"], "kernel": SHP["kernel_arcsec"],
                  "control": {"n": len(rows), "model": M, "img": nan_none(img_), "ratio": r3(float(np.exp(np.mean([np.log(r["ratio"]) for r in rows])))),
                              **{k: CP.pop(k) for k in ("call", "p_wider", "p_narrower")}, "profile": CP,
                              "wider": r3(float(np.mean([r["ratio"] > r["stars_p90"] for r in rows])))}}
for rec in tail:
    rec["width"] = r3(widths[rec["plate"]]) if rec["plate"] in widths else None
    rec["second_copy"] = SEC[rec["plate"]]["correlation"] if rec["plate"] in SEC else None
# the trail on XE492 (xe492_trail.py): the strip along the line through the track's transients, the traced trail and the
# fits, and where the transients sit
TR = json.loads((ROOT / "data/xe492_trail.json").read_text()); TZ = np.load(ROOT / "data/xe492_trail.npz")
if "--images" in sys.argv or not (OUT / "img/xe492_trail.png").exists():
    st_, sd_ = TZ["strip"], float(TZ["sd"])
    fig = plt.figure(figsize=(12, 1.2)); ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(np.nan_to_num(st_, nan=0), origin="lower", cmap="gray_r", vmin=-2 * sd_, vmax=6 * sd_, aspect="auto"); ax.axis("off")
    fig.savefig(OUT / "img/xe492_trail.png", dpi=200, facecolor="white"); plt.close(fig)
best492 = next(c_ for c_ in TT["candidates"] if c_["plate"] == "XE492")
P492 = np.c_[best492["x"], best492["y"]] / 3600; c492 = P492.mean(0); d492 = np.linalg.svd(P492 - c492)[2][0]; n492 = np.array([-d492[1], d492[0]])
site["trail"] = {"img": "img/xe492_trail.png", "along0": float(TZ["along0"]), "along1": float(TZ["along0"]) + TZ["strip"].shape[1] * float(TZ["step_arcmin"]),
                 "across": [float(TZ["across"][0]), float(TZ["across"][-1])], "s": TR["s"], "offset": TR["offset"], "quad": TR["quad"], "s_lo": TR["s_lo"],
                 "members_along": [r3(((P492[k] - c492) @ d492 - TR["s_lo"]) * 60) for k in range(len(P492))],
                 "members_across": [r1(((P492[k] - c492) @ n492) * 3600) for k in range(len(P492))],
                 **{k: TR[k] for k in ("start_arcmin", "end_arcmin", "length_arcmin", "ends_at_coverage", "rms_line", "rms_quad", "radius_deg", "sag_arcsec",
                                       "dots", "spacing_median", "spacing_16_84", "sections")}}
site["tail_tracks"] = {"tracks": tail, "levels": TT["levels"], "real": TT["real"], "copies": TT["copies_mean"], "share": TT["share_copies_ge_real"]}

# ------------------------------------------------------------------ dates
t = pd.read_csv(ROOT / "data/plate_date_corrections.csv", index_col=0, parse_dates=["header_date", "corrected_date"])
ev = []
for name, r in t.dropna(subset=["label"]).iterrows():
    for src, col in (("plate log", "log_offset"), ("asteroids", "asteroid_offset"), ("label", "label_offset")):
        v = r[col]
        if v is not None and not (isinstance(v, float) and np.isnan(v)) and v in (0, 1):
            ev.append({"plate": name, "e": int(r.label), "src": src, "off": int(v)})
site["date_evidence"] = ev
used = t.index.isin(pl.index) | t.index.isin(pd.read_parquet(ROOT / "data/plates.parquet").index)
an = t[t.index.isin(pd.read_parquet(ROOT / "data/plates.parquet").index)]
exc = []
for p in ("XE113", "XE389", "XE217", "XE128", "XE603", "XE666", "XE498", "XE098"):
    r = t.loc[p]
    exc.append({"plate": p, "e": int(r.label), "header": r.header_date.strftime("%Y-%m-%d"),
                "correct": r.corrected_date.strftime("%Y-%m-%d"), "days": int(r.offset_days)})
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
site["dust"] = {"n_dusty": len(dusty), "median_dusty": int(dusty.n.median()), "n_clear": len(clear), "median_clear": int(clear.n.median())}
# share of the 561 plates (VASCO's nights) whose band 2.0-2.75 deg from the centre is denser than their inner 1.5 deg
cs_ = cs[cs.plate != ""]; rs_ = np.maximum(cs_.xi.abs(), cs_.eta.abs())
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
site["crowding"] = [{"id": i, "b": r1(r.absb), "l": r1(r.l), "dec": r1(r.dec_c), "N": int(r.N)} for i, r in pa.iterrows()]

(OUT / "data/site.json").write_text(json.dumps(site, separators=(",", ":")))
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
