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
strict = set(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date.astype(str).str[:10])
cs = cat[cat.obs_date.astype(str).str[:10].isin(strict)]
pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
pl = pl[pl.index.isin(cs.plate.unique())].copy()
pl["N"] = cs.groupby("plate").size().reindex(pl.index).fillna(0).astype(int)
pt = pd.read_parquet(ROOT / "data/plate_shadow_geo_corr.parquet")
avg = pt.N.sum() / len(pt)
edge = np.linspace(-3.1, 3.1, 5)
sq = np.r_[np.c_[edge, np.full(5, -3.1)], np.c_[np.full(5, 3.1), edge][1:], np.c_[edge[::-1], np.full(5, 3.1)][1:], np.c_[np.full(5, -3.1), edge[::-1]][1:-1]]
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

# ------------------------------------------------------------------ bounded interval for phi (shadow_phi_interval.py)
site["phi"] = {g: json.loads((ROOT / f"data/phi_interval_{g}.json").read_text()) for g in ("geo", "topo")}

# ------------------------------------------------------------------ southern shadow plates, our own transient extraction
# (shadow_ours.py, shadow_phi_combined.py); plates not yet processed are drawn without a count
ss = pd.read_csv(ROOT / "data/southern_shadow_scans.csv").astype({"plate_number": str})
sl = pd.read_csv(ROOT / "data/shadow_plates_log.csv").astype({"plate_number": str})
ss = ss.merge(sl[["plate_number", "evening_date", "ra", "dec", "ut"]], on=["plate_number", "evening_date", "ra", "dec"])
res = pd.read_csv(ROOT / "data/shadow_ours.csv", index_col="plate") if (ROOT / "data/shadow_ours.csv").exists() else pd.DataFrame()
crowd = pd.read_csv(ROOT / "data/crowding.csv").drop_duplicates("plate", keep="last").set_index("plate").veto
south = []
for _, p in ss.iterrows():
    ra, dec = deproject(sq[:, 0], sq[:, 1], p.ra, p.dec)
    s_ = get_sun(Time(pd.Timestamp(p.ut)))
    c0 = SkyCoord((s_.ra.deg + 180) % 360 * u.deg, -s_.dec.deg * u.deg)
    circ = c0.directional_offset_by(np.radians(np.arange(0, 361, 6)) * u.rad, 8.5 * u.deg)
    rec = {"id": p.plate, "num": p.plate_number, "date": pd.Timestamp(p.ut).strftime("%Y-%m-%d"), "frac": r3(p.shadow_fraction),
           "crowd": r3(crowd.get(p.plate, np.nan)) if p.plate in crowd else None,
           "poly": [[r1(a_), r1(b_)] for a_, b_ in zip(ra, dec)], "circle": [[r1(a_), r1(b_)] for a_, b_ in zip(circ.ra.deg, circ.dec.deg)]}
    if p.plate in res.index:
        r_ = res.loc[p.plate]
        rec.update({"N": int(r_.n), "obs": int(r_.observed), "exp": r1(r_.expected), "frac": r3(r_.shadow_area)})
    south.append(rec)
site["south"] = {"plates": south}
if (ROOT / "data/phi_combined.json").exists():
    site["south"]["phi"] = json.loads((ROOT / "data/phi_combined.json").read_text())

# ------------------------------------------------------------------ two ideas of chance
pf = ROOT / "data/paste_sums_geo.npy"
sums = np.load(pf if pf.exists() else ROOT / "data/paste_sums_geo_all.npy")
n_strict, area_exp, observed = len(cs), 521.0, 343
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

site["months"] = [{"m": "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()[i - 1], "n": r1(mm.n[i]), "b": r1(mm.b[i]),
                   "shots": int(nv.get(i, 0)), "win": int(win.get(i, 0))} for i in range(1, 13)]

# ------------------------------------------------------------------ alignments
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
# alignments as track reconstruction (analysis/alignments_tracks.py): distribution of the straight-line fit quality
# log10(chi2/ndof) (1" hit resolution) for groups of 3-6 transients, real catalogue against 1,000 copies with every real
# line broken, plus what 100 added genuine tracks (1" scatter) contribute; whole plates and with the plate edges left out
tz = np.load(ROOT / "data/alignment_tracks.npz")
rb = lambda a: a.reshape(*a.shape[:-1], -1, 2).sum(-1)                 # 0.2-dex bins
lo = tz["edges"][:-1].reshape(-1, 2)[:, 0]
site["trackdist"] = {"lo": [r1(v_) for v_ in lo]}
for z, name in ((0, "all"), (1, "inner")):
    site["trackdist"][name] = []
    for k, r in enumerate((3, 4, 5, 6)):
        B = rb(tz["null"][:, k, z])
        site["trackdist"][name].append({"r": r, "real": rb(tz["real"][k, z]).tolist(), "bmean": B.mean(0).round(2).tolist(),
                                        "bsd": B.std(0).round(2).tolist(), "sig": rb((tz["signal"] - tz["real"])[k, z]).tolist()})

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
ctrl = pd.read_csv(ROOT / "VASCO-challenge/data/XE_final_random_sources.csv")
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
site["radial"] = {"r": (ann[:-1] + 0.125).round(3).tolist(), "transients": (dt / dt[ref].mean()).round(3).tolist(),
                  "lo": np.percentile(boot, 5, 0).round(3).tolist(), "hi": np.percentile(boot, 95, 0).round(3).tolist(),
                  "control": (dc / dc[ref].mean()).round(3).tolist()}
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
def save(fig, name):
    fig.savefig(OUT / "img" / name, dpi=150, bbox_inches="tight", pad_inches=0.02, facecolor="white"); plt.close(fig)

for plate, raw in (("XE356", "xe356_bot.raw"), ("XE617", "xe617_bot.raw")):
    a = np.fromfile(ROOT / "data/plate_edges" / raw, dtype=">i2").reshape(1600, 14000).astype(float)[1150:1500, 100:900]
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
    save(fig, f"asteroid_{nm}_wide.png")
print("images written")
