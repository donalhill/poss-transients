"""Figures for report/lost_light.tex. Run from the project root: .venv/bin/python report/make_figs.py"""
import contextlib
import io
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "report/figs"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "analysis"))

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, NULLC = "#0b0b0b", "#52514e", "#e4e3df", "#b9b8b2"
plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "lines.linewidth": 1.6, "legend.frameon": False, "legend.fontsize": 8,
    "font.family": "serif", "pdf.fonttype": 42,
})


def quiet_exec(path, stop=None, extra=None):
    g = {"__file__": str(path), "__name__": "fig"}
    g.update(extra or {})
    src = Path(path).read_text()
    if stop:
        src = src.split(stop)[0]
    with contextlib.redirect_stdout(io.StringIO()):
        exec(src, g)
    return g


def panel(ax, letter):
    ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom")


# ---------------- 1. Shadow ----------------
pt = pd.read_parquet(ROOT / "data/plate_shadow_geo_corr.parquet")
avg = pt.N.sum() / len(pt)
sh = pt[pt.f_area > 0].assign(exp=lambda d: d.f_area * avg, contrib=lambda d: d.n_in - d.f_area * avg).sort_values("contrib")
sums = np.load(ROOT / "data/paste_sums_geo_all.npy")
obs, area_exp = 343, 521
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.45})
ax = axs[0]
ax.barh(range(len(sh)), sh.contrib, color=BLUE, height=0.6)
ax.axvline(0, color=INK2, lw=0.8)
ax.set_yticks(range(len(sh)), [f"{p}  (N={n})" for p, n in zip(sh.index, sh.N)])
ax.set_xlabel("observed − expected in-shadow transients")
ax.set_title("Per-plate contribution to the deficit", loc="left")
ax.grid(axis="y", visible=False)
panel(ax, "a")
ax = axs[1]
ax.hist(sums, bins=np.arange(0, 2001, 40), color=NULLC, edgecolor="white", linewidth=0.5)
ax.axvline(obs, color=BLUE, lw=1.6)
ax.axvline(area_exp, color=INK2, lw=1, ls="--")
ax.text(1050, ax.get_ylim()[1] * 0.92, f"observed {obs} (blue)\np = {np.mean(sums <= obs):.2f}", color=INK, fontsize=8, va="top")
ax.text(1050, ax.get_ylim()[1] * 0.62, f"area expectation {area_exp}\n(dashed)", color=INK2, fontsize=8, va="top")
ax.set_xlim(0, 2000)
ax.set_xlabel("in-footprint transients when the 9 shadow\nfootprints are pasted onto random plates")
ax.set_ylabel("trials (of 20,000)")
ax.set_title("Footprint-paste null", loc="left")
panel(ax, "b")
fig.tight_layout()
fig.savefig(OUT / "shadow.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------- 2. Dates: evidence vs serial label ----------------
t = pd.read_csv(ROOT / "data/plate_date_corrections.csv", index_col=0, parse_dates=["header_date"])
lst = set(pd.to_datetime(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date))
def log_off(d):
    a, b = d in lst, (d - pd.Timedelta(days=1)) in lst
    return 1 if (a and not b) else (0 if (b and not a) else np.nan)
t["log_offset"] = t.header_date.apply(log_off)
fig, ax = plt.subplots(figsize=(7.0, 2.5))
for lo, hi in ((30, 637), (968, 1610)):
    ax.axvspan(lo, hi, color="#f1f0ec", lw=0)
rng = np.random.default_rng(0)
for col, name, c, dy in (("log_offset", "Palomar night list (VI/25)", NULLC, 0.0), ("asteroid_offset", "asteroid trails", BLUE, 0.18),
                         ("label_offset", "handwritten labels", ORANGE, -0.18)):
    s = t.dropna(subset=[col, "label"])
    s = s[s[col].isin([0, 1])]
    ax.scatter(s.label, s[col] + dy + rng.uniform(-0.04, 0.04, len(s)), s=9 if col == "log_offset" else 14,
               color=c, label=f"{name} ({len(s)})", lw=0, zorder=3 if col != "log_offset" else 2)
for p, txt in (("XE113", "XE113"), ("XE389", "XE389")):
    r = t.loc[p]
    ax.annotate(txt, (r.label, r.label_offset - 0.18), xytext=(r.label + 40, r.label_offset - 0.18 + (0.33 if r.label_offset == 0 else -0.33)),
                fontsize=7.5, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
ax.set_yticks([0, 1], ["header correct", "header one\nday early"])
ax.set_ylim(-0.6, 1.6)
ax.set_xlabel("plate serial label (E-number), in order of exposure")
ax.grid(axis="y", visible=False)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=3)
ax.text(330, 1.45, "1950–52", ha="center", color=INK2, fontsize=8)
ax.text(1290, 1.45, "1954–56", ha="center", color=INK2, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "dates.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------- 3. Nuclear: season and look-elsewhere ----------------
lee = quiet_exec(ROOT / "analysis/nuclear_lee.py", stop="# breadth statistic")
nights, values, N, sets, windows, members = lee["nights"], lee["values"], lee["N"], lee["sets"], lee["windows"], lee["members"]
ratios, local_p, shifts = lee["ratios"], lee["local_p"], lee["shifts"]

def nsig_dist(shift_list):
    real = np.concatenate([ratios(sets[s], 0) for s in sets])
    ens = np.array([np.concatenate([ratios(sets[s], sh) for s in sets]) for sh in shift_list])
    p_real = np.array([local_p(real[m], ens[:, m]) for m in range(len(members))])
    rk = np.empty_like(ens)
    for m in range(ens.shape[1]):
        col = np.nan_to_num(ens[:, m], nan=-np.inf)
        rk[:, m] = (np.argsort(np.argsort(-col)) + 1) / (len(col) + 1)
    return int(np.sum(p_real < 0.05)), np.sum(rk < 0.05, axis=1)

year_shift = np.array([365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-20, 21)])
n_any, d_any = nsig_dist(shifts)
n_yr, d_yr = nsig_dist(year_shift)

pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
pl["b"] = np.abs(SkyCoord(pl.ra_c.values * u.deg, pl.dec_c.values * u.deg).galactic.b.deg)
night_b = pl.groupby(pl.ut.dt.normalize()).b.mean().reindex(nights)
m = pd.DataFrame({"month": nights.month, "b": night_b.values, "n": values}).groupby("month").agg(n=("n", "mean"), b=("b", "mean"))
nt = lee["nt"]
nv = pd.to_datetime(nt[(nt.Country.str.upper() == "USA") & nt.Site.str.startswith("NTS")].Date)
nv = nv[nv <= "1957-04-28"]
tests_m = nv.dt.month.value_counts().reindex(range(1, 13), fill_value=0)

fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.8), gridspec_kw={"width_ratios": [1, 1, 1.35]})
months = list("JFMAMJJASOND")
ax = axs[0]
ax.bar(range(1, 13), m.n, color=BLUE, width=0.65)
for i, k in enumerate(tests_m.values, start=1):
    if k:
        ax.text(i, m.n.loc[i] + 10, str(k), ha="center", fontsize=7, color=INK2)
ax.set_xticks(range(1, 13), months)
ax.set_title("Transients per night (numbers:\nNevada shots in that month)", loc="left")
panel(ax, "a")
ax = axs[1]
ax.bar(range(1, 13), m.b, color=AQUA, width=0.65)
ax.set_xticks(range(1, 13), months)
ax.set_title("Mean galactic latitude |b|\nof the night's plates (deg)", loc="left")
panel(ax, "b")
ax = axs[2]
bins = np.arange(-0.5, 30.5, 1)
ax.hist(d_any, bins=bins, color=NULLC, edgecolor="white", lw=0.4, density=True, label="slid by any day")
ax.hist(d_yr, bins=bins, histtype="step", color=ORANGE, lw=1.4, density=True, label="slid by whole years (season kept)")
ax.axvline(n_any, color=INK2, lw=1.4)
ax.text(n_any - 0.5, 0.13, f"real, vs any-day\nslides: {n_any}\np = {np.mean(d_any >= n_any):.3f}", fontsize=7.5, color=INK2, ha="right")
ax.axvline(n_yr, color=ORANGE, lw=1.4)
ax.text(n_yr + 0.6, 0.30, f"real, season\nkept: {n_yr}\np = {np.mean(d_yr >= n_yr):.2f}", fontsize=7.5, color=ORANGE)
ax.set_xlabel("analyses (of 56) with local p < 0.05")
ax.set_title("Look-elsewhere: how many of 56\nanalyses light up by chance", loc="left")
ax.set_xlim(-0.5, 25.5); ax.set_ylim(0, 0.65)
ax.legend(loc="upper right", fontsize=7, bbox_to_anchor=(1.0, 1.0))
panel(ax, "c")
fig.tight_layout()
fig.savefig(OUT / "nuclear.pdf", bbox_inches="tight")
plt.close(fig)
print(f"nuclear: any-day real {n_any} p {np.mean(d_any >= n_any):.3f}; season-kept real {n_yr} p {np.mean(d_yr >= n_yr):.3f}")

# ---------------- 4. Alignments ----------------
h = np.load(ROOT / "data/alignment_triples.npz")["h"]
edges = np.arange(0, 41, 1.0); cen = edges[:-1] + 0.5
n, _ = np.histogram(h, edges)
fitm = (cen >= 4) & (cen <= 40)
c = np.polyfit(cen[fitm], n[fitm], 1)
inj = quiet_exec(ROOT / "analysis/alignments_injection.py", stop="base_ex, base_sd")
al, rng2 = inj["al"], np.random.default_rng(7)
plates_ = al.cat.plate.unique(); rows = []
for _ in range(100):
    g = al.cat[al.cat.plate == rng2.choice(plates_)]
    r = rng2.integers(3, 6); L = rng2.uniform(120, 600); th = rng2.uniform(0, np.pi)
    x0, y0 = g.x.iloc[int(rng2.integers(len(g)))], g.y.iloc[int(rng2.integers(len(g)))]
    s = np.sort(rng2.uniform(0, L, r)); perp = rng2.normal(0, 1.0, r)
    rows.append(pd.DataFrame({"plate": g.plate.iloc[0], "x": x0 + s * np.cos(th) - perp * np.sin(th), "y": y0 + s * np.sin(th) + perp * np.cos(th)}))
injcat = pd.concat([al.cat[["plate", "x", "y"]]] + rows)
h_inj = np.concatenate([al.triples(gg.x.values, gg.y.values)[0] for _, gg in injcat.groupby("plate")])
n_inj, _ = np.histogram(h_inj, edges)

angles = np.arange(0, 91, 1)
rr = np.random.default_rng(3); pas = []
for th in angles:
    k = 0
    for _ in range(600):
        s = np.sort(rr.uniform(0, 300, 4))
        x = s * np.cos(np.radians(th)) + rr.normal(0, 1, 4); y = s * np.sin(np.radians(th)) + rr.normal(0, 1, 4)
        k += abs(np.corrcoef(x, y)[0, 1]) > 0.99
    pas.append(k / 600)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={"wspace": 0.35})
ax = axs[0]
ax.step(edges[:-1], n_inj / 1000, where="post", color=ORANGE, lw=1.2, label="with 100 injected aligned groups")
ax.step(edges[:-1], n / 1000, where="post", color=BLUE, lw=1.6, label="data")
ax.plot(cen, np.polyval(c, cen) / 1000, color=INK2, ls="--", lw=1, label="linear baseline (fit 4–40″)")
ax.set_xlabel("offset of third transient from the line (″)")
ax.set_ylabel("configurations per 1″ (thousands)")
ax.set_title("Collinearity: no spike near zero offset", loc="left")
ax.legend(loc="upper left", fontsize=7)
ax.set_xlim(0, 40)
panel(ax, "a")
ax = axs[1]
ax.plot(angles, pas, color=BLUE)
ax.set_xlabel("line orientation (deg from east–west)")
ax.set_ylabel("pass rate of |r| > 0.99 cut")
ax.set_title("VASCO's alignment cut vs line direction", loc="left")
ax.set_xlim(0, 90); ax.set_ylim(0, 1.05)
ax.annotate("east–west: where GEO\nglints would line up", (0, pas[0]), xytext=(12, 0.25), fontsize=8, color=INK2,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
panel(ax, "b")
fig.tight_layout()
fig.savefig(OUT / "alignments.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------- 5. Crowding ----------------
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
p2 = pd.read_parquet(ROOT / "data/plates.parquet")
p2["N"] = cat.groupby("plate").size().reindex(p2.index).fillna(0)
gal = SkyCoord(p2.ra_c.values * u.deg, p2.dec_c.values * u.deg).galactic
p2["absb"], p2["l"] = np.abs(gal.b.deg), gal.l.deg
p2 = p2[p2.N > 0]
inner = (p2.l > 20) & (p2.l < 80) & (p2.dec_c > 0)
anti = (p2.l > 140) & (p2.l < 220) & (p2.dec_c > 0)
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axs[0]
other = ~(inner | anti)
ax.scatter(p2.absb[other], p2.N[other], s=8, color=NULLC, lw=0, label="other fields")
ax.scatter(p2.absb[anti], p2.N[anti], s=11, color=ORANGE, lw=0, label="anticentre (l 140–220°), sparse")
ax.scatter(p2.absb[inner], p2.N[inner], s=11, color=BLUE, lw=0, label="inner Galaxy (l 20–80°), dense")
ax.set_yscale("log")
ax.set_xlabel("|galactic latitude| of plate centre (deg)")
ax.set_ylabel("transients per plate")
ax.set_title("Fewer survivors where stars are crowded", loc="left")
ax.legend(loc="lower right", fontsize=7)
panel(ax, "a")
ax = axs[1]
lowb = p2.absb < 20
vals = [p2.N[inner & lowb].median(), p2.N[anti & lowb].median()]
ax.bar([0, 1], vals, color=[BLUE, ORANGE], width=0.55)
for i, v in enumerate(vals):
    ax.text(i, v + 2, f"{v:.1f}", ha="center", fontsize=8)
ax.set_xticks([0, 1], [f"inner Galaxy\n({(inner & lowb).sum()} plates)", f"anticentre\n({(anti & lowb).sum()} plates)"])
ax.set_ylabel("median transients per plate")
ax.set_title("Northern fields, |b| < 20°", loc="left")
ax.grid(axis="x", visible=False)
panel(ax, "b")
fig.tight_layout()
fig.savefig(OUT / "crowding.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------- 6. Raster evidence (asteroids, label) ----------------
import shutil
shutil.copy(ROOT / "figs/confirm_day_offsets.png", OUT / "asteroids.png")
a = np.fromfile(ROOT / "data/plate_edges/xe356_bot.raw", dtype=">i2").reshape(1600, 14000).astype(float)[1150:1500, 100:900]
fig, ax = plt.subplots(figsize=(5.0, 2.2))
lo, hi = np.percentile(a, [1, 99.7])
ax.imshow(a, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)
ax.axis("off")
fig.savefig(OUT / "label_xe356.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("figures written to", OUT)
