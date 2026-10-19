"""Build docs/data/site.json and docs/img/* for the walkthrough page from the analysis outputs.
Run from the repo root after the analysis scripts: python analysis/build_site.py"""
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
from astropy.visualization import ZScaleInterval
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
sums = np.load(ROOT / "data/shadow_paste_sums.npy")
bins = np.arange(0, 2001, 25)
site["paste"] = {"bins": bins.tolist(), "counts": np.histogram(sums, bins)[0].tolist(), "n": int(len(sums)),
                 "observed": SP["ours_geo"]["observed"], "area_expected": SP["ours_geo"]["expected"],
                 "p": SP["ours_geo"]["p_random"], "binom_sd": SP["ours_geo"]["binom_sd"]}

# ------------------------------------------------------------------ nuclear tests (nuclear.py)
NU = json.loads((ROOT / "data/nuclear.json").read_text())
NF = json.loads((ROOT / "data/nuclear_fit.json").read_text())
NC = json.loads((ROOT / "data/nuclear_checks.json").read_text())
site["perm"] = NU["perm"]
site["perm"]["corrected"] = NU["corrected"]
site["nuc_scan"] = {"grid": NU["scan"]["grid"], "rank": NU["scan"]["vasco_rank"]}
site["nuc_checks"] = {"base": NU["corrected"]["real"], "rows": [x for x in NU["checks"] if x["label"] != "All of them together"]}
NB_ = json.loads((ROOT / "data/nuclear_basic.json").read_text())
# the nuclear measurement from the ground up (nuclear_basic.py), with the full model fit (nuclear_fit.py) as a cross-check
lro = {x["run"]: x["r"] for x in NB_["leave_run_out"]}
abcd = {x["split"]: x for x in NB_["abcd"]}
# the two plates near tests that carry most of the excess, each beside a typical plate from its own run
site["nuc_plates"] = [{**q, "near": near, "xi": cat.xi[cat.plate == q["plate"]].round(3).tolist(), "eta": cat.eta[cat.plate == q["plate"]].round(3).tolist()}
                      for x in sorted(NB_["plates"][:2], key=lambda x: x["date"]) for q, near in ((x, True), (x["mate"], False))]
def label(x):
    """A chart label for one robustness check (nuclear_basic.py)."""
    c = {"2' cells holding a transient": "Counting 2′ cells holding a transient", "central 5 x 5 deg only": "Central 5 × 5° of each plate only"}.get(x["check"], x["check"])
    if "plates_dropped" in x:
        c = c.replace("without ", "without Operation ")
    return c[0].upper() + c[1:]


site["nuc"] = {"r": NB_["within_runs_expected"], "runs": NB_["runs"], "grid": NF["scan"]["grid"],
               "null": {"r": NF["fit"]["r"], "p": NF["null_runs"]["p_as_high"], "hist": NF["null_runs"]["hist"]},
               "measures": [
    {"label": "Analysis fixed in advance", "r": NF["fit"]["r"], "lo": NF["fit"]["count_lo"], "hi": NF["fit"]["count_hi"], "group": "catalogue"},
    {"label": "Step-by-step ratio", "r": NB_["within_runs_expected"], "group": "catalogue"},
    *[{"label": f"Step by step, {label(x)[0].lower()}{label(x)[1:]}", "r": x["r"], "group": "catalogue"} for x in NB_["robust"] if x["check"].startswith("without XE")],
    *[{"label": f"Runs split after {x['gap']} empty nights", "r": x["r"], "group": "choices"} for x in NF["robust"]["run_gap"] if x["gap"] != 4],
    {"label": "Counting 2′ cells holding a transient", "r": NF["robust"]["outcomes"][1]["r"], "group": "choices"},
    {"label": "Central 5 × 5° of each plate only", "r": NF["robust"]["outcomes"][2]["r"], "group": "choices"},
    *[{"label": f"Without Operation {x['series']}", "r": x["r"], "group": "choices"} for x in NF["robust"]["leave_series_out"] if x["plates_dropped"] >= 2],
    {"label": "In the central part of the plate", "r": abcd["edge against central"]["odds_ratio"], "group": "which"},
    {"label": "No neighbour within 5′", "r": abcd["crowded against isolated"]["odds_ratio"], "group": "which"},
    {"label": "Shaped like the stars", "r": abcd["shaped like stars or not"]["odds_ratio"], "group": "which"},
    {"label": "Red-plate specks", "r": NC["blue"]["red"]["r"], "lo": NC["blue"]["red"]["count_lo"], "hi": NC["blue"]["red"]["count_hi"], "group": "scans"},
    {"label": "Blue-plate specks, same nights", "r": NC["blue"]["blue"]["r"], "lo": NC["blue"]["blue"]["count_lo"], "hi": NC["blue"]["blue"]["count_hi"], "group": "scans"},
    *[{"label": f"Tracks scoring {t['score']} or better", "r": t["r"], "lo": t["count_lo"], "hi": t["count_hi"], "group": "tracks"} for t in NC["tracks"]]]}

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
site["crowding"] = [{"id": i, "b": r1(r.absb), "l": r1(r.l), "dec": r1(r.dec_c), "N": int(r.N)} for i, r in pa.iterrows()]

# ------------------------------------------------------------------ numbers quoted in the text ({n:key} in site/copy.md)
WORDS = "zero one two three four five six seven eight nine ten eleven twelve".split()
word = lambda n: WORDS[int(n)] if 0 <= int(n) < len(WORDS) else f"{int(n):,}"
odds = lambda p_: word(round(1 / p_)) if round(1 / p_) <= 10 else f"{round(1 / p_):,}"    # "one in {odds}"
num = lambda v, d=0: f"{v:,.{d}f}"
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
    "cat.transients": num(len(cat)), "cat.plates": num(len(pl)), "cat.vasco_plates": num(len(vasco_plates)),
    "cat.left_out": num(len(pl) - len(vasco_plates)), "cat.min": num(N_.min()), "cat.max": num(N_.max()),
    "cat.avg": num(avg_all), "cat.noise": num(sd_all), "cat.band_lo": num(lo_b), "cat.band_hi": num(hi_b),
    "cat.in_band": num(((N_ >= lo_b) & (N_ <= hi_b)).sum()), "cat.under30": num((N_ < 30).sum()), "cat.over500": num((N_ > 500).sum()),
    "cat.over1000": word((N_ > 1000).sum()), "cat.edge_share": num(100 * share_edge),
    "crowd.inner": num(pa[pa.absb < 5].N.median()), "crowd.outer": num(round(pa[pa.absb > 40].N.median(), -1)),
    "crowd.ratio": word(round(pa[pa.absb > 40].N.median() / pa[pa.absb < 5].N.median())),
    "dust.n_dusty": num(site["dust"]["n_dusty"]), "dust.median_dusty": num(site["dust"]["median_dusty"]),
    "dust.n_clear": num(site["dust"]["n_clear"]), "dust.median_clear": num(site["dust"]["median_clear"]),
    # Earth's shadow, pooled (shadow_pooled.py)
    "sh.plates": word(G["plates"]), "sh.plates_n": num(G["plates"]), "sh.partial": word(len(F_["plates"])),
    "sh.never": num(len(pl) - G["plates"]), "sh.observed": num(G["observed"]), "sh.worth": num(G["plates_worth"], 2),
    "sh.avg": num(G["avg_per_plate"], 1), "sh.expected": num(G["expected"]), "sh.binom_sd": num(G["binom_sd"]),
    "sh.binom_z": num(abs(G["binom_z"]), 1), "sh.p_random": num(100 * G["p_random"]), "sh.sig_random": sig(G["p_random"]),
    "sh.odds_random": odds(G["p_random"]), "sh.odds_near": odds(G["p_near"]), "sh.sig_near": sig(G["p_near"]),
    "sh.centre_obs": num(G["centre"]["observed"]), "sh.centre_exp": num(G["centre"]["expected"]), "sh.odds_centre": odds(G["centre"]["p_random"]),
    "sh.topo_plates": word(T_["plates"]), "sh.topo_obs": num(T_["observed"]), "sh.topo_exp": num(T_["expected"]),
    "sh.odds_topo": odds(T_["p_random"]), "sh.sig_topo": sig(T_["p_random"]),
    "sh.vasco_obs": num(V_["observed"]), "sh.vasco_exp": num(V_["expected"]), "sh.vasco_plates": word(V_["plates"]),
    "sh.taurus": word(sum(1 for q in pt.index[pt.f_area > 0] if 40 < pl.loc[q].ra_c < 75)),
    # Earth's shadow, plate by plate (shadow_fit.py)
    "fit.r": num(F_["fit"]["r"], 2), "fit.count_lo": num(F_["fit"]["count_lo"], 2), "fit.count_hi": num(F_["fit"]["count_hi"], 2),
    "fit.null_median": num(F_["toys"]["null"]["median"], 2), "fit.null_lo": num(F_["toys"]["null"]["q16"], 2),
    "fit.null_hi": num(F_["toys"]["null"]["q84"], 2), "fit.ci90_lo": num(F_["toys"]["interval90"][0], 2),
    "fit.ci90_hi": num(F_["toys"]["interval90"][1], 2), "fit.inj_r": num(F_["injected"]["r_true"], 2),
    "fit.inj_odds": odds(F_["injected"]["p_as_high"]),
    "fit.max_plate": max(Fp, key=lambda q: Fp[q]["r"] or 0), "fit.max_r": num(max(p_["r"] or 0 for p_ in F_["plates"])),
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
numbers.update({"sh.taurus_cap": numbers["sh.taurus"].capitalize(), "sh.south": num((~south.in_catalogue).sum()),
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
    "dates.by_label": num(src["label"]), "dates.by_runs": num(src["block rule"]), "dates.unresolved": word(src.get("unresolved", 0)),
    "dates.exceptions": word(len(DS["exceptions"])).capitalize(), "dates.gross": word(sum(abs(e_["days"]) > 1 for e_ in DS["exceptions"])).capitalize(),
    "dates.vasco_nights": num(len(sample.vasco_nights())), "dates.dropped": num((left_ & (off_ == 0)).sum()),
    "dates.shifted": num((~left_ & (off_ == 0)).sum())})
# the alignments (alignments.py, vasco_candidates.py, tail_shapes.py, xe492_trail.py)
VC, TK, CM, RC = json.loads((ROOT / "data/vasco_candidates.json").read_text()), AJ["track"], site["compare"], site["recovery"]
i5, i4 = TK["levels"].index(-5.0), TK["levels"].index(-4.0)
calls = [w_["call"] for w_ in SHP["tracks"]]
numbers.update({
    "align.vasco_transients": num(VC["n_transients"]), "align.vasco_in_cat": word(VC["n_in_catalogue"]),
    "align.copies": num(AJ["ncopies"]), "align.shift": num(AJ["shift_arcsec"]),
    "align.n5": word(TK["real"][i5]).capitalize(), "align.exp5": num(TK["copies_mean"][i5], 1), "align.n4": word(TK["real"][i4]),
    "align.n4_cap": word(TK["real"][i4]).capitalize(), "align.exp4": word(round(TK["copies_mean"][i4])),
    "align.p_global": num(100 * site["trackdist"]["tail"]["p_global"]),
    "align.rho": num(CM["rho"], 2), "align.top50": num(CM["top50"]), "align.vasco_strong": num(CM["vasco_strong"][0]),
    "align.vasco_bent": num(CM["vasco_strong"][1]), "align.both": num(CM["counts"][2]),
    "align.false_alarms": word(RC["false_alarms"]), "align.vasco_ew": num(RC["vasco"][RC["angles"].index(0.0)]),
    "align.track_ew": num(RC["track"][RC["angles"].index(0.0)]),
    "widths.wider": word(calls.count("wider")).capitalize(), "widths.like": word(calls.count("like the stars")),
    "widths.control_n": num(SHP and site["widths"]["control"]["n"]), "widths.control_pct": num(100 * (site["widths"]["control"]["ratio"] - 1)),
    "widths.defocus_km": num(5 * round(SHP["defocus_km"] / 5)),
    "trail.length": num(np.floor(TR["length_arcmin"])), "trail.sag": num(round(TR["sag_arcsec"], -1)), "trail.radius": num(TR["radius_deg"]),
    "trail.turn_km": num(SHP["defocus_km"] * np.tan(np.radians(TR["radius_deg"]))),
    "trail.spacing_lo": num(TR["spacing_16_84"][0]), "trail.spacing_hi": num(TR["spacing_16_84"][1])})
# the nuclear tests (nuclear.py, nuclear_basic.py, nuclear_fit.py, nuclear_checks.py)
pc = lambda p_: num(100 * p_)
days_from = lambda w: (f"{abs(w[0])} day{'s' if abs(w[0]) != 1 else ''} {'after' if w[0] > 0 else 'before'} a test" if w[0] == w[1] and w[0] != 0
                       else f"{w[0]} to {w[1]} days from a test".replace("-", "−"))
bysize = {x["days"]: x for x in NB_["scan"]["by_size"]}
phys = {x["window"]: x for x in NB_["physical"]}
pred = {x["feature"]: x for x in NC["prediction"]}
rob = {x["check"]: x for x in NB_["robust"]}
fs = {x["days"]: x for x in NF["scan"]["by_size"]}; nfp = {x["window"]: x for x in NF["physical"]}
tp = NB_["plates"][:2]; wo = [rob[f"without {x['plate']}"]["r"] for x in tp]; both = rob[f"without {tp[0]['plate']} and {tp[1]['plate']}"]
choice = [x["r"] for x in NB_["robust"] if not x["check"].startswith("without XE") and x.get("plates_dropped", 2) >= 2]
RS = NB_["runs_summary"]; RN = NC["rich_nights"]
longdate = lambda d: pd.Timestamp(d).strftime("%-d %B %Y")
numbers.update({
    "nuc.nights": num(ST["nights_window"]), "nuc.ratio": num(NU["perm"]["raw"]["real"], 2), "nuc.odds": odds(NU["perm"]["raw"]["p"]),
    "nuc.sig": num(NU["perm"]["raw"]["sigma"], 1), "nuc.plates_window": num(ST["plates_window"], 1), "nuc.plates_other": num(ST["plates_other"], 1),
    "nuc.per_plate": num(ST["per_plate"], 2),
    **{f"nc2.{k}_r": num(x["r"], 2) for k, x in zip(("night", "photo", "point"), NU["corrections"])},
    **{f"nc2.{k}_pv": pc(x["vasco_null"]["p"]) for k, x in zip(("night", "photo", "point"), NU["corrections"])},
    **{f"nc2.{k}_pc": pc(x["calendar_null"]["p"]) for k, x in zip(("night", "photo", "point"), NU["corrections"])},
    "nc2.point_odds": odds(NU["corrections"][2]["calendar_null"]["p"]),
    "nc2.photo_cal_median": num(NU["corrections"][1]["calendar_null"]["median"], 2), "nc2.point_cal_median": num(NU["corrections"][2]["calendar_null"]["median"], 2),
    "nc2.empty_near": word(ST["empty_window"]), "nc2.empty_other": num(ST["empty_other"]), "nc2.nights_other": num(ST["nights_other"]),
    "nc2.med_low": num(NU["pointing"]["median_low"]), "nc2.med_high": num(NU["pointing"]["median_high"]),
    "nc2.exp_pct": num(round(100 * (NU["pointing"]["expected_near_over_other"] - 1))), "nc2.spring": pc(NU["pointing"]["near_spring"]),
    "nc2.slides": num(NU["pointing"]["slides"]), "nc2.sig": num(NU["corrected"]["sigma"], 1), "nc2.odds": odds(NU["corrected"]["p"]),
    "nc3.lo": num(min(x["r"] for x in NU["checks"]), 2), "nc3.hi": num(max(x["r"] for x in NU["checks"]), 2),
    "nc3.p_lo": pc(min(x["p"] for x in NU["checks"])), "nc3.p_hi": pc(max(x["p"] for x in NU["checks"])),
    "nc4.plate": NU["scan"]["best_plate"]["plate"], "nc4.plate_date": longdate(NU["scan"]["best_plate"]["date"]),
    "nc4.plate_n": num(NU["scan"]["best_plate"]["n"]), "nc4.plate_e": num(round(NU["scan"]["best_plate"]["expected"], -1)),
    "nc4.plate_share": pc(NU["scan"]["best_plate"]["share_of_window"]),
    "nc4.rank": num(NU["scan"]["vasco_rank"]), "nc4.windows": num(NU["scan"]["windows"]),
    "nc4.rank_share": pc(NU["scan"]["share_slides_rank_as_high"]),
    "nc4.best": days_from(NU["scan"]["by_size"][2]["best"]), "nc4.best_r": num(NU["scan"]["by_size"][2]["best_r"], 1),
    "nc4.best_p": num(100 * NU["scan"]["by_size"][2]["best_p_local"], 1),
    **{f"nc4.g{x['days']}": pc(x["p_global"]) for x in NU["scan"]["by_size"]},
    "nc3.all": num(next(x["r"] for x in NU["checks"] if x["label"] == "All of them together"), 2),
    "nb.plates": num(NB_["plates_near"]), "nb.nights": num(NB_["nights_near"]), "nb.runs": word(len(NB_["runs"])),
    "nb.raw": num(NB_["raw"], 2), "nb.within": num(NB_["within_runs"], 2), "nb.r": num(NB_["within_runs_expected"], 2),
    "nb.p": pc(NB_["null_runs"]["p"]), "nb.odds": odds(NB_["null_runs"]["p"]), "nb.sig": num(NB_["null_runs"]["sigma"], 1),
    "nb.null_median": num(NB_["null_runs"]["median"], 2), "nb.null_lo": num(NB_["null_runs"]["q16"], 2), "nb.null_hi": num(NB_["null_runs"]["q84"], 2),
    "nb.p_cal": num(100 * NB_["null_calendar"]["p"], 1), "nb.upper": num(NB_["upper95"], 1), "nb.vasco_low": pc(NB_["vasco_p_as_low"]),
    "nb.runs_sep": word(sum(1 for x in NB_["runs"] if x["lo"] > 1)), "nb.runs_close": word(RS["within_10pct"]),
    "nb.run_median": num(RS["median"], 2), "nb.without_nov51": num(lro["1951-11-03"], 2), "nb.without_mar55": num(lro["1955-03-19"], 2),
    "nb.nov51": num(next(x["contrast"] for x in NB_["runs"] if x["from"] == "1951-11-03"), 1),
    "nb.mar55": num(next(x["contrast"] for x in NB_["runs"] if x["from"] == "1955-03-19"), 1),
    **{f"nb.p{k + 1}{f}": v for k, x in enumerate(tp) for f, v in (("", x["plate"]), ("_date", longdate(x["date"])), ("_n", num(x["n"])),
       ("_e", num(round(x["expected"]))), ("_ratio", num(x["ratio"], 1)), ("_mate", x["mate"]["plate"]), ("_mate_n", num(x["mate"]["n"])),
       ("_without", num(wo[k], 2)))},
    **{f"nb.p{k + 1}{m_}{f}": v for k, x in enumerate(tp) for m_, q_ in (("", x), ("_mate", x["mate"])) for f, v in
       (("_edge", pc(q_["edge"])), ("_nn", num(q_["nn_arcmin"], 1)), ("_fwhm", num(q_["star_fwhm"], 1)), ("_fwhm_pct", num(q_["star_fwhm_pct"])),
        ("_narrower", pc(q_["narrower"])))},
    "nb.crop": max(tp, key=lambda x: x["ratio"])["plate"], "nb.narrow_corr": num(NB_["shapes"]["narrower_star_width_corr"], 2),
    "nb.without_both": num(both["r"], 2), "nb.without_both_p": pc(both["p"]),
    "nb.choice_lo": num(min(choice), 2), "nb.choice_hi": num(max(choice), 2),
    "nf.r": num(NF["fit"]["r"], 2), "nf.p": num(100 * NF["null_runs"]["p_as_high"], 1), "nf.odds": odds(NF["null_runs"]["p_as_high"]),
    "nf.sig": num(NF["fit"]["sigma"], 1), "nf.p_cal": num(100 * NF["null_calendar"]["p_as_high"], 1),
    "nf.i68_lo": num(NF["interval68"][0], 2), "nf.i68_hi": num(NF["interval68"][1], 2), "nf.i90_lo": num(NF["interval90"][0], 1), "nf.i90_hi": num(NF["interval90"][1], 1),
    "nf.power": pc(NF["power_vasco"]["p_as_low"]), "nf.vasco_r": num(NF["power_vasco"]["r_true"], 2),
    "nf.rank_others": word(NF["scan"]["primary_rank"] - 1), "nf.windows": num(NF["scan"]["windows"]), "nf.rank_p": num(100 * NF["scan"]["share_trials_rank_as_high"], 1),
    **{f"nf.g{x['days']}": num(100 * x["p_global"], 0 if x["p_global"] >= 0.1 else 1) for x in NF["scan"]["by_size"]},
    "nf.g14_best": f"{fs[14]['best'][0]} to {fs[14]['best'][1]}".replace("-", "−"), "nf.g14_best_r": num(fs[14]["best_r"], 1),
    "nf.after_r": num(nfp["0 to 24 hours after a shot"]["r"], 2), "nf.after_p": num(100 * nfp["0 to 24 hours after a shot"]["p_runs"], 1),
    "nf.before_r": num(nfp["24 hours before a shot"]["r"], 2), "nf.before_p": num(100 * nfp["24 hours before a shot"]["p_runs"], 1),
    "nf.after2_r": num(nfp["24 to 48 hours after"]["r"], 2),
    "nf.gap_lo": num(min(x["r"] for x in NF["robust"]["run_gap"]), 2), "nf.gap_hi": num(max(x["r"] for x in NF["robust"]["run_gap"]), 2),
    "nf.gap_pmax": num(100 * max(x["p_runs"] for x in NF["robust"]["run_gap"]), 1),
    "nf.cells_r": num(NF["robust"]["outcomes"][1]["r"], 2), "nf.inner_r": num(NF["robust"]["outcomes"][2]["r"], 2),
    "nf.series_lo": num(min(x["r"] for x in NF["robust"]["leave_series_out"]), 2), "nf.series_hi": num(max(x["r"] for x in NF["robust"]["leave_series_out"]), 2),
    "nf.run_lo": num(min(x["r"] for x in NF["robust"]["leave_run_out"]), 2), "nf.run_hi": num(max(x["r"] for x in NF["robust"]["leave_run_out"]), 2),
    "nb.drop2_r": num(NB_["drop_two"]["r"], 2), "nb.drop2_median": num(NB_["drop_two"]["trials_median"], 2), "nb.drop2_p": num(100 * NB_["drop_two"]["p"], 1),
    "nb.plates_inf": num(sum(x["near"] for x in NB_["runs"])), "nb.runs_above": word(RS["above"]),
    "nc.tracks10_p": pc(1 - NC["tracks"][0]["p_runs"]),
    "nb.fit_r": num(NF["fit"]["r"], 2), "nb.fit_p": num(100 * NF["null_runs"]["p_as_high"], 1),
    "nb.rank": num(NB_["scan"]["primary_rank"]), "nb.windows": num(NB_["scan"]["windows"]),
    "nb.g3": pc(bysize[3]["p_global"]), "nb.g7": pc(bysize[7]["p_global"]), "nb.g14": pc(bysize[14]["p_global"]),
    "nb.g14_best": f"{bysize[14]['best'][0]} to {bysize[14]['best'][1]}".replace("-", "−"), "nb.g14_best_r": num(bysize[14]["best_r"], 1),
    "nb.after_r": num(phys["0 to 24 hours after a shot"]["r"], 2), "nb.before_r": num(phys["24 hours before a shot"]["r"], 2),
    "nb.after2_r": num(phys["24 to 48 hours after"]["r"], 2),
    "nb.central": num(abcd["edge against central"]["odds_ratio"], 2), "nb.isolated": num(abcd["crowded against isolated"]["odds_ratio"], 2),
    "nb.isolated_p": num(100 * abcd["crowded against isolated"]["p_less_signal_like"], 1),
    "nb.shapes": num(NB_["shapes"]["measured"]), "nb.shapes_plates": num(NB_["shapes"]["plates"]),
    "nb.starlike": pc(NB_["shapes"]["star_like"] / NB_["shapes"]["measured"]), "nb.wider": pc(NB_["shapes"]["wider"] / NB_["shapes"]["measured"]),
    "nb.shape_r": num(abcd["shaped like stars or not"]["odds_ratio"], 2), "nb.shape_q05": num(abcd["shaped like stars or not"]["trials_q05"], 2),
    "nb.shape_q95": num(abcd["shaped like stars or not"]["trials_q95"], 2), "nb.shape_if_star": num(NB_["shapes_if_excess"]["all_star_like"], 1),
    "nb.shape_if_none": num(NB_["shapes_if_excess"]["none_star_like"], 1),
    "nb.shape_p_more": pc(abcd["shaped like stars or not"]["p_more_signal_like"]), "nb.shape_p_less": pc(abcd["shaped like stars or not"]["p_less_signal_like"]),
    "nr.top1": longdate(RN["top_near"][0]["date"]), "nr.top1_r": num(RN["top_near"][0]["ratio"], 1),
    "nr.top2": longdate(RN["top_near"][1]["date"]), "nr.top2_r": num(RN["top_near"][1]["ratio"], 1),
    "nr.threshold": num(RN["threshold"], 1), "nr.rich": word(RN["rich"]), "nr.rich_near": word(RN["rich_near"]), "nr.p_chance": pc(RN["p_chance"]),
    "nc.plates": num(NC["blue"]["plates"]), "nc.red_r": num(NC["blue"]["red"]["r"], 2), "nc.blue_r": num(NC["blue"]["blue"]["r"], 2),
    "nc.pred_p": pc(pred["VASCO's window"]["p"]), "nc.pred_after_p": num(100 * pred["the day after"]["p"], 1),
    "nc.pred_before_p": pc(pred["the day before"]["p"]), "nc.pred_day_p": pc(pred["the test's date"]["p"]), "nc.pred_nights": num(pred["VASCO's window"]["nights_held_out"]),
    "nc.tracks10_r": num(NC["tracks"][0]["r"], 2), "nc.tracks100_r": num(NC["tracks"][1]["r"], 2),
    "nc.tracks10_near": num(NC["tracks"][0]["tracks_near"]), "nc.tracks100_near": word(NC["tracks"][1]["tracks_near"])})
site["numbers"] = numbers

(OUT / "data/site.json").write_text(json.dumps(site, separators=(",", ":"), allow_nan=False))
print("site.json", round((OUT / "data/site.json").stat().st_size / 1024), "kB")

# ------------------------------------------------------------------ images (DSS plate material, figure-sized crops)
IMGS = ("label_xe356.png", "label_xe617.png", "asteroid_isabella_wide.png", "asteroid_erigone_wide.png", "nuc_densest.png")
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
# the densest 20' patch of the richer (against expected) of the two plates that carry most of the excess near tests,
# every catalogue transient circled
q = max(NB_["plates"][:2], key=lambda x: x["ratio"])["plate"]
g = cat[cat.plate == q]; H2, xe, ye = np.histogram2d(g.xi, g.eta, bins=18, range=[[-3, 3], [-3, 3]]); i, j = np.unravel_index(H2.argmax(), H2.shape)
k = (g.xi >= xe[i]) & (g.xi < xe[i + 1]) & (g.eta >= ye[j]) & (g.eta < ye[j + 1])
f = CUT / f"{q}_densest_20.fits"
if not f.exists():
    f.write_bytes(requests.get(f"https://archive.stsci.edu/cgi-bin/dss_search?v=poss1_red&r={g[k].ra.median():.5f}&d={g[k].dec.median():.5f}"
                               f"&e=J2000&h=20&w=20&f=fits&c=none", timeout=180).content)
hd = fits.getheader(f); assert hd["REGION"].strip() == q, (q, hd["REGION"])
img, w = fits.getdata(f).astype(float), WCS(hd)
fig = plt.figure(figsize=(9, 9)); ax = fig.add_axes([0, 0, 1, 1])
lo, hi = ZScaleInterval().get_limits(img); ax.imshow(img, origin="lower", cmap="gray_r", vmin=lo, vmax=hi); ax.axis("off")
x, y = w.world_to_pixel(SkyCoord(g.ra.values * u.deg, g.dec.values * u.deg)); m = (x > 0) & (y > 0) & (x < img.shape[1]) & (y < img.shape[0])
ax.scatter(x[m], y[m], s=90, facecolors="none", edgecolors="#c0392b", lw=1.1)
L = 120 / abs(hd["PLTSCALE"] * hd["XPIXELSZ"] / 1000); bx, by = 0.04 * img.shape[1], 0.04 * img.shape[0]
ax.plot([bx, bx + L], [by, by], color="black", lw=2.5, solid_capstyle="butt")
ax.text(bx + L / 2, by + 0.02 * img.shape[0], "2′", ha="center", va="bottom", fontsize=15, color="black", path_effects=[pe.withStroke(linewidth=3, foreground="white")])
save(fig, "nuc_densest.png")
print("images written")
