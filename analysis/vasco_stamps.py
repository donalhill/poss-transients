"""Is there anything on the scan at VASCO's transient positions?

For each VASCO transient on the plate, the significance of the brightest pixel within 2.5" of its position, above the
local background (median and robust spread of an annulus 15-25 pixels out). The same for two comparison sets: faint
Gaia stars (G 19-20.5, moved to the plate epoch) and random positions. Also saves image stamps of a sample.
Positions are mapped to pixels with a cubic fit to our own detections (data/scans/own_<plate>.parquet).

python analysis/vasco_stamps.py XE320"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scan_io import scan_path

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
f = scan_path(plate)
img = fits.open(f, memmap=True)[0].data
ny, nx = img.shape
src = pd.read_parquet(ROOT / f"data/scans/own_{plate}.parquet")
ra0, dec0 = src.ra.median(), src.dec.median()


def tangent(ra, dec):
    ra, dec, a0, d0 = map(np.radians, (ra, dec, ra0, dec0))
    c = np.sin(d0) * np.sin(dec) + np.cos(d0) * np.cos(dec) * np.cos(ra - a0)
    return np.degrees(np.cos(dec) * np.sin(ra - a0) / c), np.degrees((np.cos(d0) * np.sin(dec) - np.sin(d0) * np.cos(dec) * np.cos(ra - a0)) / c)


terms = lambda x, y: np.c_[np.ones_like(x), x, y, x * x, x * y, y * y, x ** 3, x * x * y, x * y * y, y ** 3]
fit = src[src.gaia_sep < 1.0]                                            # well-measured stars
u, v = tangent(fit.ra.values, fit.dec.values)
cx = np.linalg.lstsq(terms(u, v), fit.xcentroid.values, rcond=None)[0]
cy = np.linalg.lstsq(terms(u, v), fit.ycentroid.values, rcond=None)[0]
res = np.hypot(terms(u, v) @ cx - fit.xcentroid.values, terms(u, v) @ cy - fit.ycentroid.values)
print(f"sky -> pixel fit on {len(fit):,} stars: median residual {np.median(res):.2f} px")


def to_pix(ra, dec):
    u, v = tangent(np.asarray(ra), np.asarray(dec))
    return terms(u, v) @ cx, terms(u, v) @ cy


img_sign = 1 if np.median(img[6000:6200, 6000:6200]) - img[6000:6200, 6000:6200].min() < img[6000:6200, 6000:6200].max() - np.median(img[6000:6200, 6000:6200]) else -1
Y, X = np.mgrid[-25:26, -25:26]; R = np.hypot(X, Y)


def peak_snr(x, y):
    xi, yi = int(round(x)), int(round(y))
    if not (30 <= xi < nx - 30 and 30 <= yi < ny - 30):
        return np.nan
    s = img_sign * img[yi - 25:yi + 26, xi - 25:xi + 26].astype(float)
    ann = s[(R >= 15) & (R <= 25)]
    med = np.median(ann); sig = 1.4826 * np.median(np.abs(ann - med))
    return (s[R <= 1.5].max() - med) / sig


v = pd.read_parquet(ROOT / "data/transients_with_plates.parquet"); v = v[v.plate == plate]
g = pd.read_parquet(ROOT / f"data/scans/gaia_{plate}.parquet")
dt = (pd.Timestamp("2016-01-01") - pd.Timestamp(pd.read_parquet(ROOT / "data/plates.parquet").loc[plate, "ut"])).days / 365.25
g["ra"] = g.ra - g.pmra.fillna(0) / 3.6e6 * dt / np.cos(np.radians(g.dec)); g["dec"] = g.dec - g.pmdec.fillna(0) / 3.6e6 * dt
rng = np.random.default_rng(1)
sets = {"VASCO transients": (v.ra.values, v.dec.values)}
for lo, hi in ((17, 18), (19, 20), (20, 20.5), (20.5, 21)):
    gg = g[g.phot_g_mean_mag.between(lo, hi)].sample(min(400, (g.phot_g_mean_mag.between(lo, hi)).sum()), random_state=1)
    sets[f"Gaia G {lo}-{hi}"] = (gg.ra.values, gg.dec.values)
rx, ry = rng.uniform(300, nx - 300, 400), rng.uniform(300, ny - 300, 400)
out = {}
xlo, xhi, ylo, yhi = fit.xcentroid.min(), fit.xcentroid.max(), fit.ycentroid.min(), fit.ycentroid.max()
inside = lambda x, y: (x > xlo + 50) & (x < xhi - 50) & (y > ylo + 50) & (y < yhi - 50)   # where the mapping was fitted
for name, (ra, dec) in sets.items():
    x, y = to_pix(ra, dec)
    out[name] = np.array([peak_snr(a, b) if inside(a, b) else np.nan for a, b in zip(x, y)])
rx, ry = rng.uniform(xlo + 50, xhi - 50, 400), rng.uniform(ylo + 50, yhi - 50, 400)
out["random positions"] = np.array([peak_snr(a, b) for a, b in zip(rx, ry)])
for name, s in out.items():
    s = s[np.isfinite(s)]
    print(f"{name:20s}: n {len(s):4d}  median peak {np.median(s):5.1f} sigma;  share above 5 sigma {np.mean(s > 5):.0%}, above 3 sigma {np.mean(s > 3):.0%}")

# stamps: 24 VASCO transients and 8 faint Gaia stars, 40" across
x, y = to_pix(v.ra.values, v.dec.values); keep = inside(x, y)
x, y, vs = x[keep], y[keep], out["VASCO transients"][keep]
gx, gy = to_pix(*sets["Gaia G 20-20.5"]); gk = inside(gx, gy)
gx, gy, gs = gx[gk], gy[gk], out["Gaia G 20-20.5"][gk]
fig, axes = plt.subplots(4, 8, figsize=(16, 8.4))
items = [("VASCO", a, b, s) for a, b, s in list(zip(x, y, vs))[:24]] + [("Gaia 20-20.5", a, b, s) for a, b, s in list(zip(gx, gy, gs))[:8]]
for ax, (lab, a, b, s) in zip(axes.ravel(), items):
    xi, yi = int(round(a)), int(round(b))
    st = img_sign * img[max(0, yi - 12):yi + 13, max(0, xi - 12):xi + 13].astype(float)
    lo_, hi_ = np.percentile(st, [5, 99.5])
    ax.imshow(st, origin="lower", cmap="gray_r", vmin=lo_, vmax=hi_)
    ax.add_patch(plt.Circle((12, 12), 3, fill=False, color="#e05a3a", lw=1))
    ax.set_title(f"{lab}  {s:.1f}σ", fontsize=8); ax.axis("off")
fig.tight_layout()
fig.savefig(ROOT / f"logs/stamps_{plate}.png", dpi=110)
print("stamps ->", ROOT / f"logs/stamps_{plate}.png")
