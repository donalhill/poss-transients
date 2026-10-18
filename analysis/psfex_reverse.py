"""Which stars did PSFEx use? Rebuild each tile's PSF in several ways and see which one reproduces the edge that
SPREAD_MODEL > -0.002 must leave in VASCO's catalogue (every catalogued transient passed it).
  psfex: PSFEx defaults - S/N > 20, unblended, ellipticity < 0.3, FWHM within 20% of the sample's modal FWHM, pixels
         weighted by 1 / ((rms/F)^2 + (PSF_ACCURACY p)^2) with PSF_ACCURACY = 0.01, so stars above S/N ~100 stop
         gaining weight
  snr_a_b: unweighted mean of unblended round stars with S/N between a and b
python analysis/psfex_reverse.py XE320   -> data/scans/spread_<plate>.parquet"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.background import Background2D, MedianBackground
from scipy import ndimage

from scan_io import scan_path

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
sys.argv = sys.argv[:2]
t = pd.read_parquet(ROOT / f"data/scans/sex_{plate}_cut.parquet")
f = scan_path(plate)
raw = open(f, "rb").read(2880 * 40)
end = raw.find(b"END" + b" " * 77)
h = fits.Header.fromstring(raw[:end + 80].decode("ascii"))
NX, NY = h["NAXIS1"], h["NAXIS2"]
img = np.memmap(f, dtype=">i2", mode="r", offset=(end // 2880 + 1) * 2880, shape=(NY, NX))
c = np.asarray(img[NY // 2 - 500:NY // 2 + 500, NX // 2 - 500:NX // 2 + 500], dtype=np.float32)
SIGN = 1 if c.max() - np.median(c) > np.median(c) - c.min() else -1

T, M, R, OS = 2000, 40, 7, 8
yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
BANDS = [(20, 30), (30, 50), (50, 100), (100, 200), (200, 500)]


def models(psf, fwhm):
    z = ndimage.zoom(psf, OS, order=3) / OS ** 2
    k = np.arange(-4 * OS, 4 * OS + 1) / OS
    e = np.exp(-np.hypot(*np.meshgrid(k, k)) / (fwhm / 16)); e /= e.sum()
    g = ndimage.convolve(z, e, mode="constant")
    n = psf.shape[0]
    binned = lambda a: a.reshape(n, OS, n, OS).sum(axis=(1, 3))
    phi, G = binned(z), binned(g)
    return phi / phi.sum(), G / G.sum()


out = []
for tile, tt in t.groupby("tile"):
    ty, tx = map(int, tile.split("_"))
    y0, x0 = ty * T, tx * T
    ya, yb, xa, xb = max(0, y0 - M), min(NY, y0 + T + M), max(0, x0 - M), min(NX, x0 + T + M)
    d = SIGN * np.asarray(img[ya:yb, xa:xb], dtype=np.float32)
    sub = d - Background2D(d, 64, filter_size=3, bkg_estimator=MedianBackground()).background
    rms = sigma_clipped_stats(sub, sigma=3.0)[2]
    tt = tt[tt.snr_win > 20].copy()
    cuts, offs, keep = [], [], []
    for x, y in zip(tt.x - xa, tt.y - ya):
        ix, iy = int(round(x)), int(round(y))
        ok = R <= iy < sub.shape[0] - R and R <= ix < sub.shape[1] - R
        keep.append(ok)
        cuts.append(sub[iy - R:iy + R + 1, ix - R:ix + R + 1] if ok else None); offs.append((x - ix, y - iy))
    tt["ok"] = keep
    centred = {}

    def centre(i):
        if i not in centred:
            dx, dy = offs[i]
            s = ndimage.shift(cuts[i], (-dy, -dx), order=3, mode="nearest")
            centred[i] = (s / s.sum(), s.sum())
        return centred[i]

    round_ = tt.ok & ~tt.blended & (1 - 1 / tt.elongation < 0.3) & (tt.gaia_sep < 1.5)
    psfs = {}
    # PSFEx default selection and weighting
    fw = 2 * tt.r50
    sel = round_ & fw.between(2, 10)
    hist, edges = np.histogram(fw[sel], bins=np.arange(2, 10.05, 0.1))
    mode = edges[np.argmax(hist)] + 0.05
    sel &= fw.between(mode / 1.2, mode * 1.2)
    idx = np.flatnonzero(sel.values)
    if len(idx) > 20:
        stack = np.array([centre(i)[0] for i in idx]); flux = np.array([centre(i)[1] for i in idx])
        p = np.median(stack, axis=0)
        w = 1 / ((rms / flux)[:, None, None] ** 2 + (0.01 * p[None]) ** 2)
        psfs["psfex"] = ((stack * w).sum(0) / w.sum(0), mode)
    for a, b in BANDS:
        idx = np.flatnonzero((round_ & tt.snr_win.between(a, b)).values)
        if len(idx) > 20:
            psfs[f"snr_{a}_{b}"] = (np.mean([centre(i)[0] for i in idx], axis=0), np.median(fw.values[idx]))
    res = {k: np.full(len(tt), np.nan) for k in psfs}
    for k, (psf, fwhm) in psfs.items():
        phi0, G0 = models(psf, fwhm)
        for i in np.flatnonzero(tt.ok.values):
            dx, dy = offs[i]
            phi = ndimage.shift(phi0, (dy, dx), order=3); G = ndimage.shift(G0, (dy, dx), order=3)
            x = cuts[i]
            res[k][i] = (G * x).sum() / (phi * x).sum() - (G * phi).sum() / (phi * phi).sum()
    for k, v in res.items():
        tt[f"sp_{k}"] = v
    tt["psfex_mode_fwhm"] = mode
    out.append(tt)
    print(f"tile {tile}: PSFEx-like FWHM mode {mode:.2f} px, {sel.sum()} stars", flush=True)
pd.concat(out).to_parquet(ROOT / f"data/scans/spread_{plate}.parquet")
