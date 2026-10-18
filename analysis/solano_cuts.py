"""Re-measure one plate the way SExtractor/PSFEx would, so that Solano et al. (2022)'s full list of cuts can be applied:
detection at 5 sigma of the local background rms on a 3x3-filtered image (DETECT_THRESH 5, DETECT_MINAREA 5, deblending
32 levels / 0.005), windowed S/N (SNR_WIN), the blended flag (FLAGS 2), half-maximum FWHM, bounding box (XMIN..YMAX),
and SPREAD_MODEL against a PSF built per 2000-pixel tile from stars chosen as PSFEx does (S/N > 20, FWHM within 20%
of the sample's median, ellipticity < 0.3, unblended). Each source then borrows the Gaia, Pan-STARRS, AllWISE, and
bright-star information of the nearest source in own_<plate>.parquet.
python analysis/solano_cuts.py XE320   -> data/scans/sex_<plate>.parquet"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.convolution import convolve
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog, deblend_sources, detect_sources
from scipy import ndimage
from scipy.spatial import cKDTree

from scan_io import scan_path

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
f = scan_path(plate)
raw = open(f, "rb").read(2880 * 40)
end = raw.find(b"END" + b" " * 77)
h = fits.Header.fromstring(raw[:end + 80].decode("ascii"))
data0 = (end // 2880 + 1) * 2880
NX, NY = h["NAXIS1"], h["NAXIS2"]
img = np.memmap(f, dtype=">i2", mode="r", offset=data0, shape=(NY, NX))
c = np.asarray(img[NY // 2 - 500:NY // 2 + 500, NX // 2 - 500:NX // 2 + 500], dtype=np.float32)
SIGN = 1 if c.max() - np.median(c) > np.median(c) - c.min() else -1

T, M, R = 2000, 40, 7                                   # tile, margin, half-size of measurement stamps
FILT = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], float); FILT /= FILT.sum()   # SExtractor's default.conv
OS = 8                                                  # oversampling for the SPREAD_MODEL exponential
yy, xx = np.mgrid[-R:R + 1, -R:R + 1]


def stamp(a, x, y):
    ix, iy = int(round(x)), int(round(y))
    if iy - R < 0 or ix - R < 0 or iy + R + 1 > a.shape[0] or ix + R + 1 > a.shape[1]:
        return None, 0, 0
    return a[iy - R:iy + R + 1, ix - R:ix + R + 1], x - ix, y - iy


def models(psf, fwhm):
    """PSF and PSF convolved with an exponential of scale FWHM/16, both through the same oversampled path."""
    z = ndimage.zoom(psf, OS, order=3) / OS ** 2
    k = np.arange(-4 * OS, 4 * OS + 1) / OS
    r = np.hypot(*np.meshgrid(k, k))
    e = np.exp(-r / (fwhm / 16)); e /= e.sum()
    g = ndimage.convolve(z, e, mode="constant")
    n = psf.shape[0]
    binned = lambda a: a.reshape(n, OS, n, OS).sum(axis=(1, 3))
    phi, G = binned(z), binned(g)
    return phi / phi.sum(), G / G.sum()


rows = []
for y0 in range(0, NY, T):
    for x0 in range(0, NX, T):
        ya, yb, xa, xb = max(0, y0 - M), min(NY, y0 + T + M), max(0, x0 - M), min(NX, x0 + T + M)
        d = SIGN * np.asarray(img[ya:yb, xa:xb], dtype=np.float32)
        bkg = Background2D(d, 64, filter_size=3, bkg_estimator=MedianBackground())
        sub = d - bkg.background
        rmap = bkg.background_rms                                # local noise, per 64-pixel mesh, as SExtractor's map
        rms = float(np.median(rmap))
        conv = convolve(sub, FILT)
        seg0 = detect_sources(conv, 5.0 * rmap, n_pixels=5)
        if seg0 is None:
            continue
        seg = deblend_sources(conv, seg0, n_pixels=5, n_levels=32, contrast=0.005, progress_bar=False)
        cat = SourceCatalog(sub, seg, convolved_data=conv)
        t = cat.to_table(["label", "x_centroid", "y_centroid", "bbox_xmin", "bbox_xmax", "bbox_ymin", "bbox_ymax",
                          "elongation", "max_value"]).to_pandas()
        t["r50"] = np.asarray(cat.fluxfrac_radius(0.5))
        # blended (FLAGS 2): the source's parent segment was split by the deblender
        parent = seg0.data[seg.data > 0]; child = seg.data[seg.data > 0]
        pairs = pd.DataFrame({"p": parent, "c": child}).drop_duplicates()
        nkids = pairs.groupby("p").c.nunique()
        t["blended"] = t.label.map(pairs.set_index("c").p.map(nkids)).fillna(1).values > 1
        core = (t.x_centroid >= x0 - xa) & (t.x_centroid < x0 - xa + T) & (t.y_centroid >= y0 - ya) & (t.y_centroid < y0 - ya + T)
        t = t[core].reset_index(drop=True)
        # windowed S/N, half-maximum FWHM, centred stamps
        snr, fw, cut, off = [], [], [], []
        for x, y, r50 in zip(t.x_centroid, t.y_centroid, t.r50):
            s, dx, dy = stamp(sub, x, y)
            if s is None or not np.isfinite(r50):
                snr.append(np.nan); fw.append(np.nan); cut.append(None); off.append((0, 0)); continue
            sw = max(2 * r50 / 2.355, 0.5)
            w = np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2 * sw ** 2))
            snr.append((w * s).sum() / (rmap[int(round(y)), int(round(x))] * np.sqrt((w ** 2).sum())))
            pk = s[R - 1:R + 2, R - 1:R + 2].max()
            lab = ndimage.label(s >= pk / 2)[0]
            fw.append(2 * np.sqrt((lab == lab[R, R]).sum() / np.pi) if lab[R, R] else np.nan)
            cut.append(s); off.append((dx, dy))
        t["snr_win"], t["fwhm_px"] = snr, fw
        # PSF for the tile from PSFEx-like stars: S/N > 20, unblended, ellipticity < 0.3, FWHM within 20% of median
        ok = (t.snr_win > 20) & ~t.blended & (1 - 1 / t.elongation < 0.3) & t.fwhm_px.between(2, 10)
        med = t.fwhm_px[ok].median()
        ok &= (t.fwhm_px - med).abs() < 0.2 * med
        stack, wts = [], []
        for i in np.flatnonzero(ok):
            s, (dx, dy) = cut[i], off[i]
            sh = ndimage.shift(s, (-dy, -dx), order=3, mode="nearest")
            stack.append(sh / sh.sum()); wts.append(t.snr_win[i] ** 2)
        psf = np.average(stack, axis=0, weights=wts) if len(stack) > 20 else None
        spread = np.full(len(t), np.nan)
        if psf is not None:
            phi0, G0 = models(psf, med)
            for i in range(len(t)):
                if cut[i] is None:
                    continue
                dx, dy = off[i]
                phi = ndimage.shift(phi0, (dy, dx), order=3); G = ndimage.shift(G0, (dy, dx), order=3)
                x = cut[i]
                spread[i] = (G * x).sum() / (phi * x).sum() - (G * phi).sum() / (phi * phi).sum()
        t["spread_model"] = spread
        t["n_psf"] = len(stack)
        t["tile"] = f"{y0 // T}_{x0 // T}"
        t["x"] = t.x_centroid + xa; t["y"] = t.y_centroid + ya
        rows.append(t.drop(columns=["label", "x_centroid", "y_centroid"]))
    print(f"  rows up to {min(NY, y0 + T)} of {NY}: {sum(len(r) for r in rows):,} sources", flush=True)

t = pd.concat(rows, ignore_index=True)
s = pd.read_parquet(ROOT / f"data/scans/own_{plate}.parquet")
dist, i = cKDTree(s[["xcentroid", "ycentroid"]].values).query(t[["x", "y"]].values)
for col in ("ra", "dec", "gaia_sep", "ps1_sep", "wise_sep", "near_bright"):
    t[col] = np.where(dist < 2, s[col].values[i], np.nan)
t["matched"] = dist < 2
t.to_parquet(ROOT / f"data/scans/sex_{plate}.parquet")
print(f"{len(t):,} sources, {t.matched.mean():.1%} matched to own_{plate}")
