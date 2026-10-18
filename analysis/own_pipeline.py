"""Our own transient search on one full POSS-I red plate scan, to compare with VASCO's catalogue.

1. Read the DSS scan (data/scans/<file>.fits) in tiles; subtract a smooth background; detect sources with photutils
   (3 sigma, at least 3 connected pixels, deblended) and measure centroid, FWHM and elongation.
2. Keep point-like sources: elongation < 1.3 and FWHM no wider than the stellar locus (compact spots kept).
3. Pixel -> sky with the plate's GSC-1 solution (AMDX/AMDY polynomials, as in WCSTools' dsspos).
4. Remove anything with a Gaia DR3 source within 5" (data/scans/gaia_<plate>.parquet from gaia_field.py).
5. Compare with VASCO's transients on the plate: how many we detect at all, and how many survive as Gaia-free.

python analysis/own_pipeline.py XE320 [--partial]"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.convolution import convolve
from astropy.stats import sigma_clipped_stats
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog, deblend_sources, detect_sources, make_2dgaussian_kernel
from scipy.spatial import cKDTree

warnings.simplefilter("ignore")
from scan_io import scan_path

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
f = scan_path(plate)
# read the header from the start of the file and map the pixels directly, so a partly downloaded scan can be used:
# only the complete rows present are processed (pass --partial to allow this)
raw = open(f, "rb").read(2880 * 40)
end = raw.find(b"END" + b" " * 77)
h = fits.Header.fromstring(raw[:end + 80].decode("ascii"))
data0 = (end // 2880 + 1) * 2880
NX, NY = h["NAXIS1"], h["NAXIS2"]
rows_have = min(NY, (f.stat().st_size - data0) // (2 * NX))
if rows_have < NY and "--partial" not in sys.argv:
    raise SystemExit(f"scan incomplete ({rows_have} of {NY} rows); rerun with --partial to use the rows present")
img = np.memmap(f, dtype=">i2", mode="r", offset=data0, shape=(rows_have, NX))
ny, nx = img.shape
pix_arcsec = h["PLTSCALE"] * h["XPIXELSZ"] / 1000


def plate_to_sky(px, py):
    """0-based pixel coordinates -> RA, Dec (deg), GSC-1 plate solution (WCSTools dsspos)."""
    x = px + 1 + h.get("CNPIX1", 0) - 0.5
    y = py + 1 + h.get("CNPIX2", 0) - 0.5
    xm = (h["PPO3"] - x * h["XPIXELSZ"]) / 1000.0
    ym = (y * h["YPIXELSZ"] - h["PPO6"]) / 1000.0
    a = [h.get(f"AMDX{i}", 0.0) for i in range(1, 14)]
    b = [h.get(f"AMDY{i}", 0.0) for i in range(1, 14)]
    r2 = xm * xm + ym * ym
    xi = (a[0] * xm + a[1] * ym + a[2] + a[3] * xm * xm + a[4] * xm * ym + a[5] * ym * ym + a[6] * r2 + a[7] * xm ** 3
          + a[8] * xm * xm * ym + a[9] * xm * ym * ym + a[10] * ym ** 3 + a[11] * xm * r2 + a[12] * xm * r2 * r2)
    eta = (b[0] * ym + b[1] * xm + b[2] + b[3] * ym * ym + b[4] * xm * ym + b[5] * xm * xm + b[6] * r2 + b[7] * ym ** 3
           + b[8] * ym * ym * xm + b[9] * ym * xm * xm + b[10] * xm ** 3 + b[11] * ym * r2 + b[12] * ym * r2 * r2)
    xi, eta = np.radians(xi / 3600), np.radians(eta / 3600)
    ra0 = np.radians(15 * (h["PLTRAH"] + h["PLTRAM"] / 60 + h["PLTRAS"] / 3600))
    sgn = -1 if h["PLTDECSN"].strip().startswith("-") else 1
    dec0 = np.radians(sgn * (h["PLTDECD"] + h["PLTDECM"] / 60 + h["PLTDECS"] / 3600))
    ctan, ccos = np.tan(dec0), np.cos(dec0)
    raoff = np.arctan2(xi / ccos, 1 - eta * ctan)
    ra = np.degrees(raoff + ra0) % 360
    dec = np.degrees(np.arctan(np.cos(raoff) * (eta + ctan) / (1 - eta * ctan)))
    return ra, dec


# ---- detection, tile by tile. Whether sources are high or low values is decided once, from the plate centre: deciding
# per tile fails at the edges, where the unexposed border dominates the tile
c = np.asarray(img[ny // 2 - 500:ny // 2 + 500, nx // 2 - 500:nx // 2 + 500], dtype=np.float32)
SIGN = 1 if c.max() - np.median(c) > np.median(c) - c.min() else -1
T, M = 2000, 40
kernel = make_2dgaussian_kernel(2.0, size=5)
rows = []
for y0 in range(0, ny, T):
    for x0 in range(0, nx, T):
        ya, yb, xa, xb = max(0, y0 - M), min(ny, y0 + T + M), max(0, x0 - M), min(nx, x0 + T + M)
        d = SIGN * np.asarray(img[ya:yb, xa:xb], dtype=np.float32)
        bkg = Background2D(d, 64, filter_size=3, bkg_estimator=MedianBackground())
        sub = d - bkg.background
        conv = convolve(sub, kernel)
        cnoise = sigma_clipped_stats(conv, sigma=3.0)[2]          # noise of the smoothed image, which is what is thresholded
        seg = detect_sources(conv, 3.0 * cnoise, n_pixels=3)
        if seg is None:
            continue
        seg = deblend_sources(conv, seg, n_pixels=3, n_levels=32, contrast=0.005, progress_bar=False)
        cat = SourceCatalog(sub, seg, convolved_data=conv, error=None, background=bkg.background)
        t = cat.to_table(["x_centroid", "y_centroid", "fwhm", "elongation", "segment_flux", "max_value", "area"]).to_pandas()
        t = t.rename(columns={"x_centroid": "xcentroid", "y_centroid": "ycentroid"})
        t["xcentroid"] += xa; t["ycentroid"] += ya
        core = (t.xcentroid >= x0) & (t.xcentroid < x0 + T) & (t.ycentroid >= y0) & (t.ycentroid < y0 + T)
        t["noise"] = float(cnoise)
        rows.append(t[core])
    print(f"  rows up to {min(ny, y0 + T)} of {ny}: {sum(len(r) for r in rows):,} sources", flush=True)
src = pd.concat(rows, ignore_index=True)
src["fwhm_arcsec"] = src.fwhm.astype(float) * pix_arcsec
src["ra"], src["dec"] = plate_to_sky(src.xcentroid.values, src.ycentroid.values)

# ---- Gaia, moved back to the plate's epoch with its proper motions
g = pd.read_parquet(ROOT / f"data/scans/gaia_{plate}.parquet")
pl = pd.read_parquet(ROOT / "data/plates.parquet")
if plate in pl.index:
    dt = (pd.Timestamp("2016-01-01") - pd.Timestamp(pl.loc[plate, "ut"])).days / 365.25   # Gaia DR3 epoch J2016.0
else:                                                    # plates outside VASCO's catalogue: the header's date
    dt = (pd.Timestamp("2016-01-01") - pd.Timestamp(str(h["DATE-OBS"])[:10])).days / 365.25
g["ra"] = g.ra - g.pmra.fillna(0) / 3.6e6 * dt / np.cos(np.radians(g.dec))
g["dec"] = g.dec - g.pmdec.fillna(0) / 3.6e6 * dt


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


gt = cKDTree(xyz(g.ra.values, g.dec.values))
# re-calibrate the GSC-1 solution: match bright sources to Gaia and fit a cubic correction across the plate
bx, by = (src.xcentroid.values - nx / 2) / nx, (src.ycentroid.values - ny / 2) / max(ny, 1)
terms = lambda x, y: np.c_[np.ones_like(x), x, y, x * x, x * y, y * y, x ** 3, x * x * y, x * y * y, y ** 3]
A = terms(bx, by)
bright = src.segment_flux.values >= np.quantile(src.segment_flux.values, 0.6)
for radius in (12.0, 4.0, 2.0):
    d, i = gt.query(xyz(src.ra.values, src.dec.values))
    sep = np.degrees(2 * np.arcsin(d / 2)) * 3600
    m = bright & (sep < radius)
    cosd = np.cos(np.radians(src.dec.values))
    dra = (g.ra.values[i] - src.ra.values) * cosd; ddec = g.dec.values[i] - src.dec.values
    cra = np.linalg.lstsq(A[m], dra[m], rcond=None)[0]; cde = np.linalg.lstsq(A[m], ddec[m], rcond=None)[0]
    src["ra"] = src.ra + (A @ cra) / cosd; src["dec"] = src.dec + A @ cde
    print(f"  astrometric fit on {m.sum():,} stars matched within {radius:.0f} arcsec", flush=True)
d, _ = gt.query(xyz(src.ra.values, src.dec.values))
src["gaia_sep"] = np.degrees(2 * np.arcsin(d / 2)) * 3600
bright = src[(src.gaia_sep < 3)]
print(f"{len(src):,} sources; {np.mean(src.gaia_sep < 2):.0%} have Gaia within 2\"; median offset of matched "
      f"{bright.gaia_sep.median():.2f}\" (checks the plate solution)")

# ---- star-like sources without Gaia counterparts
stellar = src[(src.gaia_sep < 2) & (src.elongation < 1.3)]
lo, hi = np.percentile(stellar.fwhm_arcsec, [2, 98])
starlike = (src.elongation < 1.3) & (src.fwhm_arcsec <= hi)          # no lower limit: compact spots are kept
cand = src[starlike & (src.gaia_sep > 5)]
print(f"stellar FWHM range {lo:.1f}-{hi:.1f}\"; star-like sources {starlike.sum():,}; with no Gaia source within 5\": {len(cand):,}")

# ---- compare with VASCO
v = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
v = v[v.plate == plate]
# keep only VASCO transients inside the scanned area (well inside it, if the scan is partial)
gy, gx = np.mgrid[20:ny - 20:20, 20:nx - 20:20]
gra, gdec = plate_to_sky(gx.ravel().astype(float), gy.ravel().astype(float))
gd, _ = cKDTree(xyz(gra, gdec)).query(xyz(v.ra.values, v.dec.values))
v = v[np.degrees(2 * np.arcsin(gd / 2)) * 3600 < 30]
for name, s in (("all detections", src), ("star-like", src[starlike]), ("our candidates (star-like, no Gaia within 5\")", cand)):
    t = cKDTree(xyz(s.ra.values, s.dec.values))
    dd, _ = t.query(xyz(v.ra.values, v.dec.values))
    sep = np.degrees(2 * np.arcsin(dd / 2)) * 3600
    print(f"VASCO's {len(v)} transients found among {name}: {np.mean(sep < 3):.0%} within 3\" (median distance {np.median(sep):.1f}\")")
out = ROOT / f"data/scans/own_{plate}.parquet"
src.assign(starlike=starlike.values, candidate=(starlike & (src.gaia_sep > 5)).values).to_parquet(out)
print("->", out)
