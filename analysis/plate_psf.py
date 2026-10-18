"""Image width (FWHM) of unsaturated stars on the DSS scans of POSS-I red plates, from the scan cutouts already on disk
(data/cutouts), against distance from the plate centre. Sets the two-hit resolution used in alignments_tracks.py.
python analysis/plate_psf.py -> data/plate_psf.csv"""
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import maximum_filter
from scipy.optimize import curve_fit

warnings.simplefilter("ignore", FITSFixedWarning)
ROOT = Path(__file__).resolve().parents[1]
plates = pd.read_parquet(ROOT / "data/plates.parquet")
gauss = lambda XY, A, x0, y0, s, b: A * np.exp(-((XY[0] - x0) ** 2 + (XY[1] - y0) ** 2) / (2 * s * s)) + b


def tangent(ra, dec, ra0, dec0):
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    c = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    return (np.degrees(np.cos(dec) * np.sin(ra - ra0) / c),
            np.degrees((np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / c))


rows, seen = [], set()
for f in sorted(glob.glob(str(ROOT / "data/cutouts/*.fits"))):
    h = fits.getheader(f); img = fits.getdata(f).astype(float)
    plate = str(h.get("REGION", "")).strip()
    if plate not in plates.index or img.ndim != 2:
        continue
    ra, dec = WCS(h).all_pix2world(img.shape[1] / 2, img.shape[0] / 2, 0)
    key = (plate, round(float(ra), 2), round(float(dec), 2))
    if key in seen:
        continue
    seen.add(key)
    x, y = tangent(float(ra), float(dec), plates.loc[plate, "ra_c"], plates.loc[plate, "dec_c"])
    pix = abs(h["PLTSCALE"] * h["XPIXELSZ"] / 1000)
    a = img - np.median(img)
    if abs(a.min()) > abs(a.max()):
        a = -a
    noise = np.std(a[a < np.percentile(a, 80)])
    peaks = (a == maximum_filter(a, 9)) & (a > 8 * noise) & (a < 0.9 * np.percentile(a, 99.9))   # no saturated stars
    Y, X = np.mgrid[-6:7, -6:7]; fw = []
    for py, px in zip(*np.where(peaks)):
        if 6 <= py < a.shape[0] - 6 and 6 <= px < a.shape[1] - 6:
            cut = a[py - 6:py + 7, px - 6:px + 7]
            try:
                p, _ = curve_fit(gauss, (X.ravel(), Y.ravel()), cut.ravel(), p0=[cut.max(), 0, 0, 1.2, 0], maxfev=4000)
                if 0.5 < abs(p[3]) < 5 and abs(p[1]) < 2 and abs(p[2]) < 2:
                    fw.append(2.355 * abs(p[3]) * pix)
            except RuntimeError:
                pass
    if len(fw) >= 10:
        rows.append({"plate": plate, "r_deg": max(abs(x), abs(y)), "n_stars": len(fw), "fwhm": np.median(fw)})
df = pd.DataFrame(rows)
df.to_csv(ROOT / "data/plate_psf.csv", index=False)
print(f"{len(df)} cutouts on {df.plate.nunique()} plates, {df.n_stars.sum():,} stars: median FWHM {df.fwhm.median():.1f}\" "
      f"(16-84% of cutouts {df.fwhm.quantile(.16):.1f}-{df.fwhm.quantile(.84):.1f}\")")
for lo, hi in ((0, 1), (1, 2), (2, 2.5), (2.5, 3.3)):
    m = (df.r_deg >= lo) & (df.r_deg < hi)
    if m.any():
        print(f"  {lo}-{hi} deg from centre: {m.sum():3d} cutouts, median FWHM {df.fwhm[m].median():.1f}\"")
