"""Do our extra candidates show in the ESO DSS-1 cutouts that VASCO searched (the compressed DSS, not the full-resolution
STScI scan)? Fetches a 3x3 arcmin cutout around each candidate passing Solano's cuts other than SPREAD_MODEL, and
measures the same windowed S/N and half-maximum FWHM in both versions of the image.
python analysis/eso_check.py XE320"""
import io
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.convolution import convolve
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import SourceCatalog, detect_sources
from scipy import ndimage

from scan_io import scan_path

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
t = pd.read_parquet(ROOT / f"data/scans/sex_{plate}_cut.parquet")
dxb, dyb = t.bbox_xmax - t.bbox_xmin, t.bbox_ymax - t.bbox_ymin
k = 1.4826
sel = ((t.gaia_sep > 5) & (t.ps1_sep > 5) & (t.near_bright == 0) & (t.snr_win > 30) & ~t.blended
       & ((t.fwhm_px - t.fm).abs() < 2 * k * t.fs) & ((t.elongation - t.em).abs() < 2 * k * t.es)
       & t.fwhm_px.between(2, 7, inclusive="neither") & (t.elongation < 1.3) & ((dxb - dyb).abs() < 2) & (dxb > 1) & (dyb > 1))
c = t[sel].reset_index(drop=True)
print(f"{len(c)} candidates ({c.vasco.sum()} VASCO)")

f = scan_path(plate)
raw = open(f, "rb").read(2880 * 40)
end = raw.find(b"END" + b" " * 77)
h = fits.Header.fromstring(raw[:end + 80].decode("ascii"))
img = np.memmap(f, dtype=">i2", mode="r", offset=(end // 2880 + 1) * 2880, shape=(h["NAXIS2"], h["NAXIS1"]))
FILT = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]], float); FILT /= FILT.sum()
cache = ROOT / "data/scans/eso"; cache.mkdir(exist_ok=True)
URL = "https://archive.eso.org/dss/dss/image"


def fetch(r):
    p = cache / f"{plate}_{r.ra:.5f}_{r.dec:.5f}.fits"
    if not p.exists():
        for attempt in range(5):
            try:
                b = requests.get(URL, params={"ra": r.ra, "dec": r.dec, "equinox": "J2000", "x": 3, "y": 3,
                                              "Sky-Survey": "DSS1", "mime-type": "download-fits"}, timeout=60).content
                fits.getheader(io.BytesIO(b)); p.write_bytes(b); break
            except Exception:
                time.sleep(5)
    return p


def measure(a, x, y):
    """Windowed S/N and half-maximum FWHM of the 5-sigma detection nearest (x, y) in a small, already positive image."""
    a = a - np.median(a)
    rms = sigma_clipped_stats(a, sigma=3.0)[2]
    seg = detect_sources(convolve(a, FILT), 5 * rms, n_pixels=5)
    if seg is None:
        return np.nan, np.nan
    cat = SourceCatalog(a, seg).to_table(columns=["x_centroid", "y_centroid"]).to_pandas()
    r50 = np.asarray(SourceCatalog(a, seg).flux_radius(0.5))
    j = np.argmin(np.hypot(cat.x_centroid - x, cat.y_centroid - y))
    cx, cy = cat.x_centroid[j], cat.y_centroid[j]
    if np.hypot(cx - x, cy - y) > 2:
        return np.nan, np.nan
    yy, xx = np.indices(a.shape)
    sw = max(2 * r50[j] / 2.355, 0.5)
    w = np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sw ** 2)) * (np.hypot(xx - cx, yy - cy) < 7)
    snr = (w * a).sum() / (rms * np.sqrt((w ** 2).sum()))
    ix, iy = int(round(cx)), int(round(cy))
    pk = a[iy - 1:iy + 2, ix - 1:ix + 2].max()
    lab = ndimage.label(a[iy - 7:iy + 8, ix - 7:ix + 8] >= pk / 2)[0]
    return snr, (2 * np.sqrt((lab == lab[7, 7]).sum() / np.pi) if lab[7, 7] else np.nan)


with ThreadPoolExecutor(4) as ex:
    paths = list(ex.map(fetch, [r for _, r in c.iterrows()]))
out = []
for r, p in zip(c.itertuples(), paths):
    e = fits.open(p)[0]
    sign = 1 if e.data.max() - np.median(e.data) > np.median(e.data) - e.data.min() else -1
    x0, y0 = e.header["CNPIX1"], e.header["CNPIX2"]
    ny, nx = e.data.shape
    s = sign * np.asarray(img[y0:y0 + ny, x0:x0 + nx], dtype=float)
    ex_ = sign * e.data.astype(float)
    out.append((*measure(s, r.x - x0, r.y - y0), *measure(ex_, r.x - x0, r.y - y0)))
c[["snr_stsci", "fwhm_stsci", "snr_eso", "fwhm_eso"]] = out
c.to_parquet(ROOT / f"data/scans/eso_check_{plate}.parquet")
for nm, m in (("VASCO", c.vasco), ("others", ~c.vasco)):
    x = c[m]
    print(f"{nm:7s} n {len(x):4d}  found in ESO cutout {x.snr_eso.notna().mean():.0%}  "
          f"S/N STScI {x.snr_stsci.median():.0f} -> ESO {x.snr_eso.median():.0f}  "
          f"ESO S/N > 30 {np.mean(x.snr_eso > 30):.0%}  FWHM {x.fwhm_stsci.median():.2f} -> {x.fwhm_eso.median():.2f} px")
