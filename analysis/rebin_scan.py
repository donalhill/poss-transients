"""Resample a 15-micron DSS scan (1.0"/pixel) to the 25-micron pixels (1.7") of the other scans, so that the same
pixel-based cuts apply. Each block of 5 x 5 pixels becomes 3 x 3, by area-weighted averaging; the header's pixel size
is scaled to match, so the plate solution still holds.
python analysis/rebin_scan.py xe926_a2ti_01_00   -> data/scans/xe926_a2ti_01_00_r25.fits"""
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

from scan_io import header

ROOT = Path(__file__).resolve().parents[1]
name = sys.argv[1]
src = ROOT / f"data/scans/{name}.fits"
h, off = header(src)
NX, NY = h["NAXIS1"], h["NAXIS2"]
img = np.memmap(src, dtype=">i2", mode="r", offset=off, shape=(NY, NX))
# output pixel j covers input [5j/3, 5(j+1)/3): weights within each block of 5 inputs -> 3 outputs, as averages
W = np.array([[1, 2 / 3, 0, 0, 0], [0, 1 / 3, 1, 1 / 3, 0], [0, 0, 0, 2 / 3, 1]]) * 0.6
nx, ny = NX // 5 * 3, NY // 5 * 3
out = np.empty((ny, nx), dtype=np.int16)
for b in range(NY // 5):
    a = np.asarray(img[5 * b:5 * b + 5, :NX // 5 * 5], dtype=np.float64)       # 5 rows
    a = np.einsum("ij,jk->ik", W, a)                                          # 3 rows
    a = np.einsum("rbj,ij->rbi", a.reshape(3, -1, 5), W).reshape(3, -1)       # 3 x (nx)
    out[3 * b:3 * b + 3] = np.round(a).astype(np.int16)
for k in ("XPIXELSZ", "YPIXELSZ"):
    h[k] = h[k] * 5 / 3
for k in ("CNPIX1", "CNPIX2"):
    h[k] = h.get(k, 0) * 3 / 5
h["NAXIS1"], h["NAXIS2"] = nx, ny
h["HISTORY"] = f"Resampled 5x5 -> 3x3 from {name}.fits (rebin_scan.py)"
fits.PrimaryHDU(out, header=h).writeto(ROOT / f"data/scans/{name}_r25.fits", overwrite=True)
print(f"{name}: {NX}x{NY} -> {nx}x{ny}")
