"""Shared helpers for the plate scans: which file to use for a plate, its header, and the plate centre in pixels."""
import os
from pathlib import Path

import pandas as pd

from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]


def scan_path(plate):
    """The 25-micron scan, or a 15-micron scan resampled to 25 microns (rebin_scan.py)."""
    d = ROOT / "data/scans"
    return (sorted(d.glob(f"{plate.lower()}_*_r25.fits")) + sorted(d.glob(f"{plate.lower()}_*_00_00.fits")))[0]


def header(path):
    with open(path, "rb") as fh:
        raw = fh.read(2880 * 40)
    end = raw.find(b"END" + b" " * 77)
    return fits.Header.fromstring(raw[:end + 80].decode("ascii")), (end // 2880 + 1) * 2880


def plate_header(plate):
    """A plate's scan header, from a small cache in data/headers (the scans themselves can be deleted after processing,
    or evicted by iCloud). Filled from the local scan if it is fully on disk, else fetched from STScI (115 kB)."""
    cache = ROOT / f"data/headers/{plate.upper()}.hdr"
    if cache.exists():
        return fits.Header.fromstring(cache.read_text())
    d = ROOT / "data/scans"
    local = sorted(d.glob(f"{plate.lower()}_*_r25.fits")) + sorted(d.glob(f"{plate.lower()}_*_00_00.fits"))
    if local and not os.stat(local[0]).st_flags & 0x40000000:         # SF_DATALESS: evicted to iCloud
        h, _ = header(local[0])
    else:
        import requests
        name = pd.read_csv(ROOT / "data/xe_scan_centres.csv").set_index("region").file
        f = name[plate.upper()] if isinstance(name[plate.upper()], str) else sorted(name[plate.upper()])[0]
        raw = requests.get(f"https://archive.stsci.edu/missions/dss/FitsArchive/XEsurvey/{f}.fits",
                           headers={"Range": "bytes=0-115199"}, timeout=60).content
        h = fits.Header.fromstring(raw[:raw.find(b"END" + b" " * 77) + 80].decode("ascii"))
        if h["NAXIS1"] > 20000:                                        # 15-micron scan, used resampled to 25 microns
            for k in ("XPIXELSZ", "YPIXELSZ"):
                h[k] = h[k] * 5 / 3
            for k in ("CNPIX1", "CNPIX2"):
                h[k] = h.get(k, 0) * 3 / 5
            h["NAXIS1"], h["NAXIS2"] = h["NAXIS1"] // 5 * 3, h["NAXIS2"] // 5 * 3
    cache.parent.mkdir(exist_ok=True); cache.write_text(h.tostring())
    return h


def centre(plate):
    """Pixel position (0-based) of the plate centre, from the plate solution's PPO3/PPO6 (WCSTools dsspos)."""
    h = plate_header(plate)
    return (h["PPO3"] / h["XPIXELSZ"] - h.get("CNPIX1", 0) - 0.5,
            h["PPO6"] / h["YPIXELSZ"] - h.get("CNPIX2", 0) - 0.5)


def pixel_arcsec(plate):
    h = plate_header(plate)
    return h["PLTSCALE"] * h["XPIXELSZ"] / 1000


def plate_to_sky(h, px, py):
    """0-based pixel coordinates -> RA, Dec (deg), from the header's GSC-1 plate solution (WCSTools dsspos)."""
    import numpy as np
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
    return np.degrees(raoff + ra0) % 360, np.degrees(np.arctan(np.cos(raoff) * (eta + ctan) / (1 - eta * ctan)))
