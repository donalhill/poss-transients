"""How deep each plate goes and how sharp its star images are, measured on the scans: the weather, the emulsion and the
processing all act on a plate's transients through these.

Two 10' DSS cutouts per plate, 1 deg from its centre towards opposite corners (the scan of the plate itself; a cutout
served from a neighbouring plate is left out). On each, sources are found at 5 sigma above the local sky and matched to
Gaia DR3 stars within 4". Depth: the share of Gaia stars with 18 < G < 20 that are found (completeness). Sharpness:
the median width (FWHM) of up to 40 isolated, unsaturated detections that match a Gaia star. Detections without a Gaia
star (specks, transients) enter neither, so a plate's own transients cannot change its depth or sharpness.

The DSS and Gaia services are slow, so plates already measured are kept.
python vasco/plate_depth.py -> data/plate_depth.csv"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from photutils.detection import DAOStarFinder
from scipy.ndimage import median_filter

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).resolve().parents[1] / "analysis")]
import vasco_sample as sample                                                  # noqa: E402
from geometry import deproject                                                 # noqa: E402
from stars import fit                                                          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = ROOT / "data/cutouts/depth"; CUT.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data/plate_depth.csv"
SIZE = 10.0                                                                    # arcminutes
OFFSETS = [(1.0, 1.0), (-1.0, -1.0)]                                           # degrees from the plate centre
PIX = 1.7                                                                      # arcseconds per pixel on the DSS scans
gaia = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag"], column_filters={"Gmag": "<20"}, row_limit=-1)


def cutout(q, k, c):
    f = CUT / f"{q}_{k}.fits"
    if not f.exists():
        r = requests.get(f"https://archive.stsci.edu/cgi-bin/dss_search?v=poss1_red&r={c.ra.deg:.5f}&d={c.dec.deg:.5f}"
                         f"&e=J2000&h={SIZE}&w={SIZE}&f=fits&c=none", timeout=300)
        r.raise_for_status(); f.write_bytes(r.content)
    h = fits.getheader(f)
    return fits.getdata(f).astype(float) if h["REGION"].strip() == q else None


def stars_in(q, k, c):
    """Gaia DR3 stars brighter than G = 20 in the cutout (RA, Dec, G), cached."""
    f = CUT / f"{q}_{k}_gaia.csv"
    if not f.exists():
        t = gaia.query_region(c, width=SIZE * u.arcmin, height=SIZE * u.arcmin, catalog="I/355/gaiadr3")
        (t[0].to_pandas() if len(t) else pd.DataFrame(columns=["RA_ICRS", "DE_ICRS", "Gmag"]))[["RA_ICRS", "DE_ICRS", "Gmag"]].to_csv(f, index=False)
    return pd.read_csv(f)


def measure(q, p):
    """Depth: the share of Gaia stars with 18 < G < 20 found within 4" by a detection. Sharpness: the median width of
    isolated, unsaturated detections that match a Gaia star. Neither counts detections without a Gaia star, so specks
    and transients on a plate cannot change them."""
    faint, found_faint, widths = 0, 0, []
    for k, (dx, dy) in enumerate(OFFSETS):
        ra, dec = deproject(np.array([dx]), np.array([dy]), p.ra_c, p.dec_c); c = SkyCoord(ra[0] * u.deg, dec[0] * u.deg)
        try:
            R = cutout(q, k, c)
            g = stars_in(q, k, c)
        except Exception:
            continue
        if R is None or g.empty:
            continue
        w = WCS(fits.getheader(CUT / f"{q}_{k}.fits"))
        im = R - median_filter(R, size=31)
        _, med, sd = sigma_clipped_stats(im)
        found = DAOStarFinder(fwhm=2.5, threshold=5 * sd)(im - med)
        if found is None:
            continue
        xc, yc = ("x_centroid", "y_centroid") if "x_centroid" in found.colnames else ("xcentroid", "ycentroid")
        x, y = np.asarray(found[xc], float), np.asarray(found[yc], float)
        gx, gy = w.world_to_pixel(SkyCoord(g.RA_ICRS.values * u.deg, g.DE_ICRS.values * u.deg))
        inside = (gx > 8) & (gy > 8) & (gx < R.shape[1] - 9) & (gy < R.shape[0] - 9)
        d = np.hypot(gx[:, None] - x[None, :], gy[:, None] - y[None, :])
        near = d.min(1) * PIX < 4 if len(x) else np.zeros(len(gx), bool)
        f_ = inside & (g.Gmag.values > 18) & (g.Gmag.values < 20)
        faint += int(f_.sum()); found_faint += int((f_ & near).sum())
        # sharpness on detections that match a Gaia star, isolated and unsaturated, brightest first
        match = (d.min(0) * PIX < 4) if len(gx) else np.zeros(len(x), bool)
        sat = 0.8 * R.max() - np.median(R)
        dd = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :]); np.fill_diagonal(dd, np.inf)
        ok = match & (dd.min(1) > 12) & (im[np.round(y).astype(int), np.round(x).astype(int)] < sat)
        for i in np.argsort(-np.asarray(found["peak"], float))[:300]:
            if ok[i] and len(widths) < 40 * (k + 1):
                fw = fit(im, x[i], y[i])
                if fw and fw["fwhm"] < 10:
                    widths.append(fw["fwhm"])
    return {"plate": q, "gaia_faint": faint, "found_faint": found_faint, "completeness": found_faint / faint if faint else np.nan,
            "fwhm_arcsec": float(np.median(widths) * PIX) if widths else np.nan}


if __name__ == "__main__":
    _, plates = sample.ours()
    known = pd.read_csv(OUT).set_index("plate") if OUT.exists() else pd.DataFrame()
    if "completeness" not in known.columns:
        known = pd.DataFrame()
    todo = [q for q in plates.index if q not in known.index]
    with ThreadPoolExecutor(8) as ex:
        new = pd.DataFrame(list(ex.map(lambda q: measure(q, plates.loc[q]), todo))).set_index("plate") if todo else pd.DataFrame()
    d = pd.concat([known, new]).loc[plates.index]
    d.to_csv(OUT, index_label="plate")
    print(d.describe().round(3).to_string())
