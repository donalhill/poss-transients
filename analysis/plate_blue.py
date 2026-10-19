"""Specks on the red and blue plates of each field: sources with no star behind them, counted on the scans themselves.

For each red plate, the same two 10' patches as plate_depth.py (1 deg from the centre towards opposite corners), on
the red plate and on the blue plate of the same field. On each: sources found at 5 sigma above the local sky, and how
many of them have no Gaia DR3 star (G < 21) within 5". The red count repeats, from the scans, what VASCO's catalogue
counts; the blue count asks whether the same plates' blue partners hold extra specks too. The blue plate's archive
date is kept, to check it was taken the same night as the red one.

The DSS and Gaia services are slow, so plates already measured are kept.
python analysis/plate_blue.py -> data/plate_blue.csv"""
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

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                  # noqa: E402
from geometry import deproject                                                 # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = ROOT / "data/cutouts/depth"; CUT.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data/plate_blue.csv"
SIZE, PIX, MATCH = 10.0, 1.7, 5.0                                              # arcmin, arcsec per pixel, arcsec
OFFSETS = [(1.0, 1.0), (-1.0, -1.0)]
gaia = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag"], column_filters={"Gmag": "<21"}, row_limit=-1)


def cutout(survey, tag, q, k, c):
    f = CUT / f"{q}_{k}{tag}.fits"
    if not f.exists():
        r = requests.get(f"https://archive.stsci.edu/cgi-bin/dss_search?v={survey}&r={c.ra.deg:.5f}&d={c.dec.deg:.5f}"
                         f"&e=J2000&h={SIZE}&w={SIZE}&f=fits&c=none", timeout=300)
        r.raise_for_status(); f.write_bytes(r.content)
    return fits.getheader(f), fits.getdata(f).astype(float)


def stars(q, k, c):
    f = CUT / f"{q}_{k}_gaia21.csv"
    if not f.exists():
        t = gaia.query_region(c, width=SIZE * u.arcmin, height=SIZE * u.arcmin, catalog="I/355/gaiadr3")
        (t[0].to_pandas() if len(t) else pd.DataFrame(columns=["RA_ICRS", "DE_ICRS", "Gmag"]))[["RA_ICRS", "DE_ICRS", "Gmag"]].to_csv(f, index=False)
    return pd.read_csv(f)


def specks(h, R, g):
    """Sources, and sources with no Gaia star within 5", away from the cutout's edge."""
    w = WCS(h); im = R - median_filter(R, size=31)
    _, med, sd = sigma_clipped_stats(im)
    found = DAOStarFinder(fwhm=2.5, threshold=5 * sd)(im - med)
    if found is None:
        return 0, 0
    xc, yc = ("x_centroid", "y_centroid") if "x_centroid" in found.colnames else ("xcentroid", "ycentroid")
    x, y = np.asarray(found[xc], float), np.asarray(found[yc], float)
    ok = (x > 8) & (y > 8) & (x < R.shape[1] - 9) & (y < R.shape[0] - 9); x, y = x[ok], y[ok]
    if g.empty:
        return len(x), len(x)
    gx, gy = w.world_to_pixel(SkyCoord(g.RA_ICRS.values * u.deg, g.DE_ICRS.values * u.deg))
    d = np.hypot(x[:, None] - gx[None, :], y[:, None] - gy[None, :]).min(1) * PIX
    return len(x), int((d > MATCH).sum())


def measure(q, p):
    rec = {"plate": q, "red_sources": 0, "red_specks": 0, "blue_sources": 0, "blue_specks": 0, "blue_plate": None, "blue_date": None, "patches": 0}
    for k, (dx, dy) in enumerate(OFFSETS):
        ra, dec = deproject(np.array([dx]), np.array([dy]), p.ra_c, p.dec_c); c = SkyCoord(ra[0] * u.deg, dec[0] * u.deg)
        try:
            hr, Rr = cutout("poss1_red", "", q, k, c)
            if hr["REGION"].strip() != q:
                continue
            hb, Rb = cutout("poss1_blue", "_blue", q, k, c)
            g = stars(q, k, c)
        except Exception:
            continue
        a, b = specks(hr, Rr, g); a2, b2 = specks(hb, Rb, g)
        rec.update(red_sources=rec["red_sources"] + a, red_specks=rec["red_specks"] + b, blue_sources=rec["blue_sources"] + a2,
                   blue_specks=rec["blue_specks"] + b2, blue_plate=hb["REGION"].strip(), blue_date=str(hb.get("DATE-OBS", ""))[:10],
                   patches=rec["patches"] + 1)
    return rec


if __name__ == "__main__":
    _, plates = sample.ours()
    known = pd.read_csv(OUT).set_index("plate") if OUT.exists() else pd.DataFrame()
    todo = [q for q in plates.index if q not in known.index]
    with ThreadPoolExecutor(8) as ex:
        new = pd.DataFrame(list(ex.map(lambda q: measure(q, plates.loc[q]), todo))).set_index("plate") if todo else pd.DataFrame()
    d = pd.concat([known, new]).loc[plates.index]
    d["red_date"] = pd.read_parquet(ROOT / "data/xe_headers.parquet").drop_duplicates("plate").set_index("plate").date.reindex(d.index).astype(str).str[:10]
    d.to_csv(OUT, index_label="plate")
    print(d.describe().round(2).to_string())
    print("blue plate taken the same archive date as the red:", float((d.blue_date == d.red_date).mean()))
