"""The image shapes of transients on the plates of the observing runs that hold plates near nuclear tests (see
NUCLEAR_DESIGN.md): is each transient shaped like the stars around it?

On every plate of the runs that hold both plates within a day of a US test and other plates, up to 60 transients picked
at random (seeded) get a 6' DSS cutout of the red plate. The transient and the stars on the cutout are each fitted with
an elliptical Gaussian (stars.py); the transient's width is divided by the median width of the stars of similar peak
brightness (0.4 to 2.5 times its own), and the same ratio is worked out for each of those stars against the others. A
transient is star-like when its ratio lies within the 10-90% range of the stars' ratios, wider above it, narrower below.

The DSS service is slow, so transients already measured are kept.
python analysis/plate_shapes.py -> data/plate_shapes.csv"""
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
from photutils.detection import DAOStarFinder
from scipy.ndimage import median_filter

sys.path.insert(0, str(Path(__file__).parent))
from nuclear_model import GAP, P, cat, runs, window                       # noqa: E402
from stars import fit                                                          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = ROOT / "data/cutouts/shapes"; CUT.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data/plate_shapes.csv"
SIZE, PER_PLATE = 6.0, 60


def measure(q, i, ra, dec):
    f = CUT / f"{q}_{i}.fits"
    try:
        if not f.exists():
            r = requests.get(f"https://archive.stsci.edu/cgi-bin/dss_search?v=poss1_red&r={ra:.6f}&d={dec:.6f}"
                             f"&e=J2000&h={SIZE}&w={SIZE}&f=fits&c=none", timeout=300)
            r.raise_for_status(); f.write_bytes(r.content)
        h = fits.getheader(f)
        if h["REGION"].strip() != q:
            return None
        R = fits.getdata(f).astype(float)
    except Exception:
        return None
    im = R - median_filter(R, size=31)
    _, med, sd = sigma_clipped_stats(im)
    tx, ty = WCS(h).world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    t = fit(im, float(tx), float(ty))
    if t is None:
        return None
    found = DAOStarFinder(fwhm=2.5, threshold=5 * sd)(im - med)
    if found is None:
        return None
    xc, yc = ("x_centroid", "y_centroid") if "x_centroid" in found.colnames else ("xcentroid", "ycentroid")
    sat = 0.8 * R.max() - np.median(R)
    st = []
    for x, y in zip(np.asarray(found[xc], float), np.asarray(found[yc], float)):
        if np.hypot(x - tx, y - ty) < 4 or im[int(round(y)), int(round(x))] > sat:
            continue
        s = fit(im, x, y)
        if s and s["fwhm"] < 10:
            st.append(s)
    peaks, widths = np.array([s["peak"] for s in st]), np.array([s["fwhm"] for s in st])
    near = (peaks > 0.4 * t["peak"]) & (peaks < 2.5 * t["peak"])
    if near.sum() < 5:
        return None
    ref = np.median(widths[near]); ratio = t["fwhm"] / ref
    # the same ratio for each similar star against the others
    own = [widths[k] / np.median(np.delete(widths[near], j)) for j, k in enumerate(np.flatnonzero(near))]
    lo, hi = np.percentile(own, [10, 90])
    return {"plate": q, "i": int(i), "ratio": float(ratio), "stars": int(near.sum()), "lo": float(lo), "hi": float(hi),
            "class": "star-like" if lo <= ratio <= hi else ("wider" if ratio > hi else "narrower")}


if __name__ == "__main__":
    run0 = runs(GAP); S0 = window()
    use = [r_ for r_ in np.unique(run0[S0 == 1]) if ((run0 == r_) & (S0 == 0)).any()]
    plates = P.index[np.isin(run0, use)]
    rng = np.random.default_rng(61)
    jobs = []
    for q in plates:
        g = cat[cat.plate == q]
        pick = rng.choice(len(g), min(PER_PLATE, len(g)), replace=False)
        jobs += [(q, int(k), float(g.ra.iloc[k]), float(g.dec.iloc[k])) for k in pick]
    known = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=["plate", "i"])
    done = set(zip(known.plate, known.i))
    todo = [j for j in jobs if (j[0], j[1]) not in done]
    with ThreadPoolExecutor(8) as ex:
        new = [r for r in ex.map(lambda j: measure(*j), todo) if r]
    d = pd.concat([known, pd.DataFrame(new)], ignore_index=True)
    d.to_csv(OUT, index=False)
    print(len(jobs), "transients picked;", len(d), "measured")
    print(d.groupby("class").size().to_string())
