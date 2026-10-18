"""Apply Solano et al. (2022)'s cuts, one at a time, to the SExtractor-style measurements from solano_cuts.py, and count
VASCO's transients kept against other candidates, inside VASCO's footprint (central 6.2 x 6.2 deg of the plate).
python analysis/solano_compare.py XE320"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from scan_io import centre, pixel_arcsec

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
t = pd.read_parquet(ROOT / f"data/scans/sex_{plate}.parquet")
v = pd.read_parquet(ROOT / "data/transients_with_plates.parquet"); v = v[v.plate == plate]
xc, yc = centre(plate)
half = 3.1 * 3600 / pixel_arcsec(plate)                     # 3.1 deg in scan pixels
t = t[t.matched & ((t.x - xc).abs() < half) & ((t.y - yc).abs() < half)].copy()


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


d, i = cKDTree(xyz(t.ra.values, t.dec.values)).query(xyz(v.ra.values, v.dec.values))
sep = np.degrees(2 * np.arcsin(d / 2)) * 3600
t["vasco"] = False
t.iloc[i[sep < 3], t.columns.get_loc("vasco")] = True
print(f"VASCO transients matched to a 5-sigma detection: {(sep < 3).sum()} of {len(v)}")

# per-tile median and MAD-sigma of FWHM and elongation, over sources passing the S/N and flag cuts
good = (t.snr_win > 30) & ~t.blended
stats = t[good].groupby("tile").agg(fm=("fwhm_px", "median"), em=("elongation", "median"))
t = t.join(stats, on="tile")
mad = lambda col, m: (t[good][col] - t[good][m]).abs().groupby(t[good].tile).median()
t["fs"] = t.tile.map(mad("fwhm_px", "fm")); t["es"] = t.tile.map(mad("elongation", "em"))

dxb, dyb = t.bbox_xmax - t.bbox_xmin, t.bbox_ymax - t.bbox_ymin
cuts = [
    ("no Gaia or Pan-STARRS within 5\"", (t.gaia_sep > 5) & (t.ps1_sep > 5)),
    ("not near a bright star", t.near_bright == 0),
    ("SNR_WIN > 30", t.snr_win > 30),
    ("FLAGS = 0 (not blended)", ~t.blended),
    ("FWHM within 2 sigma of median", (t.fwhm_px - t.fm).abs() < 2 * t.fs),
    ("elongation within 2 sigma of median", (t.elongation - t.em).abs() < 2 * t.es),
    ("2 < FWHM < 7 px", t.fwhm_px.between(2, 7, inclusive="neither")),
    ("elongation < 1.3", t.elongation < 1.3),
    ("bounding box similar, > 1 px", ((dxb - dyb).abs() < 2) & (dxb > 1) & (dyb > 1)),
    ("SPREAD_MODEL > -0.002", t.spread_model > -0.002),
]
keep = pd.Series(True, index=t.index)
print(f"{'cut':40s} {'VASCO kept':>10s} {'others':>8s}")
for name, m in cuts:
    keep &= m.fillna(False)
    print(f"{name:40s} {(keep & t.vasco).sum():10d} {(keep & ~t.vasco).sum():8,d}")
print("\nEach cut alone, applied after all the others (VASCO kept, others):")
for name, m in cuts:
    rest = pd.Series(True, index=t.index)
    for n2, m2 in cuts:
        if n2 != name:
            rest &= m2.fillna(False)
    print(f"  without {name:38s} {(rest & t.vasco).sum():4d} {(rest & ~t.vasco).sum():7,d}")
t["solano"] = keep
t.to_parquet(ROOT / f"data/scans/sex_{plate}_cut.parquet")
