"""Compare our own transient candidates on one plate with VASCO's catalogue, using Solano et al. (2022)'s cuts.

Solano et al. (2022), Sect. 2: SExtractor at 5 sigma keeping SNR_WIN > 30 and FLAGS = 0; FWHM between 2 and 7 pixels;
elongation < 1.3; no Gaia EDR3 or Pan-STARRS DR2 source within 5"; sources near bright USNO-B stars removed
(R < 12.4, or R < -0.09 d + 15.3 with d in arcsec); high proper motion stars removed. Our equivalents: photutils
detection (own_pipeline.py), peak significance in place of SNR_WIN (threshold calibrated on VASCO's transients),
Gaia G in place of USNO-B R for the bright-star rule, Gaia moved to the plate epoch, Pan-STARRS objects with 3 or more
detections. AllWISE, which Solano et al. applied at a later stage, is reported separately.
python analysis/compare_vasco.py XE320"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
s = pd.read_parquet(ROOT / f"data/scans/own_{plate}.parquet")
v = pd.read_parquet(ROOT / "data/transients_with_plates.parquet"); v = v[v.plate == plate]
pix = 1.7   # arcsec per scan pixel


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


def nearest(cat_ra, cat_dec, ra, dec):
    d, i = cKDTree(xyz(cat_ra, cat_dec)).query(xyz(ra, dec))
    return np.degrees(2 * np.arcsin(d / 2)) * 3600, i


for key in ("ps1", "wise"):
    if not (ROOT / f"data/scans/{key}_{plate}.parquet").exists():
        s[f"{key}_sep"] = np.inf; continue
    # tree on our detections, catalogue streamed through it in chunks: dense fields hold tens of millions of sources
    c = pd.read_parquet(ROOT / f"data/scans/{key}_{plate}.parquet", columns=["ra", "dec"])
    tree, best = cKDTree(xyz(s.ra.values, s.dec.values)), np.full(len(s), np.inf)
    for k in range(0, len(c), 2_000_000):
        q = c.iloc[k:k + 2_000_000]
        d, i = tree.query(xyz(q.ra.values, q.dec.values), k=4, distance_upper_bound=np.radians(30 / 3600))
        ok = np.isfinite(d)
        np.minimum.at(best, i[ok], np.degrees(2 * np.arcsin(d[ok] / 2)) * 3600)
    s[f"{key}_sep"] = best                                 # beyond 30" reads as infinite
s["snr"] = s.max_value / s.noise

# bright-star rule (Solano et al.): reject if a star with mag < 12.4 is within 90', or mag < 15.3 - 0.09 d
g = pd.read_parquet(ROOT / f"data/scans/gaia_{plate}.parquet")
gb = g[g.phot_g_mean_mag < 15.3]
tree = cKDTree(xyz(gb.ra.values, gb.dec.values))
near = tree.query_ball_point(xyz(s.ra.values, s.dec.values), np.radians(170 / 3600))   # d < 170" covers the rule
spike = np.zeros(len(s), bool)
for k, idx in enumerate(near):
    if idx:
        d = np.degrees(np.arccos(np.clip(xyz(gb.ra.values[idx], gb.dec.values[idx]) @ xyz([s.ra.values[k]], [s.dec.values[k]])[0], -1, 1))) * 3600
        spike[k] = np.any(gb.phot_g_mean_mag.values[idx] < 15.3 - 0.09 * d)
s["near_bright"] = spike

shape = (s.elongation < 1.3) & s.fwhm_arcsec.between(2 * pix, 7 * pix)
nocat = (s.gaia_sep > 5) & (s.ps1_sep > 5)
base = shape & nocat & ~s.near_bright
sepv, iv = nearest(s.ra.values, s.dec.values, v.ra.values, v.dec.values)
vrec = s.iloc[iv[sepv < 3]]
print(f"{len(s):,} detections. VASCO's {len(v)} transients: detected {np.mean(sepv < 3):.0%}; "
      f"passing shape {shape.values[iv[sepv < 3]].mean():.0%}, no Gaia/PS1 {nocat.values[iv[sepv < 3]].mean():.0%}, "
      f"not near a bright star {(~s.near_bright).values[iv[sepv < 3]].mean():.0%}, all three {base.values[iv[sepv < 3]].mean():.0%}")
print(f"VASCO transients' peak significance: 2nd pct {vrec.snr.quantile(.02):.0f}, 10th {vrec.snr.quantile(.1):.0f}, median {vrec.snr.median():.0f} sigma")
print(" peak cut | our candidates (shape, no Gaia/PS1, not near bright star) | of which no WISE | VASCO recovered")
for cut in (5, 10, 15, 20, 25, 30):
    c = base & (s.snr >= cut)
    rec = c.values[iv[sepv < 3]].sum()
    print(f"  {cut:4d}    | {c.sum():8,d}{'':48s}| {(c & (s.wise_sep > 5)).sum():8,d}         | {rec}/{len(v)}")
s.to_parquet(ROOT / f"data/scans/own_{plate}.parquet")
