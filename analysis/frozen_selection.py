"""Our transient selection, with every threshold fixed on XE320 and then applied unchanged to other plates.
Solano et al. (2022)'s cuts as reproduced in solano_cuts.py / psfex_reverse.py, with the SNR_WIN and SPREAD_MODEL thresholds
set on XE320 from where VASCO's transients stop (2% and 5% of them below). The FWHM and elongation cuts are loosened
from 2 to 3 sigma, since our measurements are noisier than SExtractor's. Plate edge: within 2.9 deg of the centre.
python analysis/frozen_selection.py XE320 --freeze     (sets the threshold, writes data/frozen_selection.json)
python analysis/frozen_selection.py XE448              (applies it)"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scan_io import centre, pixel_arcsec

ROOT = Path(__file__).resolve().parents[1]
plate = sys.argv[1].upper()
s = pd.read_parquet(ROOT / f"data/scans/spread_{plate}.parquet")
xc, yc = centre(plate)
s["r"] = np.maximum((s.x - xc).abs(), (s.y - yc).abs()) * pixel_arcsec(plate) / 3600
# remove the bias of our SPREAD_MODEL (stacking blurs the PSF slightly): the PSF's own stars should average zero
own = ((s.gaia_sep < 1.5) & ~s.blended & (1 - 1 / s.elongation < 0.3)
       & (2 * s.r50).between(s.psfex_mode_fwhm / 1.2, s.psfex_mode_fwhm * 1.2) & s.sp_psfex.notna())
s["spc_psfex"] = s.sp_psfex - s.tile.map(s[own].groupby("tile").sp_psfex.median())
k = 1.4826
dxb, dyb = s.bbox_xmax - s.bbox_xmin, s.bbox_ymax - s.bbox_ymin
stages = [
    ("inside 2.9 deg of centre", s.r < 2.9),
    ("no Gaia/Pan-STARRS within 5\", not near a bright star", (s.gaia_sep > 5) & (s.ps1_sep > 5) & (s.near_bright == 0)),
    ("not blended", ~s.blended),
    ("shape (FWHM, elongation within 3 sigma, box)", ((s.fwhm_px - s.fm).abs() < 3 * k * s.fs) & ((s.elongation - s.em).abs() < 3 * k * s.es)
     & s.fwhm_px.between(2, 7, inclusive="neither") & (s.elongation < 1.3) & ((dxb - dyb).abs() < 2) & (dxb > 1) & (dyb > 1)),
]
m = pd.Series(True, index=s.index)
rows = []
for name, c in stages:
    m &= c.fillna(False)
    rows.append((name, (m & s.vasco).sum(), (m & ~s.vasco).sum()))
cfg = ROOT / "data/frozen_selection.json"
if "--freeze" in sys.argv:
    # SNR_WIN: VASCO's catalogue stops at 30 on SExtractor's scale; ours reads higher, so take the S/N below which
    # only 2% of VASCO's transients lie. SPREAD_MODEL: keep 95% of VASCO's transients that pass everything else.
    snr = float(s.loc[m & s.vasco, "snr_win"].quantile(0.02))
    th = float(s.loc[m & s.vasco & (s.snr_win > snr), "spc_psfex"].quantile(0.05))
    cfg.write_text(json.dumps({"snr_threshold": snr, "spread_threshold": th, "set_on": plate}, indent=1))
cfg = json.loads(cfg.read_text())
snr, th = cfg["snr_threshold"], cfg["spread_threshold"]
m &= s.snr_win > snr
rows.append((f"SNR_WIN > {snr:.1f} (frozen; SExtractor's 30)", (m & s.vasco).sum(), (m & ~s.vasco).sum()))
m &= s.spc_psfex > th
rows.append((f"SPREAD_MODEL > {th:.4f} (frozen)", (m & s.vasco).sum(), (m & ~s.vasco).sum()))
v = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
nv = (v.plate == plate).sum()
print(f"{plate}: {nv} VASCO transients on the plate, {(s.vasco & (s.r < 2.9)).sum()} within 2.9 deg")
print(f"{'stage':55s} {'VASCO':>6s} {'others':>8s}")
for name, a, b in rows:
    print(f"{name:55s} {a:6d} {b:8,d}")
s.loc[m, ["ra", "dec", "x", "y", "r", "tile", "snr_win", "vasco"]].to_parquet(ROOT / f"data/scans/selected_{plate}.parquet")
nin = (s.vasco & (s.r < 2.9)).sum()
if nin:
    a, b = rows[-1][1:]
    print(f"kept {a / nin:.0%} of VASCO's within 2.9 deg; total selected / VASCO = {(a + b) / nin:.2f}")
