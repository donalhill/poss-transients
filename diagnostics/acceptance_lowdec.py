"""A check of the edge map: the shuffled plate times cross mostly low-declination plates (median +12 deg,
against +30 for the survey), whose transients are spread over the plate differently from the high-declination plates
that dominate the hosts' average map. This copies a sample and replaces its acceptance with the average map of the host
plates (those the real shadow never reaches, build_sample.py's rule) below dec 25 only, smoothed by the same rule.
usage: python diagnostics/acceptance_lowdec.py SAMPLE_DIR OUT_DIR"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import shadow as SH                                                          # noqa: E402
from shadow import HALF, smooth_map                                          # noqa: E402

src, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
for f in ("transients.parquet", "SAMPLE.json"):
    shutil.copy(src / f, out / f)
cat = pd.read_parquet(src / "transients.parquet"); pl = pd.read_parquet(src / "plates.parquet")
cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)]
wide = SH.wide("topo"); crossed = {q for q, p in pl.iterrows() if SH.cell_dark(p, wide).any()}
hosts = [g for q, g in cat.groupby("plate") if q not in crossed and pl.at[q, "dec_c"] < 25]
acc, sigma, _ = smooth_map(hosts)                                             # smoothed as build_sample.py does
pl["acc"] = [acc.astype(np.float32)] * len(pl); pl.to_parquet(out / "plates.parquet")
sel = json.loads((src / "SAMPLE.json").read_text()); sel.update({"acceptance": "average map of hosts below dec 25", "hosts": len(hosts), "smoothing_cells": sigma})
(out / "SAMPLE.json").write_text(json.dumps(sel, indent=1))
print(f"{len(hosts)} host plates below dec 25 -> {out}")
