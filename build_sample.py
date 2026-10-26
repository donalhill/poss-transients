"""A sample for the tests (read by analysis/sample.py) from VASCO's catalogue: their transients on their plates
(vasco/plates.py), each plate on its corrected date and time (dates/date_corrections.py), on the tests' grid
(analysis/shadow.py: 0.1 deg cells over the central 6 x 6 deg of each plate).
Plate times (--times calibrated, the default): every plate's mid-exposure time and exposure from the asteroid-calibrated
table (data/plate_times_all.csv, dates/plate_times.py: its own asteroid trails, else the survey's measured clock
offset). --times header: the corrected times alone, without the calibration.
The acceptance, the density of transients over the face of a plate, is the average over the plates the shadow never
reaches (each plate counted equally, transients by their weight), smoothed by a Gaussian whose width is chosen by how
well the map from nine tenths of these plates predicts where the transients lie on the other tenth (smooth_map()), the
same for every plate, and written to every plate so the fit and its null use the same map (--acceptance flat: a flat
map, as a check). With --weights rating each
transient counts by VASCO-ML's rating of how likely it is to be real; with none, once.
Writes, to OUT: transients.parquet, plates.parquet, SAMPLE.json (analysis/sample.py).
usage: python build_sample.py OUT [--weights rating|none] [--times calibrated|header] [--acceptance hosts|flat]"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / "vasco"), str(ROOT / "analysis")]
import vasco_sample                                                          # noqa: E402
import shadow as SH                                                          # noqa: E402
from shadow import HALF, smooth_map                                          # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("--weights", default="rating", choices=["rating", "none"])
ap.add_argument("--times", default="calibrated", choices=["calibrated", "header"])
ap.add_argument("--acceptance", default="hosts", choices=["hosts", "flat"])   # flat: no edge map (a check of the map)
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

cat, P = vasco_sample.ours()
P["time_source"] = "header"
if a.times == "calibrated":
    T = pd.read_csv(ROOT / "data/plate_times_all.csv", index_col=0, parse_dates=["ut"]).reindex(P.index)
    assert T.ut.notna().all(), "a plate with no calibrated time"
    P["ut"], P["exp"], P["time_source"] = T.ut, T.exp, T.time_source
if a.weights == "none":
    cat["w"] = 1.0
cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].copy()           # the tests' square

# the acceptance, from the plates the shadow never reaches (crossed.py's rule)
wide = SH.wide("topo")
crossed = {q for q, p in P.iterrows() if SH.cell_dark(p, wide).any()}
hosts = [g for q, g in cat.groupby("plate") if q not in crossed]
if a.acceptance == "hosts":
    acc, sigma, cv = smooth_map(hosts)
else:
    acc, sigma, cv = np.ones(SH.GX.size), None, {}

# depth: the share of Gaia stars with 18 < G < 20 found on the plate (vasco/plate_depth.py), the median where unmeasured
depth = pd.read_csv(ROOT / "data/plate_depth.csv").set_index("plate").completeness.reindex(P.index)
depth = depth.fillna(depth.median())
pl = pd.DataFrame({"x": P.x, "y": P.y, "z": P.z, "ut": P.ut, "ra_c": P.ra_c, "dec_c": P.dec_c, "exp": P.exp, "half": HALF,
                   "time_source": P.time_source}, index=P.index)
pl["acc"] = [acc.astype(np.float32)] * len(pl)
pl["depth"] = depth.reindex(pl.index)
pl["eq_row"] = (pl.dec_c.abs() < 3).astype(float)                          # the row centred on the equator (report/vasco_stats.py)
pl["date"] = pl.ut.dt.normalize()
pl.index.name = "plate"
cat["date"] = cat.plate.map(pl.date)
cat["row"] = cat.index.values                                                # the catalogue row (for other ratings in its order)
cat[["plate", "ra", "dec", "xi", "eta", "date", "w", "snr", "row"]].reset_index(drop=True).to_parquet(out / "transients.parquet")
pl.to_parquet(out / "plates.parquet")
(out / "SAMPLE.json").write_text(json.dumps({"weights": a.weights, "times": a.times, "grid": "shadow.py HALF=3, CELL=0.1",
                                                "acceptance": "hosts' average map" if a.acceptance == "hosts" else "flat", "hosts": len(hosts), "crossed": len(crossed),
                                                "smoothing_cells": sigma, "smoothing_cv": {str(k): round(v, 5) for k, v in cv.items()}}, indent=1))
print(f"{len(cat):,} transients on {cat.plate.nunique()} plates; {len(crossed)} plates crossed, {len(hosts)} hosts")
