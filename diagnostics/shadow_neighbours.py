"""The sky around the two near-empty plates that drive the pooled shadow counts: the plates within RADIUS deg of them
that no shadow reaches, and how many transients those hold against the pointing model's expectation (shadow_stats.py),
set beside every plate of the survey. The near-empty plates are those the shadow test uses that hold under a tenth of
their expected number (shadow_fit.json's pooled_by_plate, the registry's rule).
python diagnostics/shadow_neighbours.py -> $POSS_SAMPLE/shadow_neighbours.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
import sample                                                                # noqa: E402
import shadow_stats as ST                                                    # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
RADIUS, LOW = 9.0, 0.2                                                       # deg; a plate "sparse" below this share

cat, pl = sample.load()
mu, y, _ = ST.pointing(pl, cat); r = y / mu
crossed = set(json.loads((D / "crossed.json").read_text())["crossed"])
low = sorted(p["plate"] for p in json.loads((D / "shadow_fit.json").read_text())["pooled_by_plate"]
             if not p["set_aside"] and p["n"] < 0.1 * p["expected"])
c = SkyCoord(pl.ra_c.values * u.deg, pl.dec_c.values * u.deg)
near = set()
for q in low:
    sep = c[pl.index.get_loc(q)].separation(c).deg
    near |= {p for p, s in zip(pl.index, sep) if 0 < s < RADIUS and p not in crossed}
nb = r[sorted(near)]
out = {"low_plates": low, "radius_deg": RADIUS, "neighbours": len(nb), "neighbour_median": round(float(np.median(nb)), 4),
       "neighbour_share_low": round(float(np.mean(nb < LOW)), 4), "survey_median": round(float(np.median(r)), 4),
       "survey_share_low": round(float(np.mean(r < LOW)), 4), "low": LOW}
(D / "shadow_neighbours.json").write_text(json.dumps(out, indent=1))
print(out)
