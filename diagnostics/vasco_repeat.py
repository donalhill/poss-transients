"""VASCO's shadow count (Villarroel et al. 2025, PASP, section 8) repeated, then changed one step at a time to the
pooled count E1 of the shadow test's statistics (analysis/shadow_stats.py: the same statistic, global controls).
Each step counts the transients in the central 6 x 6 deg of their plates that fall in Earth's shadow at geostationary
height, against uniform control positions on the same plates: the ratio r = (K_T / N_T) / (K_C / N_C), K and N the
counts in shadow and in all, with VASCO's binomial significance z = (K_C/N_C - K_T/N_T) / sqrt(K_C/N_C^2 + V_T/N_T^2),
V_T = K_T for counts, and for rated counts the sum of the squared ratings in shadow (each transient still independent).
  1. VASCO's construction: the plates of their listed nights (vasco/vasco_sample.vasco()), their time for each plate
     (the header's start, read as hours and minutes), the shadow seen from the centre of the Earth with a sharp edge at
     8.69 deg at that one instant, each transient counted once, against their control positions
     (data/external/XE_final_random_sources.csv, 100 per plate).
  2. The same, the shadow seen from Palomar.
  3. Every plate of the catalogue on its corrected date and asteroid-calibrated mid-exposure time (the sample's).
  4. The shadow as the tests model it (analysis/shadow.py): moving over the exposure, its edge from the share of the
     Sun's disc the Earth hides; each transient counted by its share of the exposure in shadow, and the control's
     share in shadow from the plate's area.
  5. Each transient counted by VASCO-ML's rating.
Step 5 is E1, whose shuffled plate times give the chance of a ratio as low (from permutation_null.json of the
rated and unrated samples).
POSS_SAMPLE=DIR (the rated sample) python diagnostics/vasco_repeat.py RAW_DIR -> DIR/vasco_repeat.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "analysis"), str(ROOT / "vasco")]
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402
import vasco_sample                                                          # noqa: E402
from geometry import tangent, unit                                           # noqa: E402

D, RAW = Path(os.environ["POSS_SAMPLE"]), Path(sys.argv[1])
HALF = SH.HALF
COS = np.cos(np.radians(SH.SHADOW_DEG))
GEO, TOPO = SH.Shadow("geo"), SH.MAIN


def sharp(d, t, S):
    """In the shadow's sharp cone (radius SHADOW_DEG about the anti-sun axis at the orbit) along lines of sight d at t."""
    t = np.datetime64(np.asarray(t, dtype="datetime64[ns]"))
    ax = np.asarray(SH.antisun_axis(Time(t))).reshape(3)
    P = SH.orbit_point(d, t, S)
    return (P @ ax) / np.linalg.norm(P, axis=1) > COS


def count(cat, ctrl, times, S):
    """Transients and controls in shadow, plate by plate, each at its plate's time in `times`, the sharp shadow of S."""
    KT = KC = 0
    for q, t in times.items():
        c, k = cat[cat.plate == q], ctrl[ctrl.Name == q]
        if len(c):
            KT += int(sharp(unit(c.ra.values, c.dec.values), t, S).sum())
        if len(k):
            KC += int(sharp(unit(k.RA.values, k.Dec.values), t, S).sum())
    return KT, KC


def row(name, KT, NT, KC, NC, VT=None):
    r = (KT / NT) / (KC / NC)
    z = (KC / NC - KT / NT) / np.sqrt(KC / NC ** 2 + (KT if VT is None else VT) / NT ** 2)
    return {"step": name, "in_shadow": round(float(KT), 1), "transients": round(float(NT), 1), "controls_in_shadow": round(float(KC), 1),
            "controls": int(NC), "expected": round(float(NT * KC / NC), 1), "ratio": round(float(r), 4), "binomial_z": round(float(z), 1)}


ctrl = pd.read_csv(ROOT / "data/external/XE_final_random_sources.csv")
vcat, vpl = vasco_sample.vasco()
vcat = vcat[(vcat.xi.abs() < HALF) & (vcat.eta.abs() < HALF)]
vctrl = ctrl[ctrl.Name.isin(vpl.index)]
vtimes = vctrl.groupby("Name").UTobservation.first().map(pd.Timestamp)       # VASCO's time for each plate, as their code has it
steps = []
KT, KC = count(vcat, vctrl, vtimes, GEO); steps.append(row("vasco", KT, len(vcat), KC, len(vctrl)))
KT, KC = count(vcat, vctrl, vtimes, TOPO); steps.append(row("palomar", KT, len(vcat), KC, len(vctrl)))
cat, pl = sample.load()
cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)]
cctrl = ctrl[ctrl.Name.isin(pl.index)]
KT, KC = count(cat, cctrl, pl.ut, TOPO); steps.append(row("corrected_times", KT, len(cat), KC, len(cctrl)))
# the shadow as the tests model it: each transient by its share in shadow, each control by its plate's share of area
s = np.zeros(len(cat)); f = {}
for q, p in pl.iterrows():
    if SH.near_miss(p, SH.wide("topo")):
        continue
    sub = SH.sub_dark(p, TOPO)
    if not sub.any():
        continue
    f[q] = float(sub.mean())
    i = (cat.plate == q).values
    s[i] = SH.dark(cat.xi.values[i], cat.eta.values[i], p, TOPO)
J = len(pl)                                                                 # every plate, its controls spread evenly over it
KCf = sum(f.values()) / J
steps.append(row("moving_shadow", float(s.sum()), len(cat), KCf * len(cctrl), len(cctrl)))
steps.append(row("rated", float((cat.w * s).sum()), float(cat.w.sum()), KCf * len(cctrl), len(cctrl), float((cat.w ** 2 * s).sum())))
est = {k: json.loads((d / "permutation_null.json").read_text()) for k, d in (("rated", D), ("raw", RAW))}
out = {"steps": steps, "shuffled": {k: {"ratio": v["observed"]["E1"], "p_as_low": v["p"]["E1"], "null_median": v["null"]["E1"]["median"]} for k, v in est.items()},
       "plates_vasco": int(len(vpl)), "plates_all": int(J), "shadow_plates": sorted(f)}
(D / "vasco_repeat.json").write_text(json.dumps(out, indent=1))
for r_ in steps:
    print(r_)
print(out["shuffled"])
