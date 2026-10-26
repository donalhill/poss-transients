"""The shadow's geometry on the plates it reaches, for the page: at geostationary height, where the Earth hides the Sun's disc wholly (umbra) and not at all (penumbra); on each plate the
shadow reaches (shadow_fit.json: partly and wholly covered), how far the shadow seen from Palomar sits from where the
centre of the Earth would place it (at mid-exposure), and how far its edge moves over the exposure; how many plates'
worth of the test's area it covers, and a phrase naming the plates by how much of each it covers; of the plates it
covers by a quarter or more, the richest and the poorest in transients, at the sample's plate times.
POSS_SAMPLE=DIR python analysis/shadow_geometry.py -> DIR/shadow_geometry.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
F = json.loads((D / "shadow_fit.json").read_text())
_, pl = sample.load()
share = {p["plate"]: float(p["shadow"]) for p in F["plates"]}
# the plates left out (no usable cell sunlit throughout): their share of the test's area over the exposure, as the fit measures it
WHOLE, THIN = F["whole"], F["thin"]
share.update({q: float(SH.cell_dark(pl.loc[q], SH.MAIN).mean()) for q in list(WHOLE) + list(THIN)})
edge = [float(np.degrees(np.arcsin(SH.R_EARTH / SH.R_GEO) + s * np.arcsin(SH.R_SUN / 1.496e8))) for s in (-1, 1)]
shift, sweep = [], []
for q in share:
    p = pl.loc[q]; cg, ct = (np.array(SH.shadow_ring(p.ut, f, step=30)) for f in ("geo", "topo"))
    shift.append(float(SkyCoord(*cg * u.deg).separation(SkyCoord(*ct * u.deg)).deg.mean()))
    e = [np.array(SH.shadow_ring(np.datetime64(p.ut) + np.timedelta64(int(m * 30e9), "ns"), "topo", step=30)) for m in (-p.exp, p.exp)]
    sweep.append(float(SkyCoord(*e[0] * u.deg).separation(SkyCoord(*e[1] * u.deg)).deg.mean()))


def words(s):
    """How much of a plate, in words."""
    return "all" if s >= 0.9995 else "nearly all" if s >= 0.9 else "most" if s >= 0.6 else "about half" if s >= 0.35 else "about a quarter" if s >= 0.18 else None


groups, slivers, none = {}, 0, 0
for q, s in sorted(share.items(), key=lambda kv: -kv[1]):
    w = words(s)
    if w:
        groups.setdefault(w, []).append(q)
    elif s >= 0.005:                                                        # a sliver
        slivers += 1
    else:                                                                   # reached only outside the test's area
        none += 1
parts = []
for w in ("all", "nearly all", "most", "about half", "about a quarter"):
    if w in groups:
        names = groups[w]; lst = names[0] if len(names) == 1 else ", ".join(names[:-1]) + ("," if len(names) > 2 else "") + " and " + names[-1]
        parts.append(f"{w} of {('plates ' if len(names) > 1 else 'plate ') if not parts else ''}{lst}")
if slivers:
    parts.append(f"slivers of {['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'][slivers] if slivers < 10 else slivers} other{'s' if slivers > 1 else ''}")
W = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]
if none:
    parts.append(f"none of the usable area of {W[none] if none < 10 else none} other{'s' if none > 1 else ''}")
phrase = parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + ("," if len(parts) > 2 else "") + " and " + parts[-1]
# the plates the shadow partly covers (a quarter or more of the test's area): the richest and the poorest in transients
part = [p for p in F["plates"] if p["shadow"] >= 0.18]
rich, poor = max(part, key=lambda p: p["n"]), min(part, key=lambda p: p["n"])
out = {"umbra_deg": round(edge[0], 2), "penumbra_deg": round(edge[1], 2), "plates": len(share),
       "shift_lo": round(min(shift), 1), "shift_hi": round(max(shift), 1), "sweep_lo": round(min(sweep), 1), "sweep_hi": round(max(sweep), 1),
       "plates_worth": round(sum(share.values()), 2), "coverage": phrase,
       "whole": len(WHOLE), "thin": len(THIN), "left_out": len(WHOLE) + len(THIN), "partial": len(F["plates"]),
       "rich": {"plate": rich["plate"], "n": rich["n"], "share": rich["shadow"], "in_shadow": rich["Ow"]},
       "poor": {"plate": poor["plate"], "n": poor["n"], "share": poor["shadow"], "in_shadow": poor["Ow"]}, "shares": {q: round(s, 3) for q, s in share.items()}}
(D / "shadow_geometry.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
