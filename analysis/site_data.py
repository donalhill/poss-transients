"""The plate-level data behind the page's sky map and shadow charts, for one sample (POSS_SAMPLE): every plate's outline on the sky (its central square),
date, and transients counted in its central square; and, for each plate the shadow partly covers (shadow_fit.json's plates),
its count in shadow and the count expected there with no effect (its count times the fit's no-effect share), the shadow's edge on the sky at mid-exposure, its share of each 0.1 deg cell over the exposure, the edge at the start and
the end (and at eleven moments between, for the chart's animation), and every transient on it, at the sample's plate times.
Also the page's illustration of the shadow test's model, on an invented plate (not data) under the real shadow of the
plate that holds most transients in the test, its edge at mid-exposure: MODEL_ART artefacts drawn from the hosts' map,
and glints spread evenly, making up MODEL_SHARE of the transients, those in shadow not recorded; each point's place in
or out of shadow, the map, and the observed and expected counts in shadow, with and without the glints (with glints, the
expected share comes from the map the plates the shadow never reaches would then show: artefacts and glints together).
POSS_SAMPLE=DIR python analysis/site_data.py -> DIR/site_data.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402
from geometry import deproject, tangent                                      # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
r1, r3 = (lambda x: round(float(x), 1)), (lambda x: round(float(x), 3))
MODEL_ART, MODEL_SHARE, MODEL_SEP = 400, 0.3, 0.13
cat, pl = sample.load()
F = json.loads((D / "shadow_fit.json").read_text())
fitp = {p["plate"]: p for p in F["plates"]}
N = cat.groupby("plate").size().reindex(pl.index).fillna(0).astype(int)
edge = np.linspace(-SH.HALF, SH.HALF, 5)
sq = np.r_[np.c_[edge, np.full(5, -SH.HALF)], np.c_[np.full(5, SH.HALF), edge][1:], np.c_[edge[::-1], np.full(5, SH.HALF)][1:],
           np.c_[np.full(5, -SH.HALF), edge[::-1]][1:-1]]
plates, close = [], []
for name, p in pl.iterrows():
    ra, dec = deproject(sq[:, 0], sq[:, 1], p.ra_c, p.dec_c)
    rec = {"id": name, "N": int(N[name]), "date": p.ut.strftime("%Y-%m-%d"),
           "poly": [[r1(a), r1(b)] for a, b in zip(ra, dec)]}
    if name in fitp:                                                         # the shadow partly covers it: its edge at mid-exposure
        f = fitp[name]; cra, cdec = SH.shadow_ring(p.ut, "topo", step=6)
        rec.update({"frac": r3(f["shadow"]), "obs": r1(f["Ow"]), "exp": r1(f["nw"] * f["weighted"]), "r": f["r"],
                    "circle": [[r1(a), r1(b)] for a, b in zip(cra, cdec)]})
        c = cat[(cat.plate == name) & (cat.xi.abs() < SH.HALF) & (cat.eta.abs() < SH.HALF)]
        edges = [tangent(*SH.shadow_ring(np.datetime64(p.ut) + np.timedelta64(int(m * 30e9), "ns"), "topo"), p.ra_c, p.dec_c) for m in (-p.exp, p.exp)]
        close.append({**f, "xi": [r3(v) for v in c.xi], "eta": [r3(v) for v in c.eta],
                      "s": [round(float(v), 2) for v in SH.dark(c.xi.values, c.eta.values, p, SH.MAIN)],
                      "cells": "".join(np.base_repr(int(round(v * 35)), 36).lower() for v in SH.cell_dark(p, SH.MAIN)),
                      "edges": [[[r3(x), r3(y)] for x, y in zip(ex, ey)] for ex, ey in edges],
                      # the edge at eleven moments from the start of the exposure to its end, for the chart's animation
                      "sweep": [[[r3(x), r3(y)] for x, y in zip(*tangent(*SH.shadow_ring(np.datetime64(p.ut) + np.timedelta64(int(f_ * p.exp * 60e9), "ns"), "topo"), p.ra_c, p.dec_c))]
                                for f_ in np.linspace(-0.5, 0.5, 11)]})
    plates.append(rec)
# the illustration: an invented plate under the real shadow of the plate holding most transients in the test, its edge at
# mid-exposure (a sharp edge, so every count is whole): artefacts drawn from the hosts' map, glints spread evenly; the
# share of a plate's transients expected in shadow is the map summed over the cells in shadow
from matplotlib.path import Path as Poly                                     # noqa: E402
rich = max(close, key=lambda q: q["n"]); rng = np.random.default_rng(7)
A = np.asarray(pl.acc.iloc[0], float); A = A / A.mean()
ng = int(round(MODEL_SHARE / (1 - MODEL_SHARE) * MODEL_ART))


def spaced(n, draw, taken):
    """n points from draw(), each at least MODEL_SEP deg from every point already taken (no two transients closer than
    the catalogue's selection would keep)."""
    out = []
    while len(out) < n:
        x, y = draw()
        if all((x - u) ** 2 + (y - w_) ** 2 >= MODEL_SEP ** 2 for u, w_ in taken + out):
            out.append((x, y))
    return out


def from_map():
    c_ = rng.choice(A.size, p=A / A.sum())
    return SH.GX[c_] + rng.uniform(-SH.CELL / 2, SH.CELL / 2), SH.GY[c_] + rng.uniform(-SH.CELL / 2, SH.CELL / 2)


arts = spaced(MODEL_ART, from_map, [])
gls = spaced(ng, lambda: (rng.uniform(-SH.HALF, SH.HALF), rng.uniform(-SH.HALF, SH.HALF)), arts)
ax, ay = np.array(arts).T; gx, gy = np.array(gls).T
ring = rich["sweep"][len(rich["sweep"]) // 2]; poly = Poly(np.array(ring))
xs, ys = np.r_[ax, gx], np.r_[ay, gy]; sh = poly.contains_points(np.c_[xs, ys])
cin = poly.contains_points(np.c_[SH.GX, SH.GY]); Fm = float(A[cin].sum() / A.sum())
# with glints, the plates the shadow never reaches hold them too: their map is the artefacts' plus an even layer
Fg = (MODEL_ART * Fm + ng * float(cin.mean())) / (MODEL_ART + ng)
kept = ~sh[MODEL_ART:]; n_end = MODEL_ART + int(kept.sum()); O_end = int(sh[:MODEL_ART].sum())
model = {"plate": rich["plate"], "share": MODEL_SHARE, "ring": ring, "F": round(Fm, 4), "F_glints": round(Fg, 4), "acc": [round(float(v), 3) for v in A],
         "points": [[r3(x), r3(y), int(i >= MODEL_ART), int(s_)] for i, (x, y, s_) in enumerate(zip(xs, ys, sh))],
         "artefacts": MODEL_ART, "glints": ng, "lost": int((~kept).sum()), "art_in": O_end,
         "obs": O_end, "exp": r1(n_end * Fg), "oe": round(O_end / (n_end * Fg), 3),
         "oe_no_glints": round(O_end / (MODEL_ART * Fm), 3)}
(D / "site_data.json").write_text(json.dumps({"plates": plates, "shadow": {"plates": close, "fit": F["fit"], "half": SH.HALF}, "model": model}))
print(f"{len(plates)} plates, {len(close)} partly in shadow -> {D / 'site_data.json'}")
