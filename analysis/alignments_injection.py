"""Closure test for analysis/alignments.py: inject synthetic aligned groups (3-5 points, 1" perpendicular scatter,
random orientation, 2-10' long) on random plates and measure the recovered h < 3" excess. Baseline: linear fit to
h in 4-40" extrapolated to 0-3"."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import alignments as al  # runs the base analysis once (fast)

rng = np.random.default_rng(7)
edges = np.arange(0, 91, 1.0)
centres = edges[:-1] + 0.5


def excess(h):
    n, _ = np.histogram(h, edges)
    fit = (centres >= 4) & (centres <= 40)
    c = np.polyfit(centres[fit], n[fit], 1)
    exp = np.polyval(c, centres[:3]).sum()
    return n[:3].sum() - exp, np.sqrt(exp)


def run(cat):
    hs = [al.triples(g.x.values, g.y.values)[0] for _, g in cat.groupby("plate")]
    return np.concatenate(hs)


base_ex, base_sd = excess(np.load(al.ROOT / "data/alignment_triples.npz")["h"])
print(f"data: h<3\" excess vs linear baseline = {base_ex:+.0f} (Poisson sd {base_sd:.0f})")
plates = al.cat.plate.unique()
for ngroups in (25, 50, 100, 200):
    rows = []
    for _ in range(ngroups):
        g = al.cat[al.cat.plate == rng.choice(plates)]
        r = rng.integers(3, 6); L = rng.uniform(120, 600); th = rng.uniform(0, np.pi)
        x0, y0 = g.x.sample(1, random_state=int(rng.integers(1e9))).iloc[0], g.y.sample(1, random_state=int(rng.integers(1e9))).iloc[0]
        s = np.sort(rng.uniform(0, L, r))
        dx, dy = np.cos(th), np.sin(th)
        perp = rng.normal(0, 1.0, r)
        rows.append(pd.DataFrame({"plate": g.plate.iloc[0], "x": x0 + s * dx - perp * dy, "y": y0 + s * dy + perp * dx}))
    inj = pd.concat([al.cat[["plate", "x", "y"]]] + rows)
    ex, sd = excess(run(inj))
    print(f"injected {ngroups:3d} aligned groups: h<3\" excess {ex:+.0f}  (shift vs data {ex - base_ex:+.0f}, ~{(ex - base_ex)/sd:.1f} sd)")
