"""Within-plate vanishing fraction phi with a sky-dependent density weight (crowding selection).
Expected in-shadow fraction per plate: f = sum(template * w(|b|) * in_shadow) / sum(template * w(|b|)), with
w = exp(c |b|), c = 0 (as before), 0.0224 per degree (NB fit across plates), 0.0448 (stress test). Corrected dates."""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

os.environ["LL_DATES"] = "corrected"
sys.path.insert(0, str(Path(__file__).parent))
from shadow import ROOT, SHADOW_DEG
from shadow_test import plates, cat, deproject, separations, GX, GY, W_TMPL
from injection import phi_density

rng = np.random.default_rng(5)
for geom in ("geo", "topo"):
    pt = pd.read_parquet(ROOT / f"data/plate_shadow_{geom}_corr.parquet")
    st = pt[(pt.f_tmpl > 0) & (pt.f_tmpl < 0.97)].copy()
    cells = {}
    for name in st.index:
        p = plates.loc[name]
        ra, dec = deproject(GX, GY, p.ra_c, p.dec_c)
        ins = separations(ra, dec, np.full(len(ra), p.ut), geom) < SHADOW_DEG
        b = np.abs(SkyCoord(ra * u.deg, dec * u.deg).galactic.b.deg)
        cells[name] = (ins, b)
    print(f"[{geom}, corrected dates] {len(st)} straddling plates")
    for c in (0.0, 0.0224, 0.0448):
        f = np.array([np.sum(W_TMPL * np.exp(c * cells[n][1]) * cells[n][0]) / np.sum(W_TMPL * np.exp(c * cells[n][1])) for n in st.index])
        O, N = st.n_in.values.astype(float), st.N.values.astype(float)
        phi = phi_density(O, N, f)
        bs = [phi_density(O[i], N[i], f[i]) for i in (rng.integers(0, len(st), len(st)) for _ in range(5000))]
        print(f"   c = {c:.4f}/deg: phi = {phi:+.3f} (plate bootstrap sd {np.std(bs):.2f});  f per plate: "
              + ", ".join(f"{n} {x:.3f}" for n, x in zip(st.index, f)))
    b_c = {n: np.abs(SkyCoord(plates.loc[n].ra_c * u.deg, plates.loc[n].dec_c * u.deg).galactic.b.deg) for n in st.index}
    print("   plate-centre |b|: " + ", ".join(f"{n} {v:.0f}" for n, v in b_c.items()))
