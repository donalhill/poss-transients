"""Signal injection: remove a known fraction phi of in-shadow transients (as if sunlit glints that vanish in shadow)
and check what each estimator recovers.

Within-plate density ratio (Mantel-Haenszel style, straddling plates only):
    1 - phi_hat = sum_p O_p (1 - f_p) / sum_p (N_p - O_p) f_p
where O_p = in-shadow count, N_p = plate total, f_p = template-weighted shadow fraction. Unbiased for any f_p in (0, 1).
Fully shadowed plates carry no within-plate information and are excluded."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from shadow import ROOT

rng = np.random.default_rng(3)


def phi_density(O, N, f):
    return 1 - np.sum(O * (1 - f)) / np.sum((N - O) * f)


def run(geom, path=None, nrep=2000, nboot=20000):
    pt = pd.read_parquet(path or ROOT / f"data/plate_shadow_{geom}.parquet")
    st = pt[(pt.f_tmpl > 0) & (pt.f_tmpl < 0.97)]
    O, N, f = st.n_in.values.astype(float), st.N.values.astype(float), st.f_tmpl.values
    obs = phi_density(O, N, f)
    bs = [phi_density(O[i], N[i], f[i]) for i in (rng.integers(0, len(st), len(st)) for _ in range(nboot))]
    lo, hi = np.percentile(bs, [5, 95])
    print(f"[{geom}] {len(st)} straddling plates: phi_hat = {obs:+.3f}  (plate bootstrap sd {np.std(bs):.3f}, 90% CI {lo:+.2f} to {hi:+.2f})")
    print("   injected phi -> recovered shift (phi_hat_inj - phi_hat_obs combined as 1-(1-obs)(1-phi)), mean +/- sd")
    for phi in (0.1, 0.3, 0.5, 0.8):
        r = []
        for _ in range(nrep):
            Oi = rng.binomial(O.astype(int), 1 - phi)
            r.append(phi_density(Oi, N - (O - Oi), f))
        expect = 1 - (1 - obs) * (1 - phi)
        print(f"     {phi:.1f}: recovered {np.mean(r):+.3f} +/- {np.std(r):.3f}   (expected if unbiased {expect:+.3f})")
    return st


if __name__ == "__main__":
    for g in sys.argv[1:] or ("geo", "topo"):
        run(g)
