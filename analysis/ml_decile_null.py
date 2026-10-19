"""Null test for the decile analyses of Bruehl, Doherty, Streblyanska & Villarroel (2026).

The paper splits transients into tenths by a machine-learning score and tests each tenth one transient at a time:
shadow deficit per decile against 0.692% (binomial), top decile against the rest (two-proportion z), and the share of
top-decile transients inside the +/-1 day US nuclear window against the rest (two-proportion z). The score is built partly
from plate-level properties, so ranking by it partly ranks plates.

Here the score is replaced by a meaningless random number drawn once per plate. Everything else is the paper's own code
path: its shadow geometry (geocentric anti-sun, 8.5 deg, 08:00 UT, archive dates), its strict list of 370 dates, its
test list and window, its decile construction and its test statistics. If such random plate scores often reproduce the
published decile results, those results do not tell us that the score picks out real transients."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data/external/MLscripts"))   # VASCO's own code (fetch_external.py)
import nuclear_permutation_deciles as nuc   # noqa: E402
import shadow_deficit_deciles as shd        # noqa: E402

NREP = 1000
rng = np.random.default_rng(2026)

cat = pd.read_csv(ROOT / "data/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv")
tp = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
assert len(tp) == len(cat) and np.allclose(tp.ra.values, cat.ra.values) and np.allclose(tp.dec.values, cat.dec.values)
cat["plate"] = tp.plate.values

cat["obs_dt"] = shd.parse_obs_dates(cat)
cat["in_shadow"] = (shd.compute_shadow_seps(cat.ra.values.astype(float), cat.dec.values.astype(float), cat.obs_dt.values)
                    < shd.PENUMBRA_DEG).astype(int)
strict = set(pd.to_datetime(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv").obs_date).dt.strftime("%Y-%m-%d"))
cat["obs_date_str"] = nuc.parse_obs_dates_to_str(cat)
cat = cat[cat.obs_date_str.isin(strict)].reset_index(drop=True)
tests = nuc.load_nuclear_test_dates(ROOT / "data/external/nuclear_tests_1949_1957.csv", nuc.COUNTRY_FILTER)
window = nuc.expand_window(tests, nuc.WINDOW_DAYS) & strict
cat["in_window"] = cat.obs_date_str.isin(window).astype(int)
N, E = len(cat), shd.EXPECTED_SHADOW_FRAC
print(f"{N:,} transients on {cat.plate.nunique()} plates; in shadow {cat.in_shadow.sum()} ({100 * cat.in_shadow.mean():.3f}%, "
      f"expected {100 * E:.3f}%); in nuclear window {cat.in_window.sum()} ({100 * cat.in_window.mean():.2f}%)")

plates, pidx = np.unique(cat.plate.values, return_inverse=True)
rowid = np.arange(N)
sh, wi = cat.in_shadow.values, cat.in_window.values


def one_draw(u, per_row=False):
    order = np.lexsort((rowid, u if per_row else u[pidx]))   # rank(method="first") on the score
    dec = np.empty(N, int); dec[order] = (np.arange(N) * 10) // N
    n = np.bincount(dec, minlength=10); s = np.bincount(dec, sh, 10); w = np.bincount(dec, wi, 10)
    p_bin = stats.binom.cdf(s, n, E)                      # one-sided, fewer in shadow than expected
    z_s, p_s, *_ = shd.two_prop_z(int(s[9]), int(n[9]), int(s[:9].sum()), int(n[:9].sum()))
    z_w, p_w, *_ = nuc.two_prop_z(int(w[9]), int(n[9]), int(w[:9].sum()), int(n[:9].sum()))
    return dict(d9_deficit=100 * (s[9] / n[9] - E) / E, d9_shadow_pct=100 * s[9] / n[9],
                n_sig_deciles=int((p_bin < 0.005).sum()), z_shadow=z_s, p_shadow=p_s, z_window=z_w, p_window=p_w,
                d9_window_pct=100 * w[9] / n[9])


def summarise(R, label):
    print(f"\n{label}")
    print(f"  top-decile shadow deficit: median {R.d9_deficit.median():+.0f}%, 5-95% {R.d9_deficit.quantile(.05):+.0f}% to {R.d9_deficit.quantile(.95):+.0f}%")
    print(f"  top-decile deficit >= 55% (paper: 55%):                      {np.mean(R.d9_deficit <= -55):.3f}")
    print(f"  top decile below the rest at p <= 0.003 (paper: 0.003):      {np.mean((R.z_shadow < 0) & (R.p_shadow <= 0.003)):.3f}")
    print(f"  deciles individually significant (p < 0.005): median {R.n_sig_deciles.median():.0f}, range {R.n_sig_deciles.min()}-{R.n_sig_deciles.max()}")
    print(f"  top decile more often in nuclear window, p < 0.0001 (paper): {np.mean((R.z_window > 0) & (R.p_window < 1e-4)):.3f}")
    print(f"  top decile less often in nuclear window, p < 0.0001:         {np.mean((R.z_window < 0) & (R.p_window < 1e-4)):.3f}")


# score = plate part + transient part; "share" is the fraction of the score's variance that is shared within a plate
out = []
for share in (0.0, 0.1, 0.3, 1.0):
    R = pd.DataFrame([one_draw(np.sqrt(share) * rng.standard_normal(len(plates))[pidx] + np.sqrt(1 - share) * rng.standard_normal(N), per_row=True)
                      for _ in range(NREP)]).assign(plate_share=share)
    summarise(R, f"{NREP} random scores, {share:.0%} of their variance shared within each plate")
    out.append(R)
pd.concat(out).to_csv(ROOT / "data/ml_decile_null.csv", index=False)
