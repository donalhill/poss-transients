"""Is the broad nuclear-window association a season / galactic-latitude effect?
(a) per-night mean |b| of plates vs transient counts; month-of-year profile of |b|, counts and test dates;
(b) season-preserving null: slide the test calendar by k years +/- d days (k = +/-1..5, |d| <= 20);
(c) |b|-adjusted counts: divide each night's count by the count predicted from its plates' |b| (log-linear fit)."""
import contextlib, io
from pathlib import Path
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
g = {"__file__": str(ROOT / "analysis/nuclear_lee.py"), "__name__": "lee"}
with contextlib.redirect_stdout(io.StringIO()):
    exec((ROOT / "analysis/nuclear_lee.py").read_text().split("shifts = np.array")[0], g)
nights, values, N, sets, windows, ratios_fn = g["nights"], g["values"], g["N"], g["sets"], g["windows"], g["ratios"]

pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
pl["b"] = np.abs(SkyCoord(pl.ra_c.values * u.deg, pl.dec_c.values * u.deg).galactic.b.deg)
night_b = pl.groupby(pl.ut.dt.normalize()).b.mean().reindex(nights)
ok = night_b.notna().values
print(f"nights with plate |b|: {ok.sum()} of {len(nights)}")
rho = pd.Series(values[ok]).corr(pd.Series(night_b.values[ok]), method="spearman")
print(f"Spearman(counts per night, mean |b|) = {rho:.2f}")
m = pd.DataFrame({"month": nights.month, "b": night_b.values, "n": values})
tests = pd.to_datetime(g["nt"][(g["nt"].Country.str.upper() == "USA") & g["nt"].Site.str.startswith("NTS")].Date)
prof = m.groupby("month").agg(nights=("n", "size"), mean_count=("n", "mean"), mean_b=("b", "mean"))
prof["nevada_tests"] = tests.dt.month.value_counts().reindex(prof.index, fill_value=0)
print(prof.round(1).to_string())

# (c) |b|-adjusted values: log(count+1) ~ a + c*|b| fitted on nights with plates
x, y = night_b.values[ok], np.log1p(values[ok])
c1, c0 = np.polyfit(x, y, 1)
pred = np.full(len(values), np.nan); pred[ok] = np.expm1(c0 + c1 * x)
adj = np.where(ok, values / np.maximum(pred, 1), np.nan)

def ratios(vals, test_days, shift):
    out = []
    for offs in windows.values():
        target = (test_days[:, None] + shift + np.array(offs)[None, :]).ravel()
        lab = np.isin(N, target) & np.isfinite(vals)
        rest = ~np.isin(N, target) & np.isfinite(vals)
        out.append(vals[lab].mean() / vals[rest].mean() if lab.sum() >= 2 else np.nan)
    return np.array(out)

def family(vals, shifts):
    real = np.concatenate([ratios(vals, sets[s], 0) for s in sets])
    ens = np.array([np.concatenate([ratios(vals, sets[s], sh) for s in sets]) for sh in shifts])
    p = np.array([(np.sum(ens[:, i][np.isfinite(ens[:, i])] >= real[i]) + 1) / (np.isfinite(ens[:, i]).sum() + 1) for i in range(len(real))])
    rk = np.empty_like(ens)
    for i in range(ens.shape[1]):
        col = np.nan_to_num(ens[:, i], nan=-np.inf); rk[:, i] = (np.argsort(np.argsort(-col)) + 1) / (len(col) + 1)
    nsig = int(np.sum(p < 0.05)); breadth = np.mean(np.sum(rk < 0.05, axis=1) >= nsig)
    glob = (np.sum(np.nanmin(rk, axis=1) <= np.nanmin(p)) + 1) / (len(rk) + 1)
    return nsig, breadth, np.nanmin(p), glob

any_shift = np.array([s for s in range(-1100, 1101) if abs(s) >= 30])
year_shift = np.array([365 * k + d for k in (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5) for d in range(-20, 21)])
for label, vals in (("raw counts", values), ("|b|-adjusted", adj)):
    for nl, sh in (("any-day slides", any_shift), ("season-preserving (whole years +/-20 d)", year_shift)):
        nsig, br, pmin, glob = family(vals, sh)
        print(f"{label:13s} | {nl:40s}: {nsig:2d}/56 at p<0.05, breadth p = {br:.3f}, best local p = {pmin:.3f}, global p = {glob:.3f}")
