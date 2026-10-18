"""Nuclear test with a self-consistent date convention. The strict list (VI/25) gives each night's local evening date,
so the UT date of every night is list date + 1. Universe U = list + 1; transients (all 107,875, not pre-filtered on
header date) counted by corrected UT date; those not landing in U are reported. Compared with VASCO's construction
(header dates, list as universe). Raw counts; +/-1 window and single lags; permutation and lunar-shift nulls."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NPERM = 5000
log = pd.DatetimeIndex(pd.to_datetime(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date)).normalize()
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
corr = pd.read_csv(ROOT / "data/plate_date_corrections.csv", index_col=0)
cat["hdr"] = pd.to_datetime(cat.obs_date.astype(str).str[:10])
cat["ut"] = cat.hdr + pd.to_timedelta(cat.plate.map(corr.offset_days).fillna(0), unit="D")
nt = pd.read_csv(ROOT / "VASCO-ML/data/nuclear_tests_1949_1957.csv")
us = nt[nt.Country.str.upper() == "USA"]
T_all = pd.DatetimeIndex(pd.to_datetime(us.Date).drop_duplicates())
T_nts = pd.DatetimeIndex(pd.to_datetime(us[us.Site.str.startswith("NTS")].Date).drop_duplicates())
shifts = [int(round(29.53 * k)) for k in range(-36, 37) if k != 0]


def perm_p(v, lab, rng):
    obs = v[lab].mean() - v[~lab].mean(); n, k = len(v), lab.sum(); ge = 0
    for _ in range(NPERM):
        m = np.zeros(n, bool); m[rng.choice(n, k, replace=False)] = True
        ge += v[m].mean() - v[~m].mean() >= obs
    return (ge + 1) / (NPERM + 1)


def ratio(nights, v, T, offs):
    tgt = pd.DatetimeIndex(sorted({t + pd.Timedelta(days=k) for t in T for k in offs}))
    lab = nights.isin(tgt)
    return lab, (v[lab].mean() / v[~lab].mean() if 2 <= lab.sum() else np.nan)


def run(name, nights, counts, T):
    v = counts.reindex(nights, fill_value=0).values.astype(float)
    out = [name]
    for lbl, offs in (("win", (-1, 0, 1)), ("-1", (-1,)), ("0", (0,)), ("+1", (1,))):
        lab, r = ratio(nights, v, T, offs)
        p = perm_p(v, lab, np.random.default_rng(1))
        null = np.array([ratio(nights, v, T + pd.Timedelta(days=s), offs)[1] for s in shifts]); null = null[np.isfinite(null)]
        out.append(f"{lbl}: n={lab.sum():2d} r={r:4.2f} p={p:.3f} pL={np.mean(null >= r):.3f}")
    print(" | ".join(out))


# VASCO construction: header dates, universe = list dates
c_hdr = cat[cat.hdr.isin(log)].groupby("hdr").size()
U = log + pd.Timedelta(days=1)
c_ut = cat[cat.ut.isin(U)].groupby("ut").size()
print(f"VASCO construction: {c_hdr.sum()} transients kept; consistent UT construction: {c_ut.sum()} kept "
      f"({(~cat.ut.isin(U)).sum()} land outside list+1)")
for tname, T in (("all US", T_all), ("Nevada", T_nts)):
    run(f"header/list   {tname:7s}", log, c_hdr, T)
    run(f"UT/list+1     {tname:7s}", U, c_ut, T)

print("\nPhysically defined +/-24 h window around a ~13 UT dawn shot = UT nights {T, T+1} (post hoc, but VASCO could argue it):")
for tname, T in (("all US", T_all), ("Nevada", T_nts)):
    for name, nights, counts in (("UT/list+1", U, c_ut),):
        v = counts.reindex(nights, fill_value=0).values.astype(float)
        lab, r = ratio(nights, v, T, (0, 1))
        p = perm_p(v, lab, np.random.default_rng(1))
        null = np.array([ratio(nights, v, T + pd.Timedelta(days=s), (0, 1))[1] for s in shifts]); null = null[np.isfinite(null)]
        lab2 = lab & ~((nights >= "1955-03-15") & (nights <= "1955-03-31"))
        keep = ~((nights >= "1955-03-15") & (nights <= "1955-03-31"))
        r2 = v[lab & keep].mean() / v[~lab & keep].mean()
        p2 = perm_p(v[keep], lab[keep], np.random.default_rng(1))
        print(f"  {name} {tname}: n={lab.sum()} ratio {r:.2f} perm p {p:.3f} lunar-shift p {np.mean(null >= r):.3f}; "
              f"without Teapot Mar 1955: ratio {r2:.2f} perm p {p2:.3f}")
