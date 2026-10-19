"""Nuclear-test association, as in Bruehl et al. (nuclear_permutation_deciles.py), with header vs corrected plate dates.

Unit: observing night (UT date), value = transient count (public ML scores are degenerate, so raw counts).
Window: US test date +/- 1 day. Lag k: nights with obs date = test date + k (k < 0: plate taken BEFORE the test date).
One-sided permutation of night labels, 10,000 iterations (their method)."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NPERM = 10000

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
strict = pd.to_datetime(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv").obs_date)
nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
tests = pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date).dt.normalize().drop_duplicates()
corr = pd.read_csv(ROOT / "data/plate_date_corrections.csv", index_col=0)
plates = pd.read_parquet(ROOT / "data/plates.parquet")

cat["hdr_date"] = pd.to_datetime(cat.obs_date.astype(str).str[:10])
cat = cat[cat.hdr_date.isin(strict)].copy()
cat["cor_date"] = cat.hdr_date + pd.to_timedelta(cat.plate.map(corr.offset_days).fillna(0), unit="D")

# universe of nights: header version = their 370 strict dates; corrected = corrected dates of all plates on those dates
pl = plates.assign(hdr_date=plates.ut.dt.normalize())
pl["cor_date"] = pl.hdr_date + pd.to_timedelta(pl.index.map(corr.offset_days).fillna(0), unit="D")
pl = pl[pl.hdr_date.isin(strict)]
# date-block rule (from serial-label blocks) for the strict nights, incl. ~94 with no catalogued transients;
# nights with plates take the plate-level correction
def block_by_date(d):
    if d <= pd.Timestamp("1949-11-25") or pd.Timestamp("1953-01-09") <= d <= pd.Timestamp("1954-01-05") or d >= pd.Timestamp("1957-04-28"):
        return 0
    if pd.Timestamp("1949-12-21") <= d <= pd.Timestamp("1952-12-22") or pd.Timestamp("1954-01-27") <= d <= pd.Timestamp("1956-06-03"):
        return 1
    return np.nan
plate_shift = (pl.cor_date - pl.hdr_date).dt.days.groupby(pl.hdr_date).agg(lambda s: s.mode().iloc[0])
shift = pd.Series([plate_shift.get(d, block_by_date(d)) for d in strict], index=strict)
print("strict nights with unresolved shift:", int(shift.isna().sum()))
universes = {"header": pd.DatetimeIndex(strict).unique(),
             "corrected": pd.DatetimeIndex(strict + pd.to_timedelta(shift.fillna(0).values, unit="D")).unique()}


def perm_test(values, labels, rng):
    obs = values[labels].mean() - values[~labels].mean()
    n, k = len(values), labels.sum()
    ge = 0
    for _ in range(NPERM):
        m = np.zeros(n, bool); m[rng.choice(n, k, replace=False)] = True
        ge += values[m].mean() - values[~m].mean() >= obs
    return (ge + 1) / (NPERM + 1)


rows = []
for version, col in (("header", "hdr_date"), ("corrected", "cor_date")):
    nights = universes[version]
    counts = cat.groupby(col).size().reindex(nights, fill_value=0)
    values = counts.values.astype(float)
    for name, offsets in [("window +/-1", (-1, 0, 1))] + [(f"lag {k:+d}", (k,)) for k in range(-3, 4)]:
        target = pd.DatetimeIndex(sorted({t + pd.Timedelta(days=k) for t in tests for k in offsets}))
        labels = nights.isin(target)
        if labels.sum() < 2:
            continue
        p = perm_test(values, labels, np.random.default_rng(42))
        rows.append((version, name, int(labels.sum()), values[labels].mean() / values[~labels].mean(), p))
res = pd.DataFrame(rows, columns=["dates", "test", "n_nights", "ratio", "p_one_sided"])
res.to_csv(ROOT / "data/nuclear_header_vs_corrected.csv", index=False)
print(f"nights: header {len(universes['header'])}, corrected {len(universes['corrected'])}; US test dates {len(tests)}")
print(res.pivot(index="test", columns="dates", values=["n_nights", "ratio", "p_one_sided"]).round(3).to_string())


# ---- robustness: Nevada-only, and a lunar-month time-shift null (keeps test clustering and moon-phase relation) ----
def window_ratio(nights, values, test_dates, offsets=(-1, 0, 1)):
    target = pd.DatetimeIndex(sorted({t + pd.Timedelta(days=k) for t in test_dates for k in offsets}))
    lab = nights.isin(target)
    return (values[lab].mean() / values[~lab].mean() if lab.sum() >= 2 else np.nan), int(lab.sum())


us = nt[nt.Country.str.upper() == "USA"]
sets = {"all US (as coded)": tests,
        "Nevada only (as described)": pd.to_datetime(us[us.Site.str.startswith("NTS")].Date).dt.normalize().drop_duplicates()}
shifts = [int(round(29.53 * k)) for k in range(-36, 37) if k != 0]
print("\nlunar-month shift null (72 shifts, +/-1..36 months):")
for version, col in (("header", "hdr_date"), ("corrected", "cor_date")):
    nights = universes[version]
    values = cat.groupby(col).size().reindex(nights, fill_value=0).values.astype(float)
    for sname, td in sets.items():
        for label, offs in (("window +/-1", (-1, 0, 1)), ("lag -1", (-1,)), ("lag 0", (0,)), ("lag +1", (1,))):
            r0, n0 = window_ratio(nights, values, td, offs)
            null = np.array([window_ratio(nights, values, td + pd.Timedelta(days=s), offs)[0] for s in shifts])
            null = null[np.isfinite(null)]
            print(f"  {version:9s} {sname:27s} {label:11s} nights {n0:2d}  ratio {r0:5.2f}  "
                  f"P(shifted >= real) = {np.mean(null >= r0):.3f}  (null median {np.median(null):.2f}, n={len(null)})")
