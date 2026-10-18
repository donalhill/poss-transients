"""Global look-elsewhere for the nuclear-test association.

Family of analyses: single lags -7..+7, symmetric windows +/-1..3, after-only [0, w] and before-only [-w, 0] for w = 1..3,
and some deliberately arbitrary windows; test sets all US and Nevada only. 2 x 28 = 56 members.
Statistic per member: ratio of mean transients per night in window to outside (corrected UT nights, the 370-night universe).
Null ensemble: the real test calendar slid by every whole number of days with 30 <= |s| <= 1100 (keeps test-series and
observing-run clustering, breaks the link). Local p = rank of the real ratio within the ensemble; global p = how often
the ensemble's best (minimum) local p beats the real best. Soviet and UK test dates are shown as 'irrelevant calendar' checks."""
import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
g = {"__file__": str(ROOT / "analysis/nuclear.py"), "__name__": "nuc"}
with contextlib.redirect_stdout(io.StringIO()):
    exec((ROOT / "analysis/nuclear.py").read_text().split("rows = []")[0], g)
nights = g["universes"]["corrected"]
values = g["cat"].groupby("cor_date").size().reindex(nights, fill_value=0).values.astype(float)
nt = g["nt"]
day0 = pd.Timestamp("1949-01-01")
N = ((nights - day0).days).values


def dates(mask):
    return np.unique(((pd.to_datetime(nt[mask].Date) - day0).dt.days).values)


sets = {"all US": dates(nt.Country.str.upper() == "USA"),
        "Nevada": dates((nt.Country.str.upper() == "USA") & nt.Site.str.startswith("NTS"))}
others = {"USSR": dates(nt.Country.str.upper() == "USSR"), "UK": dates(nt.Country.str.upper() == "UK")}

windows = {f"lag {k:+d}": (k,) for k in range(-7, 8)}
windows |= {f"+/-{w}": tuple(range(-w, w + 1)) for w in (1, 2, 3)}
windows |= {f"after 0..{w}": tuple(range(0, w + 1)) for w in (1, 2, 3)}
windows |= {f"before -{w}..0": tuple(range(-w, 1)) for w in (1, 2, 3)}
windows |= {"+4..+6": (4, 5, 6), "+7..+10": (7, 8, 9, 10), "-10..-7": (-10, -9, -8, -7), "-6..-4": (-6, -5, -4)}
members = [(s, w) for s in sets for w in windows]


def ratios(test_days, shift):
    out = np.empty(len(windows))
    for i, offs in enumerate(windows.values()):
        target = (test_days[:, None] + shift + np.array(offs)[None, :]).ravel()
        lab = np.isin(N, target)
        out[i] = values[lab].mean() / values[~lab].mean() if 2 <= lab.sum() < len(N) - 2 else np.nan
    return out


shifts = np.array([s for s in range(-1100, 1101) if abs(s) >= 30])
real = np.concatenate([ratios(sets[s], 0) for s in sets])
ens = np.array([np.concatenate([ratios(sets[s], sh) for s in sets]) for sh in shifts])  # (n_shift, 56)

def local_p(r, col):
    c = col[np.isfinite(col)]
    return (np.sum(c >= r) + 1) / (len(c) + 1) if np.isfinite(r) else np.nan

p_real = np.array([local_p(real[m], ens[:, m]) for m in range(len(members))])
# each ensemble member's local p's against the rest of the ensemble (rank-based), then its minimum
ranks = np.empty_like(ens)
for m in range(ens.shape[1]):
    col = ens[:, m]
    order = np.argsort(np.argsort(-np.nan_to_num(col, nan=-np.inf)))
    ranks[:, m] = (order + 1) / (len(col) + 1)
min_ens = np.nanmin(ranks, axis=1)
min_real = np.nanmin(p_real)
global_p = (np.sum(min_ens <= min_real) + 1) / (len(min_ens) + 1)

res = pd.DataFrame([(s, w, real[i], p_real[i]) for i, (s, w) in enumerate(members)],
                   columns=["test set", "window", "ratio", "local p"]).sort_values("local p")
res.to_csv(ROOT / "data/nuclear_lee.csv", index=False)
print(f"{len(members)} analyses x {len(shifts)} slid calendars")
print(f"local p < 0.05: {int((res['local p'] < 0.05).sum())} of {len(res)}   (expected by chance ~{0.05*len(res):.1f})")
print(res.head(12).round(3).to_string(index=False))
print(f"\nbest local p = {min_real:.4f} ({res.iloc[0]['test set']}, {res.iloc[0]['window']})")
print(f"GLOBAL p (best-of-56 vs slid calendars) = {global_p:.3f}")

for name, d in others.items():
    r = ratios(d, 0)
    ens_o = np.array([ratios(d, sh) for sh in shifts[::5]])
    p = np.array([local_p(r[i], ens_o[:, i]) for i in range(len(windows))])
    k = int(np.nanargmin(p))
    print(f"{name} tests (irrelevant calendar): best local p = {p[k]:.3f} ({list(windows)[k]}, ratio {r[k]:.2f}); "
          f"{int(np.sum(p < 0.05))} of {len(windows)} windows < 0.05")

# breadth statistic: how many of the 56 analyses reach local p < 0.05, real vs slid calendars
n_sig_ens = np.sum(ranks < 0.05, axis=1)
n_sig_real = int(np.sum(p_real < 0.05))
print(f"\nbreadth: real calendar has {n_sig_real} of 56 analyses at local p < 0.05; "
      f"slid calendars with at least as many: {np.mean(n_sig_ens >= n_sig_real):.3f} "
      f"(median {np.median(n_sig_ens):.0f}, 90th pct {np.percentile(n_sig_ens, 90):.0f})")

# same, with the March 1955 (Teapot) observing run removed
keep = ~((nights >= "1955-02-15") & (nights <= "1955-05-31"))
N, values = N[keep], values[keep]
real2 = np.concatenate([ratios(sets[s], 0) for s in sets])
ens2 = np.array([np.concatenate([ratios(sets[s], sh) for s in sets]) for sh in shifts])
p2 = np.array([local_p(real2[m], ens2[:, m]) for m in range(len(members))])
print(f"without Feb-May 1955 (Teapot): {int(np.sum(p2 < 0.05))} of 56 at local p < 0.05; best local p {np.nanmin(p2):.3f} "
      f"({members[int(np.nanargmin(p2))]})")
