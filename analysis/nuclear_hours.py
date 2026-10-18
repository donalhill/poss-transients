"""Hour-resolved nuclear test association: each observing night (corrected UT dates; plate clock ~ mid-exposure)
is classed by time relative to the nearest US shot (DOE/NV-209 UT times):
    pre    : plates exposed 0-24 h BEFORE a shot (a causal effect is impossible here: placebo)
    post1  : 0-24 h after a shot
    post2  : 24-48 h after a shot
    far    : no shot within 48 h either side (baseline)
Ratio = mean transients per night in class / mean in 'far'. Null: shift all shot times by whole lunar months (+/-1..36)."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

shots = pd.read_csv(ROOT / "data/us_shot_times.csv")
shots["t"] = pd.to_datetime(shots.ut_datetime.str.replace("Z", ""), format="ISO8601").dt.tz_localize(None)
# consensus overrides where DOE/NV-209 is the lone outlier (see agent report): Union 18:10, Inca 21:26, Flathead/Blackfoot 06-11
fix = {("Union", "1954-04-25"): "1954-04-25 18:10", ("Inca", "1956-06-21"): "1956-06-21 21:26",
       ("Flathead", "1956-06-11"): "1956-06-11 18:26", ("Blackfoot", "1956-06-11"): "1956-06-11 18:26"}
for (name, date), t in fix.items():
    m = (shots.Shot_Name == name) & (shots.Date == date)
    shots.loc[m, "t"] = pd.Timestamp(t)
shots = shots[(shots.t >= "1949-11-01") & (shots.t <= "1957-05-31")]
shots["nevada"] = shots.Site.str.startswith("NTS")

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
strict = pd.to_datetime(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date)
cat = cat[pd.to_datetime(cat.obs_date.astype(str).str[:10]).isin(strict)]
pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
pl = pl[pl.ut_header.dt.normalize().isin(strict)]
cat = cat.join(pl[["ut"]], on="plate")
night_time = pl.groupby(pl.ut.dt.normalize()).ut.median()
counts = cat.groupby(cat.ut.dt.normalize()).size()

# the full 370-night universe: nights without catalogued plates take the block-rule date shift and a 07:00 UT plate time
g = {"__file__": str(ROOT / "analysis/nuclear.py"), "__name__": "nuc"}
src = (ROOT / "analysis/nuclear.py").read_text().split("rows = []")[0]
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    exec(src, g)
nights = g["universes"]["corrected"]
t_night = pd.Series([night_time.get(n, n + pd.Timedelta(hours=7)) for n in nights], index=nights)
vals = counts.reindex(nights, fill_value=0).astype(float)


def classify(shot_times):
    st = np.sort(shot_times.values.astype("datetime64[m]"))
    tn = t_night.values.astype("datetime64[m]")
    i = np.searchsorted(st, tn)
    nxt = np.where(i < len(st), (st[np.minimum(i, len(st) - 1)] - tn) / np.timedelta64(1, "h"), np.inf)
    prv = np.where(i > 0, (tn - st[np.maximum(i - 1, 0)]) / np.timedelta64(1, "h"), np.inf)
    cls = np.full(len(tn), "far", dtype=object)
    cls[(prv > 24) & (prv <= 48)] = "post2"
    cls[prv <= 24] = "post1"
    cls[(nxt <= 24) & (prv > 24)] = "pre"  # pre-shot only if not also within 24 h after another shot
    return cls


def ratios(shot_times):
    cls = classify(shot_times)
    base = vals[cls == "far"].mean()
    return {c: (vals[cls == c].mean() / base if (cls == c).sum() else np.nan, int((cls == c).sum())) for c in ("pre", "post1", "post2")}


shifts = [int(round(29.53 * k)) for k in range(-36, 37) if k != 0]
rows = []
for label, sel in (("all US", shots.nevada | ~shots.nevada), ("Nevada", shots.nevada)):
    st = shots.loc[sel, "t"]
    real = ratios(st)
    null = [ratios(st + pd.Timedelta(days=s)) for s in shifts]
    for c in ("pre", "post1", "post2"):
        r, n = real[c]
        nr = np.array([x[c][0] for x in null]); nr = nr[np.isfinite(nr)]
        rows.append((label, c, n, r, np.mean(nr >= r), np.median(nr)))
out = pd.DataFrame(rows, columns=["shots", "class", "n_nights", "ratio", "p_shift", "null_median"])
out.to_csv(ROOT / "data/nuclear_hours.csv", index=False)
print(f"nights {len(nights)}; shots in window {len(shots)} ({int(shots.nevada.sum())} Nevada); "
      f"Nevada shot hour UT median {shots[shots.nevada].t.dt.hour.median():.0f}, Pacific {shots[~shots.nevada].t.dt.hour.median():.0f}")
print(out.round(3).to_string(index=False))
