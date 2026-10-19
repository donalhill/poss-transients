"""Independent date check: VASCO's poss_i_observing_dates_strict.csv (370 nights, from the Palomar POSS-I
observation list, VizieR VI/25) uses one convention for every night. Compare each plate's header date with it:
  header in list, header-1 not  -> header is the list's (evening) date      -> 'first'
  header not in list, header-1 in -> header is the day after the list date  -> 'second'
  both in list -> undecidable; neither -> other error.
If the list uses the local evening (first) date, 'first' => our offset +1 and 'second' => offset 0."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
log = set(pd.to_datetime(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv").obs_date).dt.normalize())
c = pd.read_csv(ROOT / "data/plate_date_corrections.csv", index_col=0)
h = pd.read_parquet(ROOT / "data/xe_headers.parquet").drop_duplicates("plate").set_index("plate")
t = pd.read_parquet(ROOT / "data/transients_with_plates.parquet"); t = t[t.plate != ""]
N = t.groupby("plate").size()
c["hd"] = pd.to_datetime(c.header_date).dt.normalize()
a, b = c.hd.isin(log), (c.hd - pd.Timedelta(days=1)).isin(log)
c["log_says"] = np.select([a & ~b, ~a & b, a & b], ["first(+1)", "second(0)", "ambiguous"], "neither")
c["log_offset"] = c.log_says.map({"first(+1)": 1.0, "second(0)": 0.0})
c["analysis"] = c.index.isin(N.index)
lo, hi = min(log), max(log)
c = c[(c.hd >= lo) & (c.hd <= hi + pd.Timedelta(days=1))]
print(pd.crosstab([c.source, c.offset_days], c.log_says).to_string())
d = c[c.log_offset.notna() & c.offset_days.isin([0, 1])]
for src in ("block rule", "asteroids", "label"):
    s = d[d.source == src]
    print(f"{src:10s}: decidable {len(s)}, agree with log {int((s.log_offset == s.offset_days).sum())}, disagree {int((s.log_offset != s.offset_days).sum())}")
s = d[d.log_offset != d.offset_days]
print("disagreements:"); print(s[["label", "header_date", "offset_days", "source", "log_says", "analysis"]].sort_values("label").to_string())
print("\nXE113/XE773:"); print(c.loc[c.index.intersection(["XE113", "XE773", "XE068", "XE098"]), ["label", "header_date", "offset_days", "source", "log_says"]].to_string())
c.to_csv(ROOT / "data/dates_vs_night_list.csv")
