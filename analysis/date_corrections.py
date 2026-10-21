"""Corrected UT date and mid-exposure time for every XE plate: block rule from plate serial labels, overridden by direct
evidence (labels read off the scans; asteroid dating; Palomar's plate log, VizieR VI/25). Writes data/plate_date_corrections.csv and data/plates_corrected.parquet."""
from pathlib import Path
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

from date_check import resolve

ROOT = Path(__file__).resolve().parents[1]
h = pd.read_parquet(ROOT / "data/xe_headers.parquet").drop_duplicates("plate").set_index("plate")
# plates absent from the FitsArchive listing: serial label and DATE-OBS taken from DSS cutout headers
extra = pd.read_csv(ROOT / "data/extra_labels.csv", parse_dates=["header_date"]).set_index("plate")
h = pd.concat([h, extra.rename(columns={"header_date": "date"})[["label", "date"]]])


def block_offset(label):
    if label <= 28 or 639 <= label <= 962 or 1611 <= label <= 1616:
        return 0
    if 30 <= label <= 637 or 968 <= label <= 1610 or 1617 <= label <= 1618:
        return 1
    return np.nan  # 29, 638, 963-967, 1619 on: not pinned down


t = pd.DataFrame(index=h.index)
t["label"] = h.label
t["header_date"] = h.date
t["rule_offset"] = h.label.apply(block_offset)
lab = pd.read_csv(ROOT / "data/label_readings.csv", parse_dates=["second_date"]).set_index("plate")
lab = lab[lab.second_date.notna()]
t["label_offset"] = (lab.second_date - h.date.reindex(lab.index)).dt.days.reindex(t.index)
ast = pd.read_csv(ROOT / "data/date_check_summary.csv").set_index("plate")
assert all((pd.isna(b_) and pd.isna(r_[0])) or b_ == r_[0] for b_, r_ in
           zip(ast.best_day_offset, (resolve(*v) for v in ast[["date_obs_raw", "exp", "offsets_seen"]].values))), "stored asteroid offsets differ from date_check.resolve"
t["asteroid_offset"] = ast.best_day_offset.reindex(t.index)


# Palomar's plate log: matched on serial number and on position (the log repeats some serial numbers, once for a
# northern field and once for a southern one). The log gives the double date's first (evening) date; UT is the second.
cen = {r.plate: (r.ra, r.dec) for r in pd.read_csv(ROOT / "data/xe_plate_centres.csv").itertuples()}   # from the scan headers (PLTRA*, PLTDEC*)
log = pd.read_csv(ROOT / "data/palomar_log_vi25.csv", dtype={"plate_number": str})
log["num"] = pd.to_numeric(log.plate_number, errors="coerce")
log["pos"] = list(SkyCoord(log.ra_j2000.values, log.dec_j2000.values, unit=(u.hourangle, u.deg)))
t["log_ut_date"] = pd.NaT
for name, r in t[t.label.notna()].iterrows():
    if name not in cen:
        continue
    c = SkyCoord(*cen[name], unit="deg")
    m = log[(log.num == r.label) & log.pos.apply(lambda q: q.separation(c).deg < 2.5)]
    if len(m) == 1:
        t.loc[name, "log_ut_date"] = pd.Timestamp(m.evening_date.iloc[0]) + pd.Timedelta(days=1)
t["log_offset"] = (t.log_ut_date - t.header_date).dt.days

# two plates whose header date carries the same typo as their serial neighbours, which the evidence above corrects:
# XE730 (label 935) shares 933 and 934's header year, 1953 for 1954; XE787 (label 1011) has month 11 for 03, unlike 1010,
# 1012 and 1013. Each takes its neighbours' corrected date
NEIGHBOURS = {"XE730": "1954-01-01", "XE787": "1954-03-27"}
t["neighbour_offset"] = pd.Series({k: (pd.Timestamp(v) - t.header_date[k]).days for k, v in NEIGHBOURS.items()}).reindex(t.index)

t["offset_days"] = t.label_offset.fillna(t.asteroid_offset).fillna(t.log_offset).fillna(t.neighbour_offset).fillna(t.rule_offset)
t["source"] = np.select([t.label_offset.notna(), t.asteroid_offset.notna(), t.log_offset.notna(), t.neighbour_offset.notna(), t.rule_offset.notna()],
                        ["label", "asteroids", "Palomar log", "serial neighbours", "block rule"], "unresolved")
t["corrected_date"] = t.header_date + pd.to_timedelta(t.offset_days, unit="D")

# consistency: independent evidence vs the block rule
for a, b in (("log_offset", "label_offset"), ("log_offset", "asteroid_offset"), ("log_offset", "rule_offset")):
    m = t[a].notna() & t[b].notna()
    bad = t[m & (t[a] != t[b])]
    print(f"{a} vs {b}: {int(m.sum())} plates; agree {int(m.sum() - len(bad))}")
    if len(bad):
        print(bad[["label", "header_date", a, b]].to_string())
for src in ("label_offset", "asteroid_offset"):
    m = t[src].notna() & t.rule_offset.notna()
    bad = t[m & (t[src] != t.rule_offset)]
    print(f"{src}: {int(m.sum())} plates; agree with block rule: {int(m.sum() - len(bad))}; disagree: {len(bad)}")
    if len(bad):
        print(bad[["label", "header_date", "rule_offset", src]].to_string())
both = t.label_offset.notna() & t.asteroid_offset.notna()
print(f"label vs asteroid on the same plate: {int(both.sum())} plates, agree {int((t.label_offset[both] == t.asteroid_offset[both]).sum())}")
print(t.source.value_counts().to_string())
print("offsets:", t.offset_days.value_counts(dropna=False).to_dict())
t.to_csv(ROOT / "data/plate_date_corrections.csv")

# mid-exposure time. The header's time is the exposure's start in hours and hundredths of an hour (06:82 is 06:49), as
# Palomar's log confirms (header = log start + 8 h, to within a minute); VASCO's UTobservation reads it as hours and
# minutes. Plates without a header keep VASCO's time, read the same way, and the median exposure
pl = pd.read_parquet(ROOT / "data/plates.parquet")
pl["offset_days"] = t.offset_days.reindex(pl.index)
pl["ut_header"] = pl.ut
hh = h.hh.reindex(pl.index).fillna(pd.Series(pl.ut.dt.hour, pl.index))
mm = h.mm.reindex(pl.index).fillna(pd.Series(pl.ut.dt.minute, pl.index))
exp = h.exp.reindex(pl.index).fillna(h.exp.median())
start = pl.ut.dt.normalize() + pd.to_timedelta(hh + mm / 100, unit="h")
pl["ut"] = start + pd.to_timedelta(exp / 2, unit="min") + pd.to_timedelta(pl.offset_days.fillna(0), unit="D")
pl.to_parquet(ROOT / "data/plates_corrected.parquet")
cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
used = pl[pl.index.isin(cat.plate.unique())]
print(f"\nanalysis plates {len(used)}: offset counts {used.offset_days.value_counts(dropna=False).to_dict()}")
