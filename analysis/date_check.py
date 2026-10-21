"""Date each plate from asteroid trails: try header date -1/0/+1 day, keep offsets with clean trails.
Clean = trail duration within 0.6-1.5 x exposure, < 1" off track, elongation > 2, and the trail's mid-time within one
exposure of mid-exposure. Each trail's mid-time is stored relative to DATE-OBS read as hours and minutes (asteroid_clock);
the exposure starts at DATE-OBS read as hours and hundredths of an hour, 0.4 x the minutes field earlier.
A plate's offset is the one most asteroids agree on; a tie leaves it undated.

Usage: python analysis/date_check.py GROUP XE422 XE369 ...   (appends to data/date_check_summary.csv)"""
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/date_check_summary.csv"


def resolve(date_obs_raw, exp, offsets_seen):
    """(day offset, asteroids agreeing) from a plate's clean trails, "dd:mid_minus_hdr;..." as stored."""
    if not isinstance(offsets_seen, str) or not offsets_seen:
        return np.nan, 0
    late = 0.4 * int(re.match(r".*T\d\d:(\d\d)", date_obs_raw).group(1))     # header-as-minutes minus true start, min
    ok = [int(dd) for dd, m in (e.split(":") for e in offsets_seen.split(";")) if abs(float(m) + late - exp / 2) < exp]
    if not ok:
        return np.nan, 0
    c = pd.Series(ok).value_counts()
    return (np.nan, 0) if (c == c.max()).sum() > 1 else (int(c.idxmax()), int(c.max()))


if __name__ == "__main__":
    from asteroid_clock import plate_header, candidates, measure
    group, todo = sys.argv[1], sys.argv[2:]
    done = set(pd.read_csv(OUT).plate) if OUT.exists() else set()
    for plate in todo:
        if plate in done:
            continue
        try:
            p0 = plate_header(plate)
            cands = candidates(p0, nmax=12)
            cands = cands[cands.V <= 15.8].head(3)
        except Exception as ex:
            print(plate, "header/skybot failed", ex, flush=True)
            continue
        rows = []
        for _, a in cands.iterrows():
            for dd in (-1, 0, 1):
                p = dict(p0, plate=f"{plate}d{dd:+d}", t_hdr=p0["t_hdr"] + pd.Timedelta(days=dd), t_ha=p0["t_ha"] + pd.Timedelta(days=dd))
                try:
                    r, _ = measure(p, a)
                except Exception:
                    r = None
                if r and 0.6 * p0["exp"] <= r["trail_min"] <= 1.5 * p0["exp"] and r["off_track"] < 1.0 and r["elong"] > 2:
                    rows.append((dd, int(a.Number), r["mid_minus_hdr"]))
        offs = pd.DataFrame(rows, columns=["dd", "num", "mid_minus_hdr"])
        seen = ";".join(f"{d}:{m:+.0f}" for d, m in zip(offs.dd, offs.mid_minus_hdr))
        best, n = resolve(p0["date_obs_raw"], p0["exp"], seen)
        mid = offs[offs.dd == best].mid_minus_hdr.median() if n else np.nan
        rec = dict(plate=plate, date_obs_raw=p0["date_obs_raw"], exp=p0["exp"], n_asteroids_tried=len(cands),
                   n_clean_asteroids=offs.num.nunique(), best_day_offset=best, n_at_best=n,
                   offsets_seen=seen, mid_minus_hdr_min=mid, group=group)
        print(rec, flush=True)
        pd.DataFrame([rec]).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
