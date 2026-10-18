"""Date each plate from asteroid trails: try header date -1/0/+1 day, keep offsets with clean trails.
Clean = trail duration within 0.6-1.5 x exposure, < 1" off track, elongation > 2. Appends to data/date_check.csv."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from asteroid_clock import ROOT, plate_header, candidates, measure

OUT = ROOT / "data/date_check.csv"
done = set(pd.read_csv(OUT).plate) if OUT.exists() else set()
for plate in sys.argv[1:]:
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
    if offs.empty:
        best, n, ast, mid = np.nan, 0, 0, np.nan
    else:
        c = offs.groupby("dd").num.nunique()
        best = int(c.idxmax()); n = int(c.max()); ast = offs.num.nunique()
        mid = offs[offs.dd == best].mid_minus_hdr.median()
    rec = dict(plate=plate, date_obs_raw=p0["date_obs_raw"], exp=p0["exp"], n_asteroids_tried=len(cands),
               n_clean_asteroids=ast, best_day_offset=best, n_at_best=n,
               offsets_seen=";".join(f"{d}:{m:+.0f}" for d, m in zip(offs.dd, offs.mid_minus_hdr)), mid_minus_hdr_min=mid)
    print(rec, flush=True)
    pd.DataFrame([rec]).to_csv(OUT, mode="a", header=not OUT.exists(), index=False)
