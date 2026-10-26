"""The numbers the page quotes on the plates' dates and times, over the plates of VASCO's catalogue. From the
date records in data/: the share of plates whose archive date is a day early (to the nearest 5%), the asteroid trails
that date plates a day after the archive and those that confirm it (date_check_summary.csv), the log's agreement with
the labels and trails, the plates that do not follow their run of plate numbers (plate_date_corrections.csv), where
each plate's date and time come from (plate_times_all.csv's time_source: the date's record, then any time of day taken
from the log, then the calibration to the sky), the plates whose header and log give times of day more than 10 minutes
apart (plate_times_nominal.csv), and the calibration's survey offset, its range over the years with their own centre,
the clean trails' scatter about their plate's median (asteroid_clock_matches.csv, as dates/plate_times.py measures it),
and the precision (plate_time_posterior.csv).
usage: python report/dates_stats.py DATA_DIR PLATES.parquet OUT.json"""
import json
import sys
from pathlib import Path

import pandas as pd

data, src, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
plates = pd.read_parquet(src).index.unique()
C = pd.read_csv(data / "plate_date_corrections.csv"); C = C[C.plate.isin(plates)]
DC = pd.read_csv(data / "date_check_summary.csv"); DC = DC[DC.plate.isin(plates)]
both = C[(C.label_offset.notna() | C.asteroid_offset.notna()) & C.log_offset.notna()]
exc = C[C.rule_offset.notna() & (C.offset_days != C.rule_offset)]; d = (exc.offset_days - exc.rule_offset).abs()
res = {"plates": int(len(C)), "early_pct": 5 * round(20 * float((C.offset_days == 1).mean())),
       "asteroid_late_plates": int((DC.best_day_offset == 1).sum()), "asteroid_late_trails": int(DC[DC.best_day_offset == 1].n_at_best.sum()),
       "asteroid_ok_plates": int((DC.best_day_offset == 0).sum()), "asteroid_ok_trails": int(DC[DC.best_day_offset == 0].n_at_best.sum()),
       "log_checked": int(len(both)), "log_agree": int((both.log_offset == both.label_offset.fillna(both.asteroid_offset)).sum()),
       "exceptions": int(len(exc)), "gross": int((d > 1).sum()), "one_day": int((d == 1).sum())}
T = pd.read_csv(data / "plate_times_all.csv", index_col=0); T = T[T.index.isin(plates)]
post = pd.read_csv(data / "plate_time_posterior.csv", index_col=0); post = post[post.index.isin(plates)]
parts = T.time_source.str.split("+")
own = T.time_source.str.contains("own_asteroids")
half = (post.tau_q84 - post.tau_q16) / 2
err = post.tau_median[post.tau_median.abs() > 5].abs()      # plates whose shared written time (header and log agreeing) is wrong by more than 5 min
# plates with no time of their own in the table are dated by their correction's source
# and timed by the survey's offset, as build_sample.py does
rest = C[~C.plate.isin(T.index)]
src_ = pd.concat([parts.str[1],
                  rest.source.replace({"Palomar log": "Palomar_log", "block rule": "block_rule", "serial neighbours": "serial_neighbours"})])
ds = {k: 0 for k in ("label", "asteroids", "Palomar_log", "serial_neighbours", "block_rule")}
ds.update(src_.value_counts().to_dict())
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "dates"))
from plate_times import YEAR_MIN                                             # noqa: E402
Nom = pd.read_csv(data / "plate_times_nominal.csv", index_col=0); Nom = Nom[Nom.index.isin(plates)]
Mt = pd.read_csv(data / "asteroid_clock_matches.csv"); Mt = Mt[Mt.plate.isin(plates)]
gm = Mt.groupby("plate").dt_min; med, cnt = gm.median(), gm.size()
many = cnt[cnt >= 3].index; dev = Mt[Mt.plate.isin(many)].dt_min - Mt[Mt.plate.isin(many)].plate.map(med)
two = cnt[cnt >= 2].index; yr = pd.to_datetime(T.ut).dt.year.reindex(two)
by = med[two].groupby(yr).agg(["median", "size"]); by = by[by["size"] >= YEAR_MIN]["median"]
res["times"] = {"plates": int(len(T) + len(rest)), "date_source": ds,
                "header_log_apart": int(Nom.log_start.notna().sum()), "trail_scatter_min": round(1.4826 * float(np.median(np.abs(dev))), 2),
                "year_offset_min": [round(float(by.max()), 2), round(float(by.min()), 2)],
                "time_from_log": int(T.time_source.str.contains("log_time").sum()),
                "calibrated_own_asteroids": int(own.sum()), "calibrated_survey_offset": int((~own).sum() + len(rest)),
                "survey_offset_min": round(float(-(post.from_prior_centre - post.tau_median).median()), 2),
                "time_halfwidth68_median_min": round(float(half.median()), 2),
                "shared_record_errors": {"plates": int(len(err)), "min_abs_min": round(float(err.min()), 1), "max_abs_min": round(float(err.max()), 1)}}
out.write_text(json.dumps(res, indent=1))
print(res)
