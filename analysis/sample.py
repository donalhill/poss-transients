"""The sample the tests run on: a folder written by build_sample.py, named by POSS_SAMPLE. Every test reads the same files:
  transients.parquet  one row per transient: plate, ra, dec, xi and eta (deg, the offsets from the plate centre), date,
                      w (its weight: 1, or VASCO-ML's rating)
  plates.parquet      one row per plate (index plate): ra_c, dec_c, ut (mid-exposure), exp (min), date, half, acc (the
                      acceptance on the grid of shadow.py: the expected density of transients over the face of a plate),
                      depth, and eq_row
  SAMPLE.json         how the sample was made
load() returns (cat, plates) from POSS_SAMPLE; nights() gives Palomar's observing nights from its log."""
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load():
    d = Path(os.environ["POSS_SAMPLE"])
    return pd.read_parquet(d / "transients.parquet"), pd.read_parquet(d / "plates.parquet")


def nights():
    """Palomar's observing nights (VizieR VI/25) from the first to the last night with a catalogue plate, by UT date
    (the log gives each night's local evening date; exposures began after midnight UT)."""
    log = pd.read_csv(ROOT / "data/palomar_log_vi25.csv")
    ut = pd.DatetimeIndex(pd.to_datetime(log.evening_date)).normalize().unique() + pd.Timedelta(days=1)
    _, plates = load()
    return ut[(ut >= plates.date.min()) & (ut <= plates.date.max())].sort_values()
