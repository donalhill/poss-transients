"""The two samples used on the page (see README.md).

ours:  all plates in VASCO's catalogue and all its transients, each plate on its corrected date.
vasco: VASCO's construction, to reproduce their numbers: plates whose archive date is on their list of observing
       nights (poss_i_observing_dates_strict.csv), on archive dates.

Each returns (cat, plates): cat has one row per transient with its plate, plate-frame offsets xi and eta (deg) and
`date`, the night it was taken (UT date for ours, archive date for vasco); plates has one row per plate with its
centre, mid-exposure time `ut` and `date`."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _catalogue():
    cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
    return cat[cat.plate != ""].copy()


def ours():
    cat = _catalogue()
    plates = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
    plates = plates[plates.index.isin(cat.plate.unique())].copy()
    plates["date"] = plates.ut.dt.normalize()
    cat["date"] = cat.plate.map(plates.date)
    return cat, plates


def vasco():
    cat = _catalogue()
    listed = set(vasco_nights())
    cat["date"] = pd.to_datetime(cat.obs_date.astype(str).str[:10])
    cat = cat[cat.date.isin(listed)].copy()
    plates = pd.read_parquet(ROOT / "data/plates.parquet")
    plates = plates[plates.index.isin(cat.plate.unique())].copy()
    plates["date"] = plates.ut.dt.normalize()
    return cat, plates


def vasco_nights():
    """VASCO's list of observing nights (local evening dates, from Palomar's log)."""
    return pd.DatetimeIndex(pd.to_datetime(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv").obs_date)).normalize().unique()


def nights():
    """Palomar's observing nights (VizieR VI/25) from the first to the last night with a catalogue plate, by UT date
    (the log gives each night's local evening date; exposures began after midnight UT)."""
    log = pd.read_csv(ROOT / "data/palomar_log_vi25.csv")
    ut = pd.DatetimeIndex(pd.to_datetime(log.evening_date)).normalize().unique() + pd.Timedelta(days=1)
    _, plates = ours()
    return ut[(ut >= plates.date.min()) & (ut <= plates.date.max())].sort_values()
