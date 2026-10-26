"""The two samples of VASCO's catalogue used on the page.

ours:  all plates in VASCO's catalogue and all its transients, each plate on its corrected date.
vasco: VASCO's construction, to reproduce their numbers: plates whose archive date is on their list of observing
       nights (poss_i_observing_dates_strict.csv), on archive dates.

Each returns (cat, plates): cat has one row per transient with its plate, plate-frame offsets xi and eta (deg),
`date`, the night it was taken (UT date for ours, archive date for vasco), and `w`, VASCO-ML's machine-learning rating
of how likely it is to be real (its weight in the weighted tests); plates has one row per plate with its centre, time
`ut` and `date` (ours: corrected mid-exposure time; vasco: VASCO's UTobservation, the header's start time read as hours
and minutes)."""
import zipfile
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


@cache
def ratings():
    """VASCO-ML's rating of every transient, in the catalogue's order (checked on position and date)."""
    with zipfile.ZipFile(ROOT / "data/external/scored_catalog.zip") as z, z.open("scored_catalog.csv") as f:
        sc = pd.read_csv(f, usecols=["ra", "dec", "obs_date", "prob"])
    cat = pd.read_csv(ROOT / "data/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv", usecols=["ra", "dec", "obs_date"])
    assert np.allclose(sc.ra, cat.ra) and np.allclose(sc.dec, cat.dec) and \
        (sc.obs_date.astype(str).str[:10] == cat.obs_date.astype(str).str[:10]).all(), "ratings out of step with the catalogue"
    return sc.prob.values


def _catalogue():
    cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
    cat["w"] = ratings()[cat.index.values]
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


def labelled():
    """VASCO-ML's hand labels (1 real, 0 flaw), each matched to a catalogue transient within 2" (row; -1 if none), with
    that transient's rating and its archive date."""
    from scipy.spatial import cKDTree
    lab = pd.read_csv(ROOT / "data/external/combined_training_labels.csv")
    cat = pd.read_csv(ROOT / "data/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv", usecols=["ra", "dec", "obs_date"])
    xyz = lambda ra, de: np.c_[np.cos(np.radians(de)) * np.cos(np.radians(ra)), np.cos(np.radians(de)) * np.sin(np.radians(ra)), np.sin(np.radians(de))]
    d, i = cKDTree(xyz(cat.ra.values, cat.dec.values)).query(xyz(lab.ra.values, lab.dec.values))
    lab["row"] = np.where(np.degrees(d) * 3600 < 2, i, -1)
    lab["prob"] = np.where(lab.row >= 0, ratings()[lab.row.clip(0)], np.nan)
    lab["archive_date"] = pd.to_datetime(np.where(lab.row >= 0, cat.obs_date.astype(str).str[:10].reindex(lab.row.clip(0)).values, None))
    return lab


def vasco_nights():
    """VASCO's list of observing nights (local evening dates, from Palomar's log)."""
    return pd.DatetimeIndex(pd.to_datetime(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv").obs_date)).normalize().unique()
