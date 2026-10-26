"""Retrain VASCO-ML's rating model with VASCO's own code, on different sets of their labels, with and without its
per-night inputs, and score every transient each time.

VASCO's code (MLscripts at the commit fetch_external.py pins) and its feature cache are fetched into
data/external/vasco_ml/, checked against the files analysed. Their training reads four label files that are not in their
repository; their merged label file (combined_training_labels.csv, which records each label's file) stands in for them.
No calibrator is applied, as in their released ratings.

Label sets: all 393 labels; without the labels on nights within a day of a test on VASCO's list (archive dates, as their
construction dates them); the 250 their documentation describes (all but the steve batch); and those 250 without the
same nights. Each with the per-night inputs (plate_low_snr_frac, plate_snr_std, plate_elongation_mean, plate_n_high_snr,
computed per observing night) as VASCO compute them, and held at zero (their own fallback when SNR is missing).

Runs under VASCO-ML's pinned versions (vasco/requirements-vasco-ml.txt), not the page's; diagnostics/vasco_ratings.py
runs the rated tests on the results in the page's environment.

PYML=<python with requirements-vasco-ml.txt> make data/external/vasco_ml/ratings/doc250_no_window_without_3.npy
    -> data/external/vasco_ml/ratings/<labels>_<night>_<seed>.npy (not kept: large and remade on demand)"""
import hashlib
import sys
from pathlib import Path

import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import vasco_sample as sample                                                  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/external/vasco_ml"
ML = "418e8c89b2d3e7c5f56ab75293847ede54a5a445"
URL = f"https://raw.githubusercontent.com/dca-doherty/VASCO-ML/{ML}/"
FILES = {"MLscripts/common.py": "c26a24e132856e9d0c24e472b19a8454fdc3ca288732c72f1440de09348e289f",
         "MLscripts/pipeline_train.py": "454e6ce2f92bb2366ac8a11c18220bd56f2013a72269b47993dd8b266e940b75",
         "MLscripts/pipeline_score.py": "7678f7fe843802fc6987c1d7bf1b6c6b4039b347201dc5d20ba1ca262ae64f95",
         "MLscripts/config_red_allcandidates.yaml": "1cb16621e3ecd75707a4fae4e49993438cef29c9fd843246cb3458d2f74535e7",
         "data/feature_cache.parquet": "075112ad7a811e7d25d546ab30af8d8568c27c39f5c8fea6705140568f825c0d"}
SEEDS = {"with": range(1, 11), "without": range(1, 4)}                          # the spread over seeds is small
NIGHT = ["plate_low_snr_frac", "plate_snr_std", "plate_elongation_mean", "plate_n_high_snr"]


def fetch():
    for f, sha in FILES.items():
        p = BASE / f
        if not p.exists():
            p.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(URL + f, timeout=300) as r:
                p.write_bytes(r.read())
        assert hashlib.sha256(p.read_bytes()).hexdigest() == sha, f"{f} differs from the file analysed"
    (BASE / "MLscripts/__init__.py").touch()


def label_sets():
    """The four label sets, written as CSVs VASCO's loader reads."""
    lab = pd.read_csv(ROOT / "data/external/combined_training_labels.csv")
    m = sample.labelled()
    nt = pd.read_csv(ROOT / "data/external/nuclear_tests_1949_1957.csv")
    tests = pd.DatetimeIndex(pd.to_datetime(nt[nt.Country.str.upper() == "USA"].Date).dt.normalize().unique())
    window = set(tests) | set(tests - pd.Timedelta(days=1)) | set(tests + pd.Timedelta(days=1))
    near = m.archive_date.isin(window).values
    steve = (lab.label_file == "steve_batch_xlsx").values
    out = {}
    for name, keep in (("all", ~np.zeros(len(lab), bool)), ("no_window", ~near), ("doc250", ~steve), ("doc250_no_window", ~steve & ~near)):
        f = BASE / f"labels_{name}.csv"; lab[keep].to_csv(f, index=False); out[name] = f
    return out


def retrain(labels, night, seed, out):
    sys.path.insert(0, str(BASE / "MLscripts"))
    import pipeline_score as PS
    import pipeline_train as PT
    from common import load_config
    if night == "without":
        def flat(df, group_col="obs_date"):
            for col in NIGHT:
                df[col] = 0.0
            return df
        PT.compute_plate_features = flat; PS.compute_plate_features = flat
    cfg = load_config(BASE / "MLscripts/config_red_allcandidates.yaml")
    work = BASE / "work"
    cfg.update(seed=seed)
    cfg["paths"].update(feature_cache=str(BASE / "data/feature_cache.parquet"), output=str(work), labels=[str(labels)])
    PT.run_training(cfg)
    (work / "calibrator.joblib").unlink(missing_ok=True)
    scored = PS.run_scoring(cfg)
    assert len(scored) == len(sample.ratings())                                  # the catalogue's order
    np.save(out, scored.prob.values.astype(np.float32))


if __name__ == "__main__":
    fetch()
    (BASE / "ratings").mkdir(exist_ok=True)
    jobs = [(name, f, night, seed) for name, f in label_sets().items() for night, seeds in SEEDS.items() for seed in seeds]
    if len(sys.argv) > 1:                                                        # one job, in a fresh process (see main)
        name, f, night, seed = jobs[int(sys.argv[1])]
        retrain(f, night, seed, BASE / f"ratings/{name}_{night}_{seed}.npy")
    else:
        # each retrain in its own process: VASCO's modules are patched per run and keep state between imports
        import subprocess
        for i, (name, _, night, seed) in enumerate(jobs):
            if not (BASE / f"ratings/{name}_{night}_{seed}.npy").exists():
                subprocess.run([sys.executable, __file__, str(i)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"{len(jobs)} sets of ratings in {BASE / 'ratings'}")
