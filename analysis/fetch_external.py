"""VASCO's public inputs, fetched into data/external/ (not kept in this repository):
  - SUPERVIKTIG_HELAVASCO_validated_v4.csv: the transient catalogue (107,875 transients)
  - poss_i_observing_dates_strict.csv: the survey's observing nights (from Palomar's log, VizieR VI/25)
  - nuclear_tests_1949_1957.csv: the US nuclear tests used in VASCO's analysis
  - XE_final_random_sources.csv: random positions on each plate (the control sample)
  - MLscripts/nuclear_permutation_deciles.py, MLscripts/shadow_deficit_deciles.py: VASCO's own decile analysis code,
    run unchanged by ml_decile_null.py
python analysis/fetch_external.py"""
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
GH = "https://raw.githubusercontent.com/dca-doherty"
FILES = {
    "SUPERVIKTIG_HELAVASCO_validated_v4.csv": f"{GH}/VASCO-ML/main/data/SUPERVIKTIG_HELAVASCO_validated_v4.csv",
    "poss_i_observing_dates_strict.csv": f"{GH}/VASCO-ML/main/data/poss_i_observing_dates_strict.csv",
    "nuclear_tests_1949_1957.csv": f"{GH}/VASCO-ML/main/data/nuclear_tests_1949_1957.csv",
    "XE_final_random_sources.csv": f"{GH}/VASCO-challenge/main/data/XE_final_random_sources.csv",
    "MLscripts/nuclear_permutation_deciles.py": f"{GH}/VASCO-ML/main/MLscripts/nuclear_permutation_deciles.py",
    "MLscripts/shadow_deficit_deciles.py": f"{GH}/VASCO-ML/main/MLscripts/shadow_deficit_deciles.py",
}

if __name__ == "__main__":
    out = ROOT / "data/external"
    out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        f = out / name
        f.parent.mkdir(parents=True, exist_ok=True)
        if f.exists():
            continue
        r = requests.get(url, timeout=300)
        r.raise_for_status()
        f.write_bytes(r.content)
        print(f"{name}: {len(r.content) / 1e6:.1f} MB")
