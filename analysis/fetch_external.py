"""VASCO's public inputs, fetched into data/external/ (not kept in this repository):
  - SUPERVIKTIG_HELAVASCO_validated_v4.csv: the transient catalogue (107,875 transients)
  - poss_i_observing_dates_strict.csv: the survey's observing nights (from Palomar's log, VizieR VI/25)
  - nuclear_tests_1949_1957.csv: the US nuclear tests used in VASCO's analysis (it lacks the eleven shots of Operation
    Upshot-Knothole, 1953; used only to reproduce their result, the full list being data/us_shot_times.csv)
  - XE_final_random_sources.csv: random positions on each plate (the control sample)
Each is fetched from a fixed commit and checked against the file analysed.
python analysis/fetch_external.py"""
import hashlib
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
GH = "https://raw.githubusercontent.com/dca-doherty"
ML, CH = "418e8c89b2d3e7c5f56ab75293847ede54a5a445", "1085b8e5c9af933cbbdf8135142b04435f7de8b5"   # commits used (Oct 2026)
FILES = {   # name: (url at a fixed commit, SHA-256 of the file analysed)
    "SUPERVIKTIG_HELAVASCO_validated_v4.csv": (f"{GH}/VASCO-ML/{ML}/data/SUPERVIKTIG_HELAVASCO_validated_v4.csv",
                                               "aff3e3b55cb39b59549ca5cb536fedf3a90f5282d21da8e2f747e4550d6731f4"),
    "poss_i_observing_dates_strict.csv": (f"{GH}/VASCO-ML/{ML}/data/poss_i_observing_dates_strict.csv",
                                          "da600c77bf84184591ac4ce3abbf61e34dcf8b7434534dd1e68a2dcb58b6e7a7"),
    "nuclear_tests_1949_1957.csv": (f"{GH}/VASCO-ML/{ML}/data/nuclear_tests_1949_1957.csv",
                                    "bfdce719631682c86a6a2456a689b7ad32c56333a52560ceebf2842b748c8694"),
    "XE_final_random_sources.csv": (f"{GH}/VASCO-challenge/{CH}/data/XE_final_random_sources.csv",
                                    "73d5a633bbe9e15dffde60bb3b0b177ccd1fd36d808c95b55dc1b727eaa52bcc"),
}

if __name__ == "__main__":
    out = ROOT / "data/external"
    out.mkdir(parents=True, exist_ok=True)
    for name, (url, sha) in FILES.items():
        f = out / name
        if not f.exists():
            r = requests.get(url, timeout=300)
            r.raise_for_status()
            f.write_bytes(r.content)
            print(f"{name}: {len(r.content) / 1e6:.1f} MB")
        assert hashlib.sha256(f.read_bytes()).hexdigest() == sha, f"{name} differs from the file analysed"
