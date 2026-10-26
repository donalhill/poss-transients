"""The page's provenance: the SHA-256 of every file the page was made from: VASCO's public files, the samples built from
them, the text, and the numbers.
usage: python report/provenance.py -> docs/PROVENANCE.json"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
X = ROOT / "data/external"
files = [X / "SUPERVIKTIG_HELAVASCO_validated_v4.csv", X / "scored_catalog.zip", X / "XE_final_random_sources.csv",
         X / "poss_i_observing_dates_strict.csv", X / "combined_training_labels.csv"] + \
        [ROOT / f"data/{n}" for n in ("transients_with_plates.parquet", "plates_corrected.parquet", "plate_times_all.csv", "plate_depth.csv",
                                        "plate_times_nominal.csv", "asteroids_on_plates.parquet", "asteroid_clock.csv", "asteroid_clock_matches.csv",
                                        "plate_time_posterior.csv")] + \
        [ROOT / f"release/{s}/{n}" for s in ("vasco_rated", "vasco_raw") for n in ("transients.parquet", "plates.parquet")] + \
        [ROOT / f"report/{n}" for n in ("copy.md", "briefs.md", "registry.py", "figures.py", "numbers.json", "template.html")] + \
        [ROOT / "docs/index.html", ROOT / "docs/figures.json"]
out = {"files": {str(p.relative_to(ROOT)): sha(p) for p in files if p.exists()}}
(ROOT / "docs/PROVENANCE.json").write_text(json.dumps(out, indent=1))
print(f"docs/PROVENANCE.json: {len(out['files'])} files")
