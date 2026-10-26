"""The numbers the page quotes from the alignment null's checks (alignment_direction.py, alignment_null_tune.py,
alignment_sideband.py), for the page's copies (transients within 20" moving as one piece, every move up to 20").
POSS_SAMPLE=DIR python analysis/alignment_checks_summary.py -> DIR/alignment_checks.json"""
import json
import os
from pathlib import Path

import pandas as pd

D = Path(os.environ["POSS_SAMPLE"])
Dr = json.loads((D / "alignment_direction.json").read_text())
Sb = json.loads((D / "alignment_sideband.json").read_text())["configs"]["link20_shift20"]["band0"]
found = [f[1] for f in Dr["direction"]["found"]]                            # share giving a candidate at 1e-3 or better
A = Dr["absorption"]["by_hits"]
fnd = [sum(A[h]["found"][i] for h in A) for i in (0, 1)]; kept = [sum(A[h]["kept_by_copies"][i] for h in A) for i in (0, 1)]
P = Dr["pairs"]; cp = Dr["close_pairs"]
T = pd.read_parquet(D / "transients.parquet", columns=["plate"])
out = {"transients_searched": int((T.plate != "").sum()), "direction": {"angles": len(found), "tracks_per_angle": Dr["tracks_per_angle"], "found_min": min(found), "found_max": max(found)},
       "absorption": {"per_100_1e-4": 100 * kept[0] / fnd[0], "per_100_1e-3": 100 * kept[1] / fnd[1],
                      "excess_shown_1e-4": 1 - kept[0] / fnd[0], "excess_shown_1e-3": 1 - kept[1] / fnd[1],
                      "tracks": 3 * Dr["absorption"]["tracks_per_size"]},
       "pairs": {"within_20_catalogue": sum(P["catalogue"][:4]), "within_20_unclustered": sum(P["alone"][:4]),
                 "within_20_ratio": sum(P["catalogue"][:4]) / sum(P["alone"][:4])},
       "close_pairs": {"kept_1e-3": cp["box20_clumps"][0], "split_1e-3": cp["box20_alone20"][0],
                       "kept_1e-4": cp["box20_clumps"][1], "split_1e-4": cp["box20_alone20"][1], "copies": cp["copies"]},
       "sideband": {"real_1e-3": Sb["real"][0], "copies_1e-3": Sb["copies_mean"][0], "real_1e-2": Sb["real"][1],
                    "copies_1e-2": Sb["copies_mean"][1], "share_copies_ge_real_1e-3": Sb["share_copies_ge_real"][0]}}
(D / "alignment_checks.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
