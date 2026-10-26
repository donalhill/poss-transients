"""Join the parallel parts of the alignment fit (alignment_fit_part{3,4,5,6}.json, POSS_FIT_FLASHES) into alignment_fit.json.
POSS_SAMPLE=DIR python analysis/alignment_fit_merge.py"""
import json
import os
from pathlib import Path

D = Path(os.environ["POSS_SAMPLE"])
parts = {f: json.loads((D / f"alignment_fit_part{f}.json").read_text()) for f in (3, 4, 5, 6)}
out = {k: v for k, v in parts[3].items() if k != "flashes"}; out["flashes"] = {str(f): parts[f]["flashes"][str(f)] for f in parts}
(D / "alignment_fit.json").write_text(json.dumps(out, indent=1))
print({f: (out["flashes"][f]["s_hat"], out["flashes"][f]["p"]) for f in out["flashes"]})
