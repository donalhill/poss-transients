"""The plates Earth's shadow can reach (Palomar's shadow at its widest, at any time of the exposure, anywhere on the
plate's central square), and how far the others stay clear of it.
POSS_SAMPLE=DIR python analysis/crossed.py -> DIR/crossed.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402

d = Path(os.environ["POSS_SAMPLE"])
_, plates = sample.load()
crossed = {q for q, p in plates.iterrows() if SH.cell_dark(p, SH.wide("topo")).any()}
allp = [q for q in plates.index if q not in crossed]
clear = pd.Series({q: SH.clearance(plates.loc[q], SH.wide("topo")) for q in allp})
assert (clear > 0).all()
(d / "crossed.json").write_text(json.dumps({"crossed": sorted(crossed), "hosts": len(allp), "closest_deg": round(float(clear.min()), 2)}, indent=1))
print(f"{len(crossed)} plates crossed, {len(allp)} hosts, closest host {clear.min():.2f} deg")
