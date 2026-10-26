"""One track candidate scoring 1 in 1,000 or better, to show on the page how the track score works: its transients and the others around it (arcmin from the box centre, east to the left), each
transient's distance from the fitted line, and what goes into the score. The score is recomputed from these and must
match the search's.
The candidate, unless named: the best-scoring candidate of four transients with chi2/ndof at most 2 (of three, if none).
POSS_SAMPLE=DIR python analysis/alignment_example.py [PLATE R] -> DIR/alignment_example.json"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).parent))
want = sys.argv[1:3]
sys.argv = sys.argv[:1]
import alignments as AL                                                      # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
AJ = json.loads((D / "alignments.json").read_text())
C = AJ["track"]["candidates"]
if want:
    EX = next(c for c in C if c["plate"] == want[0] and c["r"] == int(want[1]))
else:
    EX = min([c for c in C if c["r"] == 4 and c["chi2ndof"] <= 2] or [c for c in C if c["r"] == 3 and c["chi2ndof"] <= 2], key=lambda c: c["score"])
plate, r = EX["plate"], EX["r"]
scored = AL.cat
P = np.c_[EX["x"], EX["y"]]; ctr = P.mean(0); Vt = np.linalg.svd(P - ctr)[2]
al, pp = (P - ctr) @ Vt[0], (P - ctr) @ Vt[1]
L = float(np.hypot(*(P[:, None] - P[None]).transpose(2, 0, 1)).max()); side = max(600.0, 1.5 * L)
mid = (P[al.argmin()] + P[al.argmax()]) / 2
g = scored[scored.plate == plate]; Q = np.c_[g.x.values, g.y.values] - mid
N = int(((np.abs(Q) <= side / 2).all(1)).sum()); n = N / (side / 60) ** 2; w = float(np.abs(pp).max())
lam = np.pi * N * n * (2 * n * max(w, 1.0) / 60) ** (r - 2) * (L / 60) ** r / (r * math.factorial(r - 2))
assert abs(lam / EX["score"] - 1) < 0.01, (lam, EX["score"])
chi = float((pp ** 2).sum())                                                 # sigma = 1"
tol = 0.005 + 2 * np.abs(pp).sum() * 0.005 * np.sqrt(2) / (r - 2)           # the stored chi2/ndof and positions are rounded to 0.01
assert abs(chi / (r - 2) - EX["chi2ndof"]) < tol, (chi / (r - 2), EX["chi2ndof"], tol)
view = (np.abs(Q) <= 0.8 * side).all(1)
out = {"plate": plate, "side": round(side / 60, 3), "r": r, "N": N, "L": round(L / 60, 2), "w": round(w, 2), "score": EX["score"],
       "chi2": round(chi, 2), "p_chi2": round(float(chi2.sf(chi, r - 2)), 3), "sigma": 1.0, "hmax": AL.HMAX, "occupied": AL.BOX,
       "lmax": round(AL.LMAX / 60), "p_min": round(float(chi2.sf(AL.CHI2MAX[r - 2], r - 2)), 3),
       "members": np.round(np.c_[-(P[:, 0] - mid[0]), P[:, 1] - mid[1]] / 60, 3).tolist(),
       "others": np.round(np.c_[-Q[view, 0], Q[view, 1]] / 60, 3).tolist(),
       "res": np.round(np.c_[(al - al.min()) / 60, pp], 3).tolist()}
(D / "alignment_example.json").write_text(json.dumps(out))
print({k: v for k, v in out.items() if k not in ("members", "others", "res")})
