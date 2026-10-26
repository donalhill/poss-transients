"""The shadow test's result with the plate times' uncertainty: every plate's time drawn from its asteroid posterior
(time_posterior/shadow_fit_only.json) moves O / E; the p-value is given at the real O / E and at O / E one standard
deviation either side, and the 95% upper limit on the glint share is taken at O / E moved 1.645 standard deviations
towards a shortfall (the conservative side), from the replicas' distributions of O / E at each glint share
(replica_draws.npz, replicas.py's primary choice of host).
usage: python diagnostics/shadow_total.py DIR   (writes DIR/shadow_total.json)"""
import json
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
R = json.loads((d / "replicas.json").read_text()); T = json.loads((d / "time_posterior/shadow_fit_only.json").read_text())["time_posterior"]
X = np.load(d / "replica_draws.npz"); null, OEg, G0 = X["null_oe"], X["oe_by_g0"], X["g0"]
r, s = R["real"]["OE"], T["sd"]
p = lambda v: round(float((np.sum(null <= v) + 1) / (len(null) + 1)), 4)


def limit(v):
    b = np.minimum.accumulate((OEg >= v).mean(0))
    return round(float(np.interp(0.05, b[::-1], G0[::-1])), 4) if b[-1] < 0.05 else None


out = {"OE": r, "Z": R["real"]["Z"], "oe_sd_time": s, "z_sd_time": T["z_sd"], "p": p(r), "p_range": [p(r - s), p(r + s)],
       "limit95_stat": limit(r), "limit95": limit(r - 1.645 * s)}
(d / "shadow_total.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
