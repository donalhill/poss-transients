"""Earth's shadow on the real plates: the shadow test's statistics (shadow_stats.py) on the plates the shadow partly
covers. Hypothesis: some transients are glints, which vanish while their location is in Earth's shadow; the rest are
artefacts, which the shadow does not touch. The shadow follows shadow.py's model, seen from Palomar: a glint is a flash
at one unknown instant of the exposure, and its location is in shadow for a share s of the exposure (the shadow slides
across a plate during an exposure; its edge follows the share of the Sun's disc the Earth hides).
The primary statistic is Z (shadow_stats.py), the transients observed in shadow against those expected from each
plate's own total and the hosts' map, quoted as the ratio O / E; its null and the limit on glints come from replicas of
the real shadows on other plates (replicas.py). A plate is left out when no cell is sunlit throughout the exposure
(nothing sunlit to compare with) or when less than MIN_SHADOW of its area is in shadow. Each plate's own ratio is shown
with its 68% range from counting alone (O plus or minus the standard deviation of a count drawn from the map, over E).
Also: the ratio fit and the glint-share fit as checks, the pooled counts E1 and E2, and each plate left out in turn.
POSS_FIT_ONLY=1 with POSS_TIME_JITTER_FILE (each plate's time drawn from its posterior, data/plate_time_draws.parquet):
the statistics alone for each draw -> $POSS_SAMPLE/shadow_fit_only.json.
python analysis/shadow_fit.py -> $POSS_SAMPLE/shadow_fit.json"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402
import shadow_stats as ST                                                    # noqa: E402
from shadow import HALF                                                      # noqa: E402

D = Path(os.environ.get("POSS_SAMPLE", "."))
MIN_SHADOW = float(os.environ.get("POSS_MIN_SHADOW", "0.01"))                # a plate needs this share of its area in shadow


def crossed_plates():
    """Plates some shadow reaches (at its widest, at any time of the exposure: crossed.py)."""
    return set(json.loads((D / "crossed.json").read_text())["crossed"])


def acceptance(pl):
    """The map (mean 1) and the largest s / A it allows (1 / its minimum)."""
    A = np.asarray(pl.acc.iloc[0], float); A = A / A.mean()
    return A, 1.0 / A.min()


def cell_of(xi, eta):
    m = int(round(SH.GX.size ** 0.5))
    return (np.clip(np.floor((eta + HALF) / SH.CELL).astype(int), 0, m - 1) * m
            + np.clip(np.floor((xi + HALF) / SH.CELL).astype(int), 0, m - 1))


def usable(sub):
    """The fit's rule for a shadow on a plate: some of it in shadow, a cell sunlit throughout, 1% of the area in shadow."""
    if not sub.any():
        return "none"
    if not (sub <= 0).all(1).any():
        return "whole"
    return "thin" if sub.mean() < MIN_SHADOW else "ok"


def plate(c, p, sub, A):
    """One plate's record for shadow_stats: its transients c under the shadow of plate p (sub: cells x sub-points)."""
    s = SH.dark(c.xi.values, c.eta.values, p, SH.MAIN)
    return dict(ST.geometry(sub, A), s=s, a=A[cell_of(c.xi.values, c.eta.values)], u=c.w.values.astype(float))


def real_plates(cat, pl, aside=False):
    """The plates the shadow partly covers at their own times, with the plates left out by rule; with aside, also the
    records of those left out (the pooled counts take every plate the shadow reaches), flagged "aside"."""
    A, _ = acceptance(pl)
    P, whole, thin = [], [], []
    for q in sorted(crossed_plates() & set(pl.index)):
        p = pl.loc[q]; sub = SH.sub_dark(p, SH.MAIN); rule = usable(sub)
        if rule in ("whole", "thin"):
            (whole if rule == "whole" else thin).append(q)
        c = cat[cat.plate == q]
        if rule == "none" or c.empty or (rule != "ok" and not aside):
            continue
        P.append(dict(plate(c, p, sub, A), plate=q, sub=sub, n=len(c), aside=rule != "ok"))
    return P, whole, thin


def measure(cat, pl, extra=True):
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)]
    A, qmax = acceptance(pl)
    Pall, whole, thin = real_plates(cat, pl, aside=True); P = [p for p in Pall if not p["aside"]]
    st = ST.stats(P, fits=extra, qmax=qmax)
    out = {"plates": [], "whole": whole, "thin": thin, "fit": {"r": round(st["OE"], 4), "Z": round(st["Z"], 4), "O": round(st["O"], 2), "E": round(st["E"], 2),
           "n": round(float(sum(p["u"].sum() for p in P)), 1), "transients": int(sum(p["n"] for p in P))}}
    sd = np.sqrt(sum(float((p["u"] ** 2).sum()) * (p["Es2"] - p["F"] ** 2) for p in P))
    out["fit"].update(count_lo=round((st["O"] - sd) / st["E"], 4), count_hi=round((st["O"] + sd) / st["E"], 4))
    if not extra:
        return out
    mu, y, alpha = ST.pointing(pl, cat)
    for p in Pall:
        p["mu"] = float(mu[p["plate"]])
    pooled = ST.stats(Pall, totals=(float(y.sum()), len(pl)))
    out["fit"].update(ratio_fit=round(st["r"], 4), glint_share_fit=round(st["g"], 4), E1=round(pooled["E1"], 4), E2=round(pooled["E2"], 4),
                      hybrid=round(ST.fit_hybrid(P, Pall, alpha), 4), dispersion=round(alpha, 4))
    for p in P:
        O, E, _ = ST.score([p]); sd1 = np.sqrt(float((p["u"] ** 2).sum()) * (p["Es2"] - p["F"] ** 2))
        out["plates"].append({"plate": p["plate"], "shadow": round(p["S"], 3), "weighted": round(p["F"], 3), "n": p["n"], "O": round(float(p["s"].sum()), 1),
                              "nw": round(float(p["u"].sum()), 1), "Ow": round(O, 1), "E": round(E, 1), "r": round(O / E, 3),
                              "lo": round(max(O - sd1, 0) / E, 3), "hi": round((O + sd1) / E, 3), "Z": round(ST.score([p])[2], 3),
                              "ratio_fit": round(ST.fit_r([p]), 3), "mu": round(p["mu"], 1)})
    out["pooled_by_plate"] = [{"plate": p["plate"], "n": round(float(p["u"].sum()), 1), "expected": round(p["mu"], 1), "area_in_shadow": round(p["S"], 4),
                               "in_shadow": round(float(p["u"] @ p["s"]), 1), "expected_in_shadow": round(p["mu"] * p["F"], 1), "set_aside": p["aside"]}
                              for p in Pall]
    out["totals"] = {"N": round(float(y.sum()), 1), "J": len(pl)}
    print(pd.DataFrame(out["plates"]).to_string(index=False))
    print("all plates:", out["fit"])
    out["leave_one_out"] = {p["plate"]: {k: round(v, 4) for k, v in ST.stats([P[i] for i in range(len(P)) if i != j]).items() if k in ("OE", "Z")}
                            for j, p in enumerate(P)}
    print("each plate left out:", out["leave_one_out"])
    return out


def fit_only():
    """The statistics with every plate's time drawn from its posterior (POSS_TIME_JITTER_FILE: a table of draws, minutes,
    one column per plate; POSS_TIME_JITTER_N draws) -> $POSS_SAMPLE/shadow_fit_only.json"""
    cat, pl = sample.load()
    m = measure(cat, pl, extra=False)
    J = pd.read_parquet(os.environ["POSS_TIME_JITTER_FILE"]); oe, z = [], []
    for k in range(min(len(J), int(os.environ.get("POSS_TIME_JITTER_N", "200")))):
        dt = pd.to_timedelta(J.iloc[k].reindex(pl.index).fillna(0.0).to_numpy(), unit="min")
        f = measure(cat, pl.assign(ut=pl.ut + dt), extra=False)["fit"]; oe.append(f["r"]); z.append(f["Z"])
    out = {"r": m["fit"]["r"], "Z": m["fit"]["Z"], "time_posterior": {"draws": oe, "z_draws": z, "mean": round(float(np.mean(oe)), 4),
           "sd": round(float(np.std(oe, ddof=1)), 4), "z_sd": round(float(np.std(z, ddof=1)), 4)}}
    (D / "shadow_fit_only.json").write_text(json.dumps(out, indent=1))
    print(out["r"], out["time_posterior"]["sd"], out["time_posterior"]["z_sd"])


if __name__ == "__main__" and os.environ.get("POSS_FIT_ONLY") == "1":
    fit_only()
elif __name__ == "__main__":
    out = measure(*sample.load())
    (D / "shadow_fit.json").write_text(json.dumps(out))
