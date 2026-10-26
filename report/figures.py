"""The data behind every figure of the page, each from named analysis outputs, with the SHA-256 of every file it reads:
the page draws its charts from this file alone, so a chart cannot drift from the numbers in the text. Writes
report/figures.json: {"figures": {id: data}, "sources": {id: {file: sha256}}}. The briefs name each caption's figure
(field "figure"), and report/check.py holds them to this file.
usage: python report/figures.py"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VS = ROOT / "release/vasco_rated"                                            # the tests' sample (build_sample.py)


def hist(x, lo, hi, n=60):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    h, e = np.histogram(np.clip(x, lo, hi), bins=n, range=(lo, hi))
    return {"edges": np.round(e, 4).tolist(), "counts": h.tolist()}


def shadow_null(d):
    """O / E on the replicas with no glints and with glints at the 95% limit (replicas.py, the primary choice of host), in
    80 bins over the middle 99.8% of both, the few beyond counted in the end bins."""
    R = json.loads((d / "replicas.json").read_text()); X = np.load(d / "replica_draws.npz"); T = json.loads((d / "shadow_total.json").read_text())
    r = R["real"]["OE"]; m = R["matches"]["count_dec"]; null, lim = X["null_oe"], X["limit_oe"]
    pool = np.r_[null, lim, r]; e = np.linspace(np.quantile(pool, 0.001), np.quantile(pool, 0.999), 81)
    h = lambda v: np.histogram(np.clip(np.asarray(v, float), e[0], e[-1] - 1e-9), bins=e)[0].tolist()
    return {"r": round(r, 3), "p_low": T["p"], "null": m["null"]["OE"], "plates": len(json.loads((d / "shadow_fit.json").read_text())["plates"]),
            "edges": np.round(e, 4).tolist(), "counts": h(null), "log": False, "replicas": True,
            "xlabel": "Transients in shadow, observed over expected",
            "limit": {"g0": round(float(X["limit_g0"]), 3), "share_as_high": round(float(np.mean(lim >= r)), 3), "counts": h(lim)}}


def alignments(npz):
    """The track score's distributions (straight candidates only, arbitrated within each size), the catalogue against the
    copies, for tracks of 3 to 6 transients: the distributions the fit uses."""
    from scipy.stats import norm
    d = np.load(npz); C, n = d["null_l"], d["real_l"]                     # in 0.2-dex bins
    C, n, edges = C.reshape(*C.shape[:-1], -1, 2).sum(-1), n.reshape(*n.shape[:-1], -1, 2).sum(-1), d["bins"][::2]; N = len(C); mu = C.mean(0)
    # how unusual each bin's count is among the copies: the share holding as many (or as few), as a signed normal deviate;
    # with N copies it reaches at most about 3.1
    p_hi, p_lo = ((C >= n).sum(0) + 1) / (N + 1), ((C <= n).sum(0) + 1) / (N + 1)
    z = np.where(n > mu, norm.isf(np.minimum(p_hi, 0.5)), -norm.isf(np.minimum(p_lo, 0.5)))
    return {"edges": np.round(edges, 4).tolist(), "real": n.tolist(),
            "copies_mean": np.round(mu, 3).tolist(),
            "copies_q16": np.quantile(C, 0.16, axis=0).tolist(), "copies_q84": np.quantile(C, 0.84, axis=0).tolist(),
            "z": np.round(z, 2).tolist(), "points": [3, 4, 5, 6]}


def align_direction(f):
    """The share of artificial glint tracks (4-6 hits, 2-10', 1" scatter) the search finds as a candidate scoring 1e-3 or
    better, at each fixed angle from east-west, with its binomial range."""
    D = json.loads(f.read_text()); k = D["tracks_per_angle"]; p = np.array([v[1] for v in D["direction"]["found"]])
    return {"angle": D["direction"]["angle_deg"], "found": np.round(p, 4).tolist(), "err": np.round(np.sqrt(p * (1 - p) / k), 4).tolist(), "tracks": k}


def align_method(f):
    """One track candidate and what goes into its score (analysis/alignment_example.py), for the method diagram."""
    return json.loads(f.read_text())


def nuclear(f):
    j = json.loads(f.read_text())
    return {"ratio": j["ratio"], "calendars": hist(j["calendar_draws"], 0.6, 1.4, 40), "random_nights": hist(j["random_draws"], 0.6, 1.4, 40),
            "n_calendars": j["calendars"], "p_calendar": j["calendar_null"]["p"], "p_random": j["random_nights_null"]["p"],
            "sigma": j["sigma_calendar"], "plates": j["plates"], "plates_near": j["plates_near"]}


def site(d):
    """site_data.json: every plate's outline, date, and transients, and the shadow's plates."""
    return json.loads((d / "site_data.json").read_text())


def sky(d):
    """Every plate's outline, date, and transients."""
    return {"plates": [{k: p[k] for k in ("id", "N", "date", "poly")} for p in site(d)["plates"]]}


def sky_shadow(d):
    """The same, with the shadow's edge at mid-exposure on each plate it partly covers."""
    return {"plates": site(d)["plates"]}


def shadow_plates(d):
    """Each plate the shadow partly covers, up close: its shadow share by cell, the edge at the
    start and end, and its transients."""
    S = site(d)["shadow"]; j = json.loads((d / "shadow_fit.json").read_text())
    return {"half": S["half"], "plates": S["plates"], "left": sorted(j["whole"]), "thin": sorted(j.get("thin", []))}   # the plates left out, by reason


def shadow_forest(d):
    """Each plate's own ratio with its counting range and its date, and the common fit."""
    j = json.loads((d / "shadow_fit.json").read_text()); day = {q["id"]: q["date"] for q in site(d)["plates"]}
    return {"plates": [dict(p, date=day.get(p["plate"])) for p in j["plates"]], "fit": j["fit"], "left": sorted(j["whole"]), "thin": sorted(j.get("thin", []))}


def clock(matches, posterior, plates):
    """The asteroid clock: every asteroid's offset along its track (minutes from the written time), stacked over the
    sample's plates (plates: its plates table), the survey's offset, and the plates whose shared written time was wrong by
    more than 5 minutes."""
    R = pd.read_csv(matches); P = pd.read_csv(posterior, index_col=0)
    mine = set(pd.read_parquet(plates).index); R = R[R.plate.isin(mine)]; P = P[P.index.isin(mine)]
    d = R.dt_min[R.dt_min.abs() <= 10]                                       # half-minute bins, so the peak sits on the survey offset
    h, e = np.histogram(d, bins=np.arange(-10, 10.01, 0.5))
    far = R[R.dt_min.abs() > 10].groupby("plate").dt_min.apply(lambda x: sorted(round(float(v), 1) for v in x))   # beyond the axis
    off = P[P.tau_median.abs() > 5].sort_values("tau_median")
    return {"edges": e.tolist(), "counts": h.tolist(), "survey": round(float((P.tau_median - P.from_prior_centre).median()), 2),
            "beyond": [{"plate": q, "dt": v} for q, v in far.items()],
            "outliers": [{"plate": q, "dt": round(float(r.tau_median), 1), "lo": round(float(r.tau_q16), 1), "hi": round(float(r.tau_q84), 1), "n": int(r.on_line)} for q, r in off.iterrows()]}


def dates_table(corrections, plates):
    """The archive's dates by run of plate numbers, and the plates that differ from their run, from each plate's corrected
    date (plate_date_corrections.csv) against its header's: a run is a stretch of serial numbers whose archive dates are
    all correct or all one day early; a plate is an exception when its offset is neither, or when it alone (or with one
    other) interrupts a run. Runs are found over every dated plate; only the sample's plates (plates: its plates table)
    are listed as exceptions."""
    C = pd.read_csv(corrections, parse_dates=["corrected_date"])
    mine = set(pd.read_parquet(plates).index)
    C = C[C.corrected_date.notna()].copy(); C["serial"] = pd.to_numeric(C.label, errors="coerce")
    C = C.dropna(subset=["serial"]).sort_values("serial").reset_index(drop=True)
    C["final"] = C.corrected_date; C["off"] = (C.final - pd.to_datetime(C.header_date)).dt.days
    run = (C.off != C.off.shift()).cumsum(); g = C.groupby(run).off.agg(["first", "size"])
    prev, nxt = g["first"].shift(), g["first"].shift(-1)
    exc_run = ~g["first"].isin([0, 1]) | ((g["size"] < 3) & (prev == nxt))
    C["exc"] = run.map(exc_run)
    K = C[~C.exc]; r2 = (K.off != K.off.shift()).cumsum()
    m = lambda t: t.strftime("%b %Y")
    runs = [{"serials": (f"up to E-{int(k.serial.max())}" if n == 0 else f"E-{int(k.serial.min())} to E-{int(k.serial.max())}"),
             "taken": m(k.final.min()) if m(k.final.min()) == m(k.final.max()) else f"{m(k.final.min())} – {m(k.final.max())}",
             "archive": "Correct" if k.off.iloc[0] == 0 else "One day early"} for n, (_, k) in enumerate(K.groupby(r2))]
    def note(r):
        rest = int(K.off[K.serial < r.serial].iloc[-1]) if (K.serial < r.serial).any() else int(K.off.iloc[0])
        if r.off == 1 and rest == 0: return "One day early; the rest of its run correct"
        if r.off == 0 and rest == 1: return "Correct; the rest of its run one day early"
        return f"Archive date out by {abs(int(r.off))} days"
    exc = [{"plate": r.plate, "serial": f"E-{int(r.serial)}", "archive": r.header_date, "correct": r.final.strftime("%Y-%m-%d"), "note": note(r)}
           for _, r in C[C.exc].iterrows() if r.plate in mine]
    return {"runs": runs, "exceptions": exc}


STEPS = {"vasco": "VASCO’s construction: their nights and times, seen from the centre of the Earth, a sharp edge at one instant",
         "palomar": "Seen from Palomar", "times": "Every plate on its corrected date and calibrated time",
         "moving": "The shadow moving over the exposure, its edge from the Sun’s disc", "rated": "Each transient counted by its rating"}


def vasco_repeat(f):
    j = json.loads(f.read_text())
    rows = [[STEPS[k], f"{s['in_shadow']:,.0f}", f"{s['expected']:,.0f}", f"{s['ratio']:.2f}", f"{s['binomial_z']:.1f}"]
            for k, s in zip(STEPS, j["steps"])]
    p = j["shuffled"]["rated"]["p_as_low"]
    rows.append(["Judged against shuffled plate times", "", "", f"{j['shuffled']['rated']['ratio']:.2f}",
                 f"1 in {round(1 / p):,} shuffles as low"])
    return {"head": ["Step", "In shadow", "Expected", "Ratio", "Standard deviations (binomial)"], "rows": rows}


def catalogue(part):
    """A chart of VASCO's catalogue before any test, from report/vasco_stats.json's "fig" (the counts and fit the text quotes)."""
    def build(f):
        F = json.loads(f.read_text())["fig"]
        return {"spread": lambda: {"plates": [{"N": p["N"]} for p in F["plates"]], "mean": F["mean"], "lo": F["lo"], "hi": F["hi"],
                                   "in_band": json.loads(f.read_text())["in_band"]},
                "edge_map": lambda: F["map"], "edge": lambda: F["radial"],
                "crowd": lambda: {"plates": [{k: p[k] for k in ("id", "N", "b")} for p in F["plates"]], "crowd_fit": F["crowd_fit"]},
                "pointing": lambda: {"plates": [{k: p[k] for k in ("id", "N", "E")} for p in F["plates"]], "calibration": F["calibration"]}}[part]()
    return build


FIGS = {
    "fig.sky": (sky, [VS]),
    "fig.sky_shadow": (sky_shadow, [VS]),
    "fig.shadow_plates": (shadow_plates, [VS]),
    "fig.glint_model": (lambda d: site(d)["model"], [VS]),
    "fig.shadow_forest": (shadow_forest, [VS]),
    "fig.clock": (clock, [ROOT / "data/asteroid_clock_matches.csv", ROOT / "data/plate_time_posterior.csv", VS / "plates.parquet"]),
    "fig.dates_table": (dates_table, [ROOT / "data/plate_date_corrections.csv", VS / "plates.parquet"]),
    "fig.shadow_null": (shadow_null, [VS]),
    "fig.alignments": (alignments, [VS / "alignments.npz"]),
    "fig.align_direction": (align_direction, [VS / "alignment_direction.json"]),
    "fig.align_method": (align_method, [VS / "alignment_example.json"]),
    "fig.nuclear": (nuclear, [VS / "nuclear.json"]),
    "fig.vasco_repeat": (vasco_repeat, [VS / "vasco_repeat.json"]),
    **{f"fig.vasco_{k}": (catalogue(k), [ROOT / "report/vasco_stats.json"]) for k in ("spread", "edge_map", "edge", "crowd", "pointing")},
}
NEEDS = {"fig.shadow_null": ["replicas.json", "replica_draws.npz", "shadow_fit.json", "shadow_total.json"],     # files read inside folder inputs
         "fig.sky": ["site_data.json"], "fig.sky_shadow": ["site_data.json"], "fig.shadow_plates": ["site_data.json"], "fig.glint_model": ["site_data.json"],
         "fig.shadow_forest": ["shadow_fit.json", "site_data.json"]}


def needs(fid, ins):
    return [p / n for p in ins for n in NEEDS[fid]] if fid in NEEDS else list(ins)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


if __name__ == "__main__":
    figs, sources = {}, {}
    for fid, (fn, ins) in FIGS.items():
        need = needs(fid, ins)
        missing = [str(p.relative_to(ROOT)) for p in need if not p.exists()]
        if missing:
            raise SystemExit(f"{fid}: inputs not yet produced: {missing}")
        figs[fid] = fn(*ins)
        sources[fid] = {str(p.relative_to(ROOT)): sha(p) for p in need}
    out = ROOT / "report/figures.json"
    out.write_text(json.dumps({"figures": figs, "sources": sources}))
    print(f"{out.relative_to(ROOT)}: {len(figs)} figures")
