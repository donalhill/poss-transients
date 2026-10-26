"""Checks on the alignment search (alignments.py), on the sample (POSS_SAMPLE):
1. Direction: single artificial glint tracks (4-6 hits at least 8" apart, 2-10' long, 1" scatter, as
   alignments.with_tracks) added to a random plate at a fixed angle from east-west (the plate's tangent plane, x east),
   the plate scored with and without the track; a track is found when the plate gains a track-score candidate at 1e-3
   (or 1e-4) or better in any size.
2. Absorption: what the copies take in of a real track. For single tracks in any direction (4, 5 and 6 hits), the
   candidates the track adds at 1e-4 and 1e-3 or better, and what a copy of its plate made the page's way gains when
   the track's own transients are moved as the copies move them (as alignment_fit.library).
3. Pairs: transients per separation (5" bins to 60") in the catalogue, in copies made the page's way (clumps kept), and
   in copies with every transient moved alone by up to 2' (no two closer than 8"), and moved alone by up to 20".
4. Close pairs, signal or background: the copies' mean number of track-score candidates at 1e-3 and 1e-4 or better,
   arbitrated across sizes, with copies keeping the clumps, copies moving every transient alone by up to 2', and copies
   moving every transient alone by up to 20" (the clumps split, the larger-scale clustering kept), for boxes of 10" and
   20" (the page's). Equal means that the background does not depend on whether close pairs are features of the plates
   or flashes of one object.
POSS_SAMPLE=DIR python analysis/alignment_direction.py -> DIR/alignment_direction.json"""
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
sys.argv = sys.argv[:1]
import alignments as AL                                                      # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
cat = AL.cat[["plate", "x", "y"]]
by_plate = {q: g for q, g in cat.groupby("plate")}
plates = np.array(sorted(by_plate))
ANG = (0, 1, 2, 3, 5, 10, 20, 30, 45, 60, 70, 80, 85, 87, 88, 89, 90)
NDIR, NABS, NCOPY = (20, 20, 2) if os.environ.get("QUICK") else (600, 1500, 100)
LV = np.searchsorted(AL.BINS, [-4.0, -3.0]) - 1                            # Hl bins wholly at or below 1e-4, 1e-3


def add_track(rng, angle, hits):
    """One track on a random plate, drawn as alignments.with_tracks draws it."""
    q = rng.choice(plates); r = int(rng.integers(hits[0], hits[1] + 1)); L = rng.uniform(120, 600)
    th = rng.uniform(0, np.pi) if angle is None else np.radians(angle) * rng.choice([-1, 1])
    x0, y0 = rng.uniform(-2.9 * 3600 + L, 2.9 * 3600 - L, 2)
    s = np.sort(rng.uniform(0, L, r))
    while np.min(np.diff(s)) < AL.MINSEP:
        s = np.sort(rng.uniform(0, L, r))
    perp, along = rng.normal(0, 1.0, r), rng.normal(0, 1.0, r)
    return q, pd.DataFrame({"plate": q, "x": x0 + (s + along) * np.cos(th) - perp * np.sin(th),
                            "y": y0 + (s + along) * np.sin(th) + perp * np.cos(th)}, index=range(30_000_000, 30_000_000 + r))


def gained(base, withtr):
    """Candidates gained at 1e-4 and at 1e-3 or better, summed over sizes; and whether any size gained one."""
    d = withtr - base; c = np.cumsum(d, axis=1)
    return [float(c[:, i].sum()) for i in LV], [bool((c[:, i] > 0).any()) for i in LV]


def direction(job):
    angle, seed = job; rng = np.random.default_rng(seed); out = []
    for _ in range(NDIR):
        q, tr = add_track(rng, angle, (4, 6)); g = by_plate[q]
        out.append(gained(AL.score(g)[4], AL.score(pd.concat([g, tr]))[4])[1])
    return angle, np.mean(out, axis=0).tolist()


def absorption(job):
    hits, seed = job; rng = np.random.default_rng(seed); L, A = [], []
    for i in range(NABS):
        q, tr = add_track(rng, None, (hits, hits)); g = by_plate[q]
        L.append(gained(AL.score(g)[4], AL.score(pd.concat([g, tr]))[4])[0])
        cp = AL.clump_copy(g, 7_000_000 + 10 * hits + i); tc = AL.clump_copy(tr, 8_000_000 + 10 * hits + i)
        A.append(gained(AL.score(cp)[4], AL.score(pd.concat([cp, tc]))[4])[0])
    return hits, np.sum(L, 0).tolist(), np.sum(A, 0).tolist()


def copy_counts(job):
    box, kind, seed = job; AL.BOX = box
    df = AL.clump_copy(cat, seed) if kind == "clumps" else AL.clump_copy(cat, seed, link=1e-3, shift=120.0 if kind == "alone" else AL.SHIFT, sep=AL.MINSEP)
    t = AL.tally_t(AL.score(df)[6])
    return box, kind, t


def pair_hist(df):
    e = np.arange(0, 60.001, 5); off = pd.factorize(df.plate)[0] * 1e7
    d = cKDTree(np.c_[df.x.values + off, df.y.values]).sparse_distance_matrix(cKDTree(np.c_[df.x.values + off, df.y.values]), 60.0, output_type="ndarray")["v"]
    return (np.histogram(d[d > 0], e)[0] // 2).tolist()


if __name__ == "__main__":
    AL.score(by_plate[plates[0]])                                         # compile once before forking
    with Pool(8) as pool:
        dirs = pool.map(direction, [(a, 500 + i) for i, a in enumerate(ANG)])
        absb = pool.map(absorption, [(h, 900 + h) for h in (4, 5, 6)])
        cc = pool.map(copy_counts, [(b, k, 3000 + s) for b in (10.0, 20.0) for k in ("clumps", "alone", "alone20") for s in range(NCOPY)])
    out = {"tracks_per_angle": NDIR, "levels": ["1e-4", "1e-3"],
           "direction": {"angle_deg": [a for a, _ in dirs], "found": [f for _, f in dirs]},
           "absorption": {"tracks_per_size": NABS, "by_hits": {str(h): {"found": l, "kept_by_copies": a} for h, l, a in absb}},
           "pairs": {"edges_arcsec": list(range(0, 61, 5)), "catalogue": pair_hist(cat),
                     "copies": np.mean([pair_hist(AL.clump_copy(cat, 3000 + s)) for s in range(5)], 0).round(1).tolist(),
                     "alone": np.mean([pair_hist(AL.clump_copy(cat, 3000 + s, link=1e-3, shift=120.0, sep=AL.MINSEP)) for s in range(5)], 0).round(1).tolist(),
                     "alone20": np.mean([pair_hist(AL.clump_copy(cat, 3000 + s, link=1e-3, shift=AL.SHIFT, sep=AL.MINSEP)) for s in range(5)], 0).round(1).tolist()},
           "close_pairs": {"levels_log10": list(AL.LEVELS_T), "copies": NCOPY,
                           **{f"box{int(b)}_{k}": np.mean([t for b_, k_, t in cc if b_ == b and k_ == k], 0).round(2).tolist() for b in (10.0, 20.0) for k in ("clumps", "alone", "alone20")}}}
    (D / "alignment_direction.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "pairs"}, indent=1))
