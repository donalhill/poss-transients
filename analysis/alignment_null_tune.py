"""The copies' grouping distance and move size, chosen on what the catalogue shows without its lines (alignments.clump_copy: transients within
LINK of one another move as one piece, rotated and shifted by up to SHIFT; lone transients shifted by up to SHIFT), on
the sample (POSS_SAMPLE). The real catalogue is never searched for lines here.
For each (LINK, SHIFT, BOX), BOX the length within which transients along a track count once (set to LINK, the
distance over which close pairs are unusually common), against the catalogue:
- pairs of transients by separation, 5" bins to 60" (5 copies), and the mismatch over 10-60": sum of
  (copies - catalogue)^2 / catalogue;
- connected triples: pairs of neighbours within 60" of a transient, summed (a three-point count, almost none of them lines);
- what the copies take in of a real track: single artificial tracks (4-6 hits, 2-10', 1" scatter, any direction) on
  their plate, the candidates they add at 1e-4 and 1e-3 or better, and what a copy of the plate gains when the track's
  own transients move as the copies move them (alignment_fit.library's construction);
- the copies' mean number of track-score candidates at 1e-3 and 1e-4 or better, arbitrated across sizes (40 copies).
POSS_SAMPLE=DIR python analysis/alignment_null_tune.py -> DIR/alignment_null_tune.json"""
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
from alignment_direction import add_track, gained, by_plate                      # noqa: E402

D = Path(os.environ["POSS_SAMPLE"])
cat = AL.cat[["plate", "x", "y"]]
CONFIGS = ((20.0, 20.0), (20.0, 12.0), (25.0, 20.0), (25.0, 25.0), (30.0, 20.0), (30.0, 30.0))   # (LINK = BOX, SHIFT)
NABS, NCOPY, NPAIR = 600, 40, 5


def pairs_triples(df):
    off = pd.factorize(df.plate)[0] * 1e7; t = cKDTree(np.c_[df.x.values + off, df.y.values])
    d = t.sparse_distance_matrix(t, 60.0, output_type="ndarray")["v"]
    nb = np.array([len(v) - 1 for v in t.query_ball_point(t.data, 60.0)])
    return (np.histogram(d[d > 0], np.arange(0, 60.001, 5))[0] // 2).tolist(), int((nb * (nb - 1) // 2).sum())


def setup(link):
    AL.LINK, AL.BOX = link, link


def absorption(job):
    (link, shift), seed = job; setup(link); rng = np.random.default_rng(seed); L, A = [], []
    for i in range(NABS):
        q, tr = add_track(rng, None, (4, 6)); g = by_plate[q]
        L.append(gained(AL.score(g)[4], AL.score(pd.concat([g, tr]))[4])[0])
        cp = AL.clump_copy(g, 7_000_000 + i, link=link, shift=shift); tc = AL.clump_copy(tr, 8_000_000 + i, link=link, shift=shift)
        A.append(gained(AL.score(cp)[4], AL.score(pd.concat([cp, tc]))[4])[0])
    return (link, shift), np.sum(L, 0).tolist(), np.sum(A, 0).tolist()


def copy_job(job):
    (link, shift), seed = job; setup(link)
    df = AL.clump_copy(cat, seed, link=link, shift=shift)
    out = {"tally": [float(v) for v in AL.tally_t(AL.score(df)[6])]}
    if seed < 3000 + NPAIR:
        out["pairs"], out["triples"] = pairs_triples(df)
    return (link, shift), out


if __name__ == "__main__":
    AL.score(by_plate[sorted(by_plate)[0]])
    with Pool(8) as pool:
        ab = pool.map(absorption, [(c, 700 + i) for i, c in enumerate(CONFIGS)])
        cj = pool.map(copy_job, [(c, 3000 + k) for c in CONFIGS for k in range(NCOPY)])
    cp, ct = pairs_triples(cat)
    res = {}
    for s in CONFIGS:
        J = [o for s_, o in cj if s_ == s]; P = np.mean([o["pairs"] for o in J if "pairs" in o], 0); T = np.mean([o["triples"] for o in J if "triples" in o])
        l, a = next((l_, a_) for s_, l_, a_ in ab if s_ == s)
        res[f"link{int(s[0])}_shift{int(s[1])}"] = {"pairs": P.round(1).tolist(), "pair_mismatch_10_60": round(float(np.sum((P[2:] - np.array(cp[2:])) ** 2 / np.array(cp[2:]))), 1),
                            "triples": round(float(T)), "found": l, "kept_by_copies": a,
                            "kept_per_100_found": [round(100 * a_ / l_, 1) for a_, l_ in zip(a, l)],
                            "copies_mean_1e-3_1e-4": np.mean([o["tally"] for o in J], 0).round(2).tolist()[:2]}
    out = {"copies": NCOPY, "tracks": NABS, "catalogue": {"pairs": cp, "triples": ct}, "configs": res}
    (D / "alignment_null_tune.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
