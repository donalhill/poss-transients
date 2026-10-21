"""Checks on the alignment null (alignments.py): copies should remove real lines from the catalogue, and nothing else.
Held out of `make`, because it repeats the slow search; run with `make null-checks`, and its output is kept in the
repository.

1. Nothing else: pairs of transients at each separation, the real catalogue against copies made several ways (10 copies
   each): the page's copies, other move sizes, a Gaussian move with no upper limit, other distances for moving transients
   together, every transient moved on its own by up to 2' (keeping each plate's density pattern but none of its
   clustering, no two closer than the plates resolve), and transients given positions taken from other plates (each
   plate keeping its count).
2. Lines removed: 2,000 artificial tracks added to the catalogue (alignments.with_tracks), copied, and counted as still
   straight (chi2 probability at least 1%, with 1" per transient), for each move size.
3. The results: the counts the page quotes (track score at each level; VASCO's score, arbitrated, at 1e-3 or better),
   real catalogue against 1,000 copies for each way of making them, with the page's seeds.
4. Leaks: for each real track scoring 1e-3 or better (after arbitration), the share of copies holding a candidate
   (before arbitration) made of three or more of its transients and scoring at least as well. A null that kept a real
   track would show it here; a close pair turned back onto its own third transient by chance shows at a low rate.
5. Absorption: 100 artificial tracks added to the catalogue, which is then copied the page's way; the copies' counts
   against those of the plain catalogue's copies. A null that absorbed signal would rise with it.
6. Box size: the track score with boxes of 10", 15" and 20" (the page's), each against copies keeping the clumps (the
   page's) and copies keeping none (every transient moved alone by up to 2'), 1,000 of each; and the share of added
   tracks found as 4-transient tracks at a pass mark of 3 chance ones (the page's copies), diagonal, east-west and in any
   direction, as on the page. A box too small leaves close pairs counted twice, and the two kinds of copy disagree.

python analysis/alignment_null_checks.py [boxes] -> data/alignment_null_checks.json ("boxes" recomputes section 6 alone,
in place)"""
import json
import sys
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import chi2

import alignments as AL

ROOT = Path(__file__).resolve().parents[1]
cat = AL.cat
# (move in arcsec, distance for moving together in arcsec, Gaussian move with this rms per axis and no upper limit,
# closest a transient may land to one outside its clump: None for the moving-together distance)
PAGE = (AL.SHIFT, AL.LINK, False, None)
ALONE = (120.0, 1e-6, False, AL.MINSEP)
SETTINGS = [PAGE, (10.0, AL.LINK, False, None), (40.0, AL.LINK, False, None), (10.0, AL.LINK, True, None), (AL.SHIFT, 8.0, False, None),
            (AL.SHIFT, 40.0, False, None), ALONE]
SEP = np.array([0, 8, 12, 16, 20, 30, 40, 60, 120, 300, 900], float)


def copy(df, seed, shift, link, gauss, sep):
    return AL.clump_copy(df, seed, link=link, shift=shift, gauss=gauss, sep=sep)


def other_plates(df, seed):
    """Every plate keeps its count; each transient takes the position of one on another plate."""
    rng = np.random.default_rng(seed); p = df.plate.values
    j = rng.integers(0, len(df), len(df))
    while (bad := p[j] == p).any():
        j[bad] = rng.integers(0, len(df), bad.sum())
    return df.assign(x=df.x.values[j], y=df.y.values[j])


def pairs(df):
    off = pd.factorize(df.plate)[0] * 1e7                             # keeps plates apart in one tree
    t = cKDTree(np.c_[df.x.values + off, df.y.values])
    return np.diff(t.count_neighbors(t, SEP) - len(df)) / 2


def straight(x, y):
    s = np.linalg.svd(np.c_[x - x.mean(), y - y.mean()], compute_uv=False)
    return chi2.sf(s[-1] ** 2, len(x) - 2) >= 0.01


REAL = []                                                              # member ids of the real tracks, set before forking


INJ = AL.with_tracks(seed=7)                                            # the catalogue with 100 artificial tracks


def counts(job):
    setting, seed = job
    base = INJ if setting == "injected" else cat
    r = AL.score(base if seed is None else copy(base, seed, *(PAGE if setting == "injected" else setting)))
    held = [({i for i in a["ids"]}, a["mu"]) for a in r[6]]
    return AL.tally_t(r[6]), AL.tally(r[5])[3], [any(len(t & h) >= 3 and mu <= m for h, mu in held) for t, m in REAL]


name = lambda st: {"shift_arcsec": st[0], "link_arcsec": st[1] if st[1] >= 1 else 0, "gaussian": st[2], "alone": st == ALONE}

BOXES = (10.0, 15.0, AL.BOX)
SIGS = [(a_, 7 + k) for a_ in (45.0, 0.0, None) for k in range(3)]    # the page's added tracks (alignments.SIGNALS)


def box_job(job):
    kind, seed = job
    df = (AL.clump_copy(cat, seed) if kind == "page" else copy(cat, seed, *ALONE) if kind == "alone"
          else AL.with_tracks(seed=seed, angle=kind))
    r = AL.score(df)
    return AL.tally_t(r[6]), r[4][1]                                    # counts at each level; 4-transient score histogram


def boxes():
    rows = []
    for box in BOXES:
        AL.BOX = box                                                    # read when scoring, so set before forking
        r = AL.score(cat); rt, rh = AL.tally_t(r[6]), r[4][1]
        with get_context("fork").Pool(8) as pool:
            P = pool.map(box_job, [("page", 3000 + s) for s in range(AL.NSHUF)])
            A = pool.map(box_job, [("alone", 3000 + s) for s in range(AL.NSHUF)])
            S = pool.map(box_job, [(a_, k) for a_, k in SIGS])
        Tp, Ta = np.array([p_[0] for p_ in P]), np.array([p_[0] for p_ in A])
        i_ = np.flatnonzero(np.cumsum(np.mean([p_[1] for p_ in P], 0)) <= 3)[-1]          # pass mark: 3 chance 4-transient tracks
        found = lambda a_: round(float(np.mean([np.cumsum(S[j][1] - rh)[i_] for j, (b_, _) in enumerate(SIGS) if b_ == a_])), 1)
        rows.append({"box_arcsec": box, "real": rt, "kept_mean": Tp.mean(0).round(3).tolist(), "alone_mean": Ta.mean(0).round(3).tolist(),
                     "kept_share_ge": [round(float(np.mean(Tp[:, j] >= rt[j])), 3) for j in range(len(rt))],
                     "alone_share_ge": [round(float(np.mean(Ta[:, j] >= rt[j])), 3) for j in range(len(rt))],
                     "found_diagonal": found(45.0), "found_east_west": found(0.0), "found_any": found(None)})
        print(rows[-1], flush=True)
    AL.BOX = BOXES[-1]
    return rows


if __name__ == "__main__" and sys.argv[1:] == ["boxes"]:
    f = ROOT / "data/alignment_null_checks.json"; out = json.loads(f.read_text()); out["boxes"] = boxes()
    f.write_text(json.dumps(out, indent=1)); raise SystemExit
if __name__ == "__main__":
    out = {"separation_arcsec": SEP.tolist(), "pairs_real": pairs(cat).tolist(), "pairs": [], "survival": [], "results": []}
    real = np.array(out["pairs_real"])
    for st in SETTINGS:
        c = np.mean([pairs(copy(cat, 3000 + s, *st)) for s in range(10)], 0)
        out["pairs"].append({**name(st), "copies_over_real": np.round(c / np.maximum(real, 1), 3).tolist()})
    c = np.mean([pairs(other_plates(cat, 3000 + s)) for s in range(10)], 0)
    out["pairs"].append({"other_plates": True, "copies_over_real": np.round(c / np.maximum(real, 1), 3).tolist()})
    print("pairs done", flush=True)

    n0 = len(cat)
    for st in [s_ for s_ in SETTINGS if s_[1] == AL.LINK or s_ == ALONE]:
        k = tot = 0
        for seed in range(20):
            df = AL.with_tracks(seed=seed).reset_index(drop=True)
            moved = copy(df, 9000 + seed, *st).iloc[n0:]
            track = (df.plate.iloc[n0:] != df.plate.iloc[n0:].shift()).cumsum().values   # each added track, in order
            for g in np.unique(track):
                m = track == g; tot += 1; k += straight(moved.x.values[m], moved.y.values[m])
        out["survival"].append({**name(st), "tracks": tot, "still_straight": round(k / tot, 4)})
    print("survival done", flush=True)

    r = AL.score(cat); rt, rv = AL.tally_t(r[6]), AL.tally(r[5])[3]
    out["injected_track"] = AL.tally_t(AL.score(INJ)[6])
    strong = AL.arbitrated(r[6], AL.KEEP); REAL.extend(({i for i in a["ids"]}, a["mu"]) for a in strong)
    out["levels"], out["real_track"], out["real_vasco_1e-3"] = list(AL.LEVELS_T), rt, rv
    with get_context("fork").Pool(8) as pool:
        for st in SETTINGS:
            R = pool.map(counts, [(st, 3000 + s) for s in range(AL.NSHUF)])
            T = np.array([r_[0] for r_ in R]); V = np.array([r_[1] for r_ in R]); K = np.array([r_[2] for r_ in R])
            if st == PAGE:
                out["leaks"] = [{"plate": a["plate"], "r": a["r"], "score": a["mu"], "copies_holding_3": round(float(K[:, j].mean()), 4)}
                                for j, a in enumerate(strong)]
            out["results"].append({**name(st), "track_mean": T.mean(0).round(3).tolist(),
                                   "track_share_ge": [round(float(np.mean(T[:, j] >= rt[j])), 3) for j in range(len(rt))],
                                   "vasco_mean": round(float(V.mean()), 2), "vasco_share_ge": round(float(np.mean(V >= rv)), 3)})
            print(out["results"][-1], flush=True)
        R = pool.map(counts, [("injected", 3000 + s) for s in range(AL.NSHUF)])
        out["absorption"] = {"copies_of_injected_mean": np.array([r_[0] for r_ in R]).mean(0).round(3).tolist(),
                             "copies_of_plain_mean": out["results"][0]["track_mean"]}
        print(out["absorption"], flush=True)
    out["boxes"] = boxes()
    (ROOT / "data/alignment_null_checks.json").write_text(json.dumps(out, indent=1))
