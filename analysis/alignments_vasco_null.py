"""Which null is fair for VASCO's score? For each candidate null: (1) does it reproduce the real catalogue in the bulk
(scores 0.01 to 1, where chance groupings dominate), by size, split into straight (chi2/ndof <= 2) and loose candidates;
(2) does it destroy 150 injected lines like VASCO's candidates (4-6 transients, 2-10' long, 1" scatter); (3) what it gives
in the tail (lines, each counted once, scoring 1e-3 or better). Nulls: every transient moved independently by up to
3, 5, 10 or 20"; or clumps (friends of friends within 20, 30 or 60") rotated and moved by up to 20" as rigid units.

python analysis/alignments_vasco_null.py -> data/alignments_vasco_null.json (100 copies per null, 20 with injected lines)"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial import cKDTree
import alignments as A                                      # noqa: E402
cat = A.cat; B, C = A.BINS, A.EDGES
kb0, kb1 = int(round((-2 - B[0]) / 0.1)), len(B) - 1                 # bulk: log10 mu in [-2, 0]
cs = C[:-1] + 0.05

def shift_all(df, rad, rng):
    rr = rad * np.sqrt(rng.uniform(0, 1, len(df))); a = rng.uniform(0, 2 * np.pi, len(df))
    return df.assign(x=df.x.values + rr * np.cos(a), y=df.y.values + rr * np.sin(a))

def clumps_rigid(df, link, rad, rng):
    """Friends-of-friends clumps (link arcsec) moved as rigid units (alignments.clump_copy), with the clump sizes."""
    lab, G = A.clumps(df, link)
    return A.clump_copy(df, int(rng.integers(1 << 31)), link, rad), np.bincount(lab, minlength=G)


def inject(seed, n=150):
    """Artificial lines like VASCO's candidates: 4-6 transients, 2-10' long, 1" scatter, members at least 10" apart."""
    rng = np.random.default_rng(seed); rows = []; plates = cat.plate.unique(); nid = 10_000_000
    for _ in range(n):
        g = cat[cat.plate == rng.choice(plates)]
        r = int(rng.integers(4, 7)); L = rng.uniform(120, 600); th = rng.uniform(0, np.pi)
        x0, y0 = rng.uniform(-8000, 8000, 2)
        s = np.sort(rng.uniform(0, L, r))
        if np.min(np.diff(s)) < 10: continue
        perp = rng.normal(0, 1.0, r); along = rng.normal(0, 1.0, r)
        rows.append(pd.DataFrame({"plate": g.plate.iloc[0], "x": x0 + (s + along) * np.cos(th) - perp * np.sin(th),
                                  "y": y0 + (s + along) * np.sin(th) + perp * np.cos(th)}, index=range(nid, nid + r))); nid += r
    return pd.concat([cat[["plate", "x", "y"]]] + rows)

NULLS = {"shift 3": ("s", 3), "shift 5": ("s", 5), "shift 10": ("s", 10), "shift 20": ("s", 20),
         "clumps 20+20": ("c", 20), "clumps 30+20": ("c", 30), "clumps 60+20": ("c", 60)}

def make(df, name, rng):
    kind, v = NULLS[name]
    return shift_all(df, v, rng) if kind == "s" else clumps_rigid(df, v, 20, rng)[0]

def summarise(H2, keep):
    post = H2[:, 1]                                                  # (size, mu, chi2) passing VASCO's screen
    bulk = post[:, kb0:kb1]
    straight = bulk[:, :, cs <= np.log10(2)].sum((1, 2)); loose = bulk[:, :, cs > np.log10(2)].sum((1, 2))
    L3 = A.arbitrated(keep, 1e-3)
    inj = sum(any(t >= 10_000_000 for t in a["ids"]) for a in L3)
    return np.r_[straight, loose, len(L3), sum(a["chi2ndof"] <= 2 for a in L3), inj]

def job(args):
    name, seed, with_inj = args
    rng = np.random.default_rng(seed)
    df = inject(77) if with_inj else cat
    _, _, H2, _, _, keep, _ = A.score(df if name == "real" else make(df, name, rng))
    return name, with_inj, summarise(H2, keep)

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for link in (20, 30, 60):
        _, sz = clumps_rigid(cat, link, 20, rng)
        print(f"link {link}\": groups of 2+ {np.sum(sz > 1)}, transients in them {sz[sz > 1].sum()}, largest {sz.max()}")
    jobs = [("real", 0, False), ("real", 0, True)] + [(n, 100 + s, False) for n in NULLS for s in range(100)] + [(n, 500 + s, True) for n in NULLS for s in range(20)]
    with Pool(8) as pool: R = pool.map(job, jobs)
    real = [v for n, w, v in R if n == "real" and not w][0]; real_inj = [v for n, w, v in R if n == "real" and w][0]
    print(f"\nreal: bulk straight by size {real[:4]}, loose {real[4:8]}; lines<=1e-3 {real[8]}, straight {real[9]}")
    print(f"real + 150 injected lines: lines<=1e-3 containing injected transients {real_inj[10]}")
    print(f"\n{'null':14s} {'bulk straight: real/null (z)':44s} {'bulk loose: real/null (z)':44s} {'lines 1e-3':12s} {'straight':10s} {'injected surviving':18s}")
    res = {"real": real.tolist(), "real_injected_found": int(real_inj[10]), "nulls": {}}
    for n in NULLS:
        V = np.array([v for nn, w, v in R if nn == n and not w]); Vi = np.array([v for nn, w, v in R if nn == n and w])
        m, sd = V.mean(0), V.std(0)
        fz = lambda i: f"{real[i] / m[i]:4.2f}({(real[i] - m[i]) / sd[i]:+4.1f})"
        print(f"{n:14s} {' '.join(fz(i) for i in range(4)):44s} {' '.join(fz(i) for i in range(4, 8)):44s} "
              f"{m[8]:5.1f} p{np.mean(V[:, 8] >= real[8]):.3f}  {m[9]:4.1f} p{np.mean(V[:, 9] >= real[9]):.3f}  {Vi[:, 10].mean():5.1f} of {real_inj[10]}")
        res["nulls"][n] = {"mean": m.round(2).tolist(), "sd": sd.round(2).tolist(), "injected_surviving": round(float(Vi[:, 10].mean()), 1),
                           "p_lines": round(float(np.mean(V[:, 8] >= real[8])), 3), "p_straight": round(float(np.mean(V[:, 9] >= real[9])), 3)}
    (Path(__file__).resolve().parents[1] / "data/alignments_vasco_null.json").write_text(json.dumps(res))
