"""The fraction of transients that vanish inside Earth's shadow (phi), from VASCO's plates, from the southern plates we
processed ourselves, from VASCO's nine plates re-measured with our extraction (shadow_ours.py), and from all of ours
together. Same estimator and interval construction as
shadow_phi_interval.py: phi_hat = 1 - sum O (1 - f) / sum (N - O) f over plates partly in shadow (0 < f < 0.97), where
O counts transients in the shadowed part, N on the plate, and f is the share expected in the shadowed part.
Sampling model, using nothing beyond the shadow plates themselves: each footprint is pasted on one of the shadow plates
of the same set (VASCO's nine for VASCO's catalogue; all of ours for our transients), in a random orientation, using
only sky where nothing could have vanished, and each transient inside is removed with probability phi. The spread of
transients across a plate (the template for the expected share f) also comes from the same plates.
python analysis/shadow_phi_combined.py [--fc]   -> data/phi_combined.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(11)
NREP = 20000
FC = "--fc" in sys.argv          # also the Feldman-Cousins scan over the vanishing fraction (slow); default: ratio only

D8 = [lambda m: m, np.rot90, lambda m: np.rot90(m, 2), lambda m: np.rot90(m, 3), np.fliplr, np.flipud, np.transpose,
      lambda m: np.rot90(m, 2).T]


def deproject(xi, eta, ra0, dec0):
    """Plate-frame offsets (deg) -> RA, Dec (deg), gnomonic, as in shadow_paste_null.py."""
    xi, eta, ra0, dec0 = map(np.radians, (xi, eta, ra0, dec0))
    rho = np.hypot(xi, eta); c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


# ---- VASCO's catalogue on its nine shadow plates (corrected dates, nights in VASCO's list), on the same kind of grid
def vasco_grid(half=3.6, step=0.05):
    sys.path.insert(0, str(ROOT / "analysis"))
    from astropy.time import Time
    from shadow import SHADOW_DEG, antisun_geocentric
    cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
    strict = set(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date.astype(str).str[:10])
    cat = cat[cat.obs_date.astype(str).str[:10].isin(strict)]
    pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
    sh = pd.read_parquet(ROOT / "data/plate_shadow_geo_corr.parquet")
    names = sh.index[sh.f_area > 0].tolist()
    g = np.arange(-half + step / 2, half, step); GX, GY = np.meshgrid(g, g); ng = len(g)
    r = np.maximum(np.abs(GX), np.abs(GY)); use = r < half
    masks, cells, n = [], [], []
    for q in names:
        p = pl.loc[q]
        ra, dec = deproject(GX.ravel(), GY.ravel(), p.ra_c, p.dec_c)
        masks.append((antisun_geocentric(ra, dec, Time(np.full(ra.size, pd.Timestamp(p.ut).to_datetime64()))) < SHADOW_DEG).reshape(ng, ng))
        c = cat[cat.plate == q]
        ix, iy = np.floor((c.xi.values + half) / step).astype(int), np.floor((c.eta.values + half) / step).astype(int)
        ok = (ix >= 0) & (ix < ng) & (iy >= 0) & (iy < ng)
        cells.append(np.ravel_multi_index((iy[ok], ix[ok]), (ng, ng))); n.append(int(ok.sum()))
    RB = np.linspace(0, half, 11)
    cnt = sum(np.histogram(r.ravel()[c], RB)[0] for c in cells)
    W = (cnt / np.histogram(r[use], RB)[0] / len(names))[np.clip(np.digitize(r, RB) - 1, 0, len(RB) - 2)]
    return dict(plates=np.array(names), sets=np.array(["vasco"] * len(names)), n=np.array(n), W=W,
                masks=np.array(masks), usable=np.array([use] * len(names)), cells=np.array(cells, dtype=object))


# ---- our own transients (shadow_ours.py): VASCO's nine shadow plates re-measured, and the 22 southern plates
OURS = None
if (ROOT / "data/shadow_ours_grid.npz").exists() and "--vasco-only" not in sys.argv:
    s = np.load(ROOT / "data/shadow_ours_grid.npz", allow_pickle=True)
    OURS = {k: s[k] for k in ("plates", "sets", "n", "W", "masks", "usable", "cells")}


def paste_set(G, which):
    """Observed (O, N, f) on the plates of G partly in shadow, and NREP trials with their footprints pasted on G's plates.
    Hosts lend only sky where nothing could have vanished: a host's own shadowed part is left out, while a plate wholly
    in shadow is kept whole, since any loss would thin it evenly and leave the spread of transients across it unchanged."""
    def oNf(mask, q, use):
        m, u = mask.ravel(), use.ravel()
        c = G["cells"][q][u[G["cells"][q]]]
        return m[c].sum(), len(c), (G["W"] * mask * use).sum() / (G["W"] * use).sum()

    P = len(G["plates"])
    obs = np.array([oNf(G["masks"][q], q, G["usable"][q]) for q in range(P)], dtype=float)
    host = [G["usable"][q] if obs[q, 2] >= 0.97 else G["usable"][q] & ~G["masks"][q] for q in range(P)]
    strad = np.flatnonzero((obs[:, 2] > 0) & (obs[:, 2] < 0.97) & np.isin(G["sets"], which))
    print(f"{'+'.join(which)}: plates partly in shadow {', '.join(G['plates'][strad])}")
    O_, N_, f_ = (np.empty((NREP, len(strad))) for _ in range(3))
    for r in range(NREP):
        for c, (k, q) in enumerate(zip(strad, rng.choice(P, len(strad), replace=False))):
            O_[r, c], N_[r, c], f_[r, c] = oNf(D8[rng.integers(8)](G["masks"][k]), q, host[q])
    return dict(O=O_, N=N_, f=f_, obs=tuple(obs[strad].T)), G["plates"][strad].tolist()


def ratio_summary(r, null, alt):
    """Observed ratio against its spread when nothing vanishes (null) and when 32% vanish (alt)."""
    bins = np.round(np.arange(0.2, 1.8001, 0.02), 3)
    h0, _ = np.histogram(np.clip(null, bins[0], bins[-1] - 1e-9), bins)
    h1, _ = np.histogram(np.clip(alt, bins[0], bins[-1] - 1e-9), bins)
    sd = float(np.std(null))
    from scipy.stats import norm
    p_null = float(np.mean(np.abs(null - np.median(null)) >= abs(r - np.median(null))))   # two-sided
    p_alt = float(np.mean(alt >= r))                                                       # one-sided, alt lies below
    return {"obs": round(float(r), 3), "sd": round(sd, 3), "null_median": round(float(np.median(null)), 3),
            "p_null": round(p_null, 4), "sigma_null": round(float(norm.isf(max(p_null, 1e-6) / 2)), 2),
            "p_alt": round(p_alt, 5), "sigma_alt": round(float(norm.isf(max(p_alt, 1 / len(alt)))), 2),
            "bins": bins.tolist(), "null_hist": (h0 / h0.sum()).round(5).tolist(), "alt_hist": (h1 / h1.sum()).round(5).tolist()}


def interval(sets, label):
    """Each set gives its own phi_hat per trial; sets are combined weighted by the inverse variance of their phi_hat
    when nothing vanishes (a noisy set would otherwise dilute a precise one)."""
    keep = np.all([np.sum((d["N"] - d["O"]) * d["f"], axis=1) > 0 for d in sets], axis=0)
    den = [np.sum((d["N"][keep] - d["O"][keep]) * d["f"][keep], axis=1) for d in sets]
    est = lambda d, Ok, dn: 1 - np.sum(Ok * (1 - d["f"][keep]), axis=1) / dn
    w = np.array([1 / np.var(est(d, d["O"][keep], dn)) for d, dn in zip(sets, den)]); w /= w.sum()
    x0 = float(sum(wk * (1 - np.sum(d["obs"][0] * (1 - d["obs"][2])) / np.sum((d["obs"][1] - d["obs"][0]) * d["obs"][2]))
                   for wk, d in zip(w, sets)))
    oO = np.concatenate([d["obs"][0] for d in sets])
    phis = np.round(np.arange(0, 1.0001, 0.0025), 4) if FC else np.array([0.0, 0.32])
    edges = np.linspace(-2.5, 1.0, 701)
    xc = 0.5 * (edges[1:] + edges[:-1])
    P = np.empty((len(phis), len(xc)))
    U = [[rng.random(d["O"][keep].shape) for d in sets] for _ in range(4)]
    keep_xs = {}
    for i, ph in enumerate(phis):
        xs = np.concatenate([sum(wk * est(d, binom.ppf(u, d["O"][keep].astype(int), 1 - ph), dn)
                                 for wk, d, dn, u in zip(w, sets, den, uu)) for uu in U])
        if ph in (0.0, 0.32):
            keep_xs[ph] = xs
        h, _ = np.histogram(np.clip(xs, edges[0], edges[-1] - 1e-9), edges)
        h = gaussian_filter1d(h.astype(float), 3.0, mode="constant") if ph < 1 else h.astype(float)
        P[i] = h / h.sum()
    if not FC:
        res = {"plates": int(len(oO)), "weights": [round(float(v), 3) for v in w], "raw_statistic": round(x0, 3),
               "ratio": ratio_summary(1 - x0, 1 - keep_xs[0.0], 1 - keep_xs[0.32])}
        rr = res["ratio"]
        print(f"{label:34s} {res['plates']:2d} plates  ratio {rr['obs']:.3f} +- {rr['sd']:.3f} (from 1: {rr['sigma_null']:.1f} sigma; "
              f"from a third vanishing: {rr['sigma_alt']:.1f} sigma)", flush=True)
        return res
    j0 = int(np.clip(np.searchsorted(edges, x0) - 1, 0, len(xc) - 1))
    best = P.max(axis=0)

    def fc(cl):
        accept = np.zeros_like(P, dtype=bool)
        for i in range(len(phis)):
            r = P[i] / np.where(best > 0, best, 1)
            order = np.argsort(-r)
            accept[i, order[:np.searchsorted(np.cumsum(P[i][order]), cl) + 1]] = True
        inside = phis[accept[:, j0]]
        return float(inside.min()), float(min(inside.max() + (phis[1] - phis[0]) / 2, 1))

    like = P[:, j0]
    post = like / np.trapezoid(like, phis)                    # flat prior on the physical range [0, 1]
    cdf = np.concatenate([[0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(phis))])
    # expected 90% FC upper limit if no transients vanished: the median limit over measurements drawn with phi = 0
    acc = np.zeros_like(P, dtype=bool)
    for i in range(len(phis)):
        order = np.argsort(-(P[i] / np.where(best > 0, best, 1)))
        acc[i, order[:np.searchsorted(np.cumsum(P[i][order]), 0.9) + 1]] = True
    up_of_x = np.array([phis[acc[:, j]].max() if acc[:, j].any() else 0 for j in range(len(xc))])
    expected = float(up_of_x[np.searchsorted(np.cumsum(P[0]), 0.5)])
    res = {"plates": int(len(oO)), "weights": [round(float(v), 3) for v in w], "best": float(phis[np.argmax(like)]),
           "fc90": round(fc(0.9)[1], 3), "fc95": round(fc(0.95)[1], 3),
           "bayes90": round(float(np.interp(0.9, cdf, phis)), 3), "expected_fc90_if_none": round(expected, 3),
           "p_as_low_if_0.32": round(float(P[np.argmin(abs(phis - 0.32))][: j0 + 1].sum()), 4),
           "raw_statistic": round(x0, 3),
           # the same measurement as a ratio: density in shadow / density in sunlight (1 = no effect)
           "ratio": ratio_summary(1 - x0, 1 - keep_xs[0.0], 1 - keep_xs[0.32]),
           "likelihood": {"phi": phis[::4].tolist(), "rel": [round(v, 4) for v in (like / like.max())[::4]]}}
    rr = res["ratio"]
    print(f"{label:24s} ratio {rr['obs']:.3f} +- {rr['sd']:.3f} (from 1: {rr['sigma_null']:.1f} sigma; "
          f"from a third vanishing: {rr['sigma_alt']:.1f} sigma)")
    print(f"{label:24s} {res['plates']:2d} plates  best fit {res['best']:.2f}  90% upper limit: FC {res['fc90']:.3f}, "
          f"Bayes {res['bayes90']:.3f}, expected if none vanish {res['expected_fc90_if_none']:.3f}  "
          f"P(this low | 0.32) {res['p_as_low_if_0.32']:.4f}")
    return res


vasco, pv = paste_set(vasco_grid(), ["vasco"])
out = {"vasco": interval([vasco], "VASCO's plates, VASCO's catalogue"), "plates": {"vasco": pv}}
if OURS is not None:
    for key, which, label in (("south", ["south"], "Southern plates, ours"), ("north", ["north"], "VASCO's plates, ours"),
                              ("combined", ["south", "north"], "All plates, ours")):
        if np.isin(OURS["sets"], which).any():
            d, out["plates"][key] = paste_set(OURS, which)
            out[key] = interval([d], label)
(ROOT / "data/phi_combined.json").write_text(json.dumps(out))
