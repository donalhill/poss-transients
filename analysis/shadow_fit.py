"""Earth's shadow, measured plate by plate on VASCO's catalogue.

Hypothesis: on each plate, transients are r times as frequent while the place they appear is in Earth's shadow as while
it is in sunlight (r = 1: no effect). The shadow follows shadow.py's model: a transient is a flash at one unknown instant
of the exposure, and the place it appears is in shadow for a share s of the exposure (seen from Palomar the shadow slides
across a plate during an exposure; its edge follows the share of the Sun's disc the Earth hides). So a place with share s
has rate b (1 - (1 - r) s), b the rate with no effect, and each transient has its own share s_i, computed at its place.
Each transient counts by VASCO-ML's rating of how likely it is to be real (sample.py), in every share, density and toy
below; the likelihood uses each plate's effective number of transients (see measure). Also the fit in each tenth of the
ratings (VASCO-ML's decile edges, sample.decile), each transient counted once, as their deciles count them; and every
band at once, with a free signal purity in each (joint_bands).
Fit: given a plate's transients and F, the share of them expected in shadow with no effect (the plate's share in shadow
integrated over its face, weighted by the average density of transients there on plates no shadow reaches), the
likelihood of where they lie is
  ln L(r) = sum_i ln(1 - (1 - r) s_i) - n ln(1 - (1 - r) F),
which with every s_i 0 or 1 is the binomial count in shadow against sunlight. Conditioning on n removes how rich or poor
each plate is. One common r is fitted by maximum likelihood over every plate the shadow covers for part of its face or
part of the exposure, each weighing as much as its shares make it (a plate in shadow throughout carries nothing, and is
set aside); each plate's own r_i is shown with its 68% likelihood interval. Every share, in the data and in every toy and
check, is the model's own at the transient's place.

The likelihood treats transients as spread evenly and independently within a plate; they are not (clumps, busy edges,
holes). So the fit is calibrated with toys:
  - from the fitted model (parametric): transients placed at random over the plate's face as the model says, checking
    that the fit returns what went in, with the coverage it claims;
  - from the data: each real shadow (its exact shares, place, and sweep on the plate) laid on another catalogue plate no
    shadow reaches, with a similar number of transients; the toy's transients are that plate's real ones, each with the
    share of the shadow at its place. An effect r < 1 is injected by removing each transient with chance (1 - r) s, r > 1
    by keeping each with chance (1 + (r - 1) s) / r. Fitted exactly like the data: the no-effect distribution, and a
    Neyman interval on r.
Cross-checks: host plates chosen differently (any plate; the nearest on the sky); the real shadow moved at random over its
own plate (turned, flipped, and shifted, similar coverage) with the transients fixed; each plate left out in turn; and the
shadow model's systematics (shadow.variants). The shadow is seen from Palomar (shadow.MAIN); every catalogue plate on its
corrected date and time (sample.ours); the central 6 x 6 deg of each plate.
python analysis/shadow_fit.py -> data/shadow_fit.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit, prange

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                # noqa: E402
import shadow as SH                                                          # noqa: E402
from shadow import HALF, ROOT, cell_counts, turn, unit                       # noqa: E402

rng = np.random.default_rng(11)
T, TBELT, TCLOSE = 20000, 5000, 4000
GRID = np.exp(np.linspace(np.log(0.05), np.log(20), 1201))
RBELT = np.round(np.arange(0.1, 2.501, 0.025), 3)
LAT = -2 * HALF + (2 * np.arange(int(4 * 30 * HALF)) + 1) / 60               # sub-point lattice, -6 to 6 deg (1/30 deg)
PLAT = -HALF + (2 * np.arange(int(2 * 30 * HALF)) + 1) / 60                  # the same over a plate's central square


def ell(s, u, grid=GRID):
    """u ln(1 - (1 - r) s) on the grid of r, for points with shares s and weights u (points x grid)."""
    return u[:, None] * np.log(1 - (1 - grid[None, :]) * s[:, None])


def curve(s, u, F, grid=GRID):
    """A plate's log-likelihood on the grid of r: sum_i u_i ln(1 - (1 - r) s_i) - n ln(1 - (1 - r) F)."""
    m = s > 0
    return ell(s[m], u[m], grid).sum(0) - u.sum() * np.log(1 - (1 - grid) * F)


def segment_sums(E, starts, n):
    """Row sums of E over n consecutive segments beginning at starts (an empty segment sums to zero)."""
    out = np.zeros((n,) + E.shape[1:]); ne = np.flatnonzero(np.diff(np.r_[starts, len(E)]) > 0)
    if len(ne):
        out[ne] = np.add.reduceat(E, starts[ne], axis=0)
    return out


@njit(cache=True)
def _mix(x):
    """splitmix64: a well-mixed 64-bit integer from x."""
    x = x + np.uint64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


@njit(cache=True)
def _uniform(seed, t, j, i):
    """A uniform number in [0, 1) fixed by (seed, toy, plate, transient): the same whatever order threads run in."""
    x = _mix(_mix(_mix(np.uint64(seed) + np.uint64(t)) + np.uint64(j)) + np.uint64(i))
    return np.float64(x >> np.uint64(11)) * (1.0 / 9007199254740992.0)


@njit(cache=True)
def _loglik(th, S, W, C, F):
    r = np.exp(th); L = 0.0
    for i in range(len(S)):
        L += W[i] * np.log(1.0 - (1.0 - r) * S[i])
    for j in range(len(F)):
        L -= C[j] * np.log(1.0 - (1.0 - r) * F[j])
    return L


@njit(cache=True)
def maximise(S, W, C, F):
    """The r of highest likelihood, sum_i W_i ln(1 - (1 - r) S_i) - sum_j C_j ln(1 - (1 - r) F_j), for r from 0.05 to 20:
    bracketed on 81 values of ln r, then golden-section search to 1e-8 in ln r."""
    lo, hi, n = np.log(0.05), np.log(20.0), 81
    step = (hi - lo) / (n - 1); bi, bL = 0, -np.inf
    for k in range(n):
        L = _loglik(lo + k * step, S, W, C, F)
        if L > bL:
            bi, bL = k, L
    a, b = lo + max(bi - 1, 0) * step, lo + min(bi + 1, n - 1) * step
    g = 0.6180339887498949; x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = _loglik(x1, S, W, C, F), _loglik(x2, S, W, C, F)
    while b - a > 1e-8:
        if f1 < f2:
            a, x1, f1 = x1, x2, f2; x2 = a + g * (b - a); f2 = _loglik(x2, S, W, C, F)
        else:
            b, x2, f2 = x2, x1, f1; x1 = b - g * (b - a); f1 = _loglik(x1, S, W, C, F)
    th = 0.5 * (a + b)
    return np.exp(th if _loglik(th, S, W, C, F) >= bL else lo + bi * step)


@njit(parallel=True, cache=True)
def _toy_fits(choice, fs, fw, olo, ohi, ok_, oF, r_true, seed):
    """Fitted r of each toy (rows of choice: the option taken on each plate). An option is a set of transients with
    shares fs and weights fw (fs[olo[g]:ohi[g]]), a scale k (ok_[g]) and a no-effect share F (oF[g]). An effect is put in
    by removing each transient with chance 1 - (1 - (1 - r) s) / max(1, r)."""
    T, J = choice.shape; out = np.empty(T)
    for t in prange(T):
        n = 0
        for j in range(J):
            n += ohi[choice[t, j]] - olo[choice[t, j]]
        S = np.empty(n); W = np.empty(n); C = np.empty(J); F = np.empty(J); m = 0
        for j in range(J):
            g = choice[t, j]; k = ok_[g]; nk = 0.0
            for i in range(olo[g], ohi[g]):
                s = fs[i]
                if r_true != 1.0:
                    q = (1.0 - (1.0 - r_true) * s) / max(1.0, r_true)
                    if q < 1.0 and _uniform(seed, t, j, i) >= q:
                        continue
                nk += fw[i]
                if s > 0:
                    S[m] = s; W[m] = k * fw[i]; m += 1
            C[j] = k * nk; F[j] = oF[g]
        out[t] = maximise(S[:m], W[:m], C, F)
    return out


def options(sets):
    """Flatten the options of every plate (sets[j]: a list of (shares, weights, k, F)) for _toy_fits; returns the flat
    arrays and each plate's first option's index."""
    flat = [o for opts in sets for o in opts]
    size = np.array([len(o[0]) for o in flat]); olo = np.r_[0, np.cumsum(size)[:-1]]
    return (np.concatenate([o[0] for o in flat]), np.concatenate([o[1] for o in flat]), olo, olo + size,
            np.array([o[2] for o in flat], float), np.array([o[3] for o in flat], float), np.r_[0, np.cumsum([len(o) for o in sets])[:-1]])


def toy_fits(opts, picks, r_true, size):
    """Fitted r of size toys: on each plate, an option drawn from picks[j] (indices among that plate's options)."""
    fs, fw, olo, ohi, ok_, oF, base = opts
    choice = np.stack([base[j] + rng.choice(pk, size) for j, pk in enumerate(picks)], 1).astype(np.int64)
    return _toy_fits(choice, fs, fw, olo, ohi, ok_, oF, float(r_true), int(rng.integers(2 ** 62)))


def best(L, interval=False, grid=GRID):
    """The r of highest likelihood for each curve (last axis on the grid); with interval, its 68% likelihood range."""
    b = grid[L.argmax(-1)]
    if not interval:
        return b
    ok = 2 * (L.max(-1, keepdims=True) - L) <= 1
    return b, np.where(ok, grid, np.inf).min(-1), np.where(ok, grid, -np.inf).max(-1)


def crossed_plates():
    """Plates some shadow reaches (either frame, at its widest, at any time of the exposure: shadow_pooled.py)."""
    return set(json.loads((ROOT / "data/shadow_pooled.json").read_text())["hosts"]["crossed"])


def toys(P, pick, r_true, size):
    """Fitted r of data-driven toys: for each plate, a host drawn by pick(p) (indices into the host list) per toy, with
    the plate's shadow laid on it; an effect put in by removing the host's transients at random, each by its own share in
    shadow."""
    return toy_fits(P[0]["opts"], [pick(p) for p in P], r_true, size)


def summary(d):
    q = np.quantile(d, [0.05, 0.16, 0.5, 0.84, 0.95])
    return {k: round(float(v), 3) for k, v in zip(("q05", "q16", "median", "q84", "q95"), q)}


def neyman(P, pick, r_obs):
    """Central Neyman intervals (68%, 90%) on the true r from data-driven toys at each injected r."""
    rows = []
    for rt in RBELT:
        d = toys(P, pick, rt, TBELT)
        rows.append([rt, *np.quantile(d, [0.05, 0.16, 0.5, 0.84, 0.95])])
    B = np.array(rows)
    # each end where the toys' quantile crosses the observed r, interpolated between the two injected values either side
    # (the lower end at the first crossing of the upper quantile, the upper end at the last crossing of the lower one)

    def edge(q, first):
        s = np.flatnonzero(np.diff(np.sign(B[:, q] - r_obs)) != 0)
        if not len(s):                                                       # no crossing: the end of the injected range
            return float(B[0, 0] if first else B[-1, 0])
        i = s[0] if first else s[-1]
        return float(np.interp(r_obs, B[i:i + 2, q], B[i:i + 2, 0])) if B[i, q] < B[i + 1, q] else float(B[i, 0])
    return B, (edge(4, True), edge(2, False)), (edge(5, True), edge(1, False))


def fit_inputs(P):
    """maximise's inputs for the real plates: their transients' shares and weights (k u), k n and F per plate."""
    S = np.concatenate([p["s"][p["s"] > 0] for p in P]); W = np.concatenate([p["k"] * p["u"][p["s"] > 0] for p in P])
    return S, W, np.array([p["k"] * p["nw"] for p in P]), np.array([p["f"] for p in P])


def hosts_of(cat, crossed):
    """The plates no shadow reaches, each plate's transients as a group."""
    return {q: gq for q, gq in cat[~cat.plate.isin(crossed)].groupby("plate")}


def measure(cat, pl, weighted, full, S=SH.MAIN, belt=True):
    """The fit on the plates the shadow (model S) partly covers, each transient counted by its rating (weighted) or once.
    Shares, the density map and every toy use the same weights. The likelihood takes each plate's effective number of
    transients, k = sum(w) / sum(w^2) times its summed ratings, so counting statistics are those of the ratings; with no
    weights k = 1. full: also the cross-checks (other hosts, the shadow moved over its own plate, each plate left out).
    belt: the Neyman belt, intervals, limit, and injection (without: the fit and its no-effect toys)."""
    cat = cat[(cat.xi.abs() < HALF) & (cat.eta.abs() < HALF)].assign(u=lambda d: d.w if weighted else 1.0)
    crossed = crossed_plates(); hosts = hosts_of(cat, crossed)
    # acceptance: the average density of transients (by weight) over the face of a plate, each host plate counted equally;
    # the share in shadow is weighted by it, so edges and corners count for what they hold
    ACC = np.mean([cell_counts(gq.xi.values, gq.eta.values, gq.u.values) / gq.u.sum() for gq in hosts.values()], 0)
    ACC = ACC / ACC.mean()
    ACCs = np.repeat(np.repeat(ACC.reshape(len(PLAT) // 3, -1), 3, 0), 3, 1)  # at each sub-point (rows: eta, columns: xi)
    H = list(hosts); Hc = np.array([len(hosts[h]) for h in H]); Hn = np.array([hosts[h].u.sum() for h in H])
    HK = np.array([hosts[h].u.sum() / (hosts[h].u ** 2).sum() for h in H]); hb = np.r_[0, np.cumsum(Hc)]
    Hxyz = unit(pl.loc[H].ra_c.values, pl.loc[H].dec_c.values)
    hxi, heta, hu = (np.concatenate([hosts[h][c_].values for h in H]) for c_ in ("xi", "eta", "u"))

    P, whole = [], []
    for q in sorted(crossed & set(pl.index)):
        p = pl.loc[q]; sub = SH.sub_dark(p, S)
        if not sub.any():
            continue
        if (sub > 1 - 1e-9).all():                                           # in shadow throughout: nothing to compare within it
            whole.append(q); continue
        c = cat[cat.plate == q]
        if c.empty:
            continue
        F = float(sub.mean(1) @ ACC / ACC.sum()); s = SH.dark(c.xi.values, c.eta.values, p, S); u = c.u.values
        hs = SH.dark(hxi, heta, p, S)                                       # the shadow laid on every host's transients
        NF = np.log(1 - (1 - GRID) * F)
        cen = Hxyz @ unit(p.ra_c, p.dec_c).reshape(3)
        n_w = u.sum(); sim = np.flatnonzero((Hn > n_w / 1.5) & (Hn < n_w * 1.5))
        P.append(dict(plate=q, row=p, f=F, sub=sub, area=float(sub.mean()), n=len(c), O=float(s.sum()), nw=float(n_w), Ow=float(u @ s),
                      k=float(n_w / (u ** 2).sum()), s=s, u=u, xi=c.xi.values, eta=c.eta.values, NF=NF, hs=hs, similar=sim if len(sim) >= 20 else np.argsort(np.abs(Hn - n_w))[:20], near=np.argsort(-cen)[:12]))
    # every host as an option on every plate, for the toys
    opts = options([[(p["hs"][hb[h]:hb[h + 1]], hu[hb[h]:hb[h + 1]], HK[h], p["f"]) for h in range(len(H))] for p in P])
    for p in P:
        p["opts"] = opts; del p["hs"]
    data = lambda keep: maximise(*fit_inputs([P[j] for j in keep]))
    Ld = np.array([p["k"] * curve(p["s"], p["u"], p["f"]) for p in P])     # plates x grid, for likelihood ranges
    r_hat = float(data(range(len(P)))); _, c_lo, c_hi = (float(v) for v in best(Ld.sum(0), True))
    out = {"weighted": weighted, "plates": [], "whole": whole,
           "fit": {"r": round(r_hat, 3), "count_lo": round(c_lo, 3), "count_hi": round(c_hi, 3), "O": round(float(sum(p["Ow"] for p in P)), 1),
                   "n": round(float(sum(p["nw"] for p in P)), 1), "transients": int(sum(p["n"] for p in P))}}
    for j, p in enumerate(P):
        r1, lo, hi = (float(v) for v in best(Ld[j], True))
        if p["Ow"] == 0:
            r1, lo = 0.0, 0.0
        # the rated count in shadow expected at the plate's sunlit rate: the observed one over r_i (with none in shadow,
        # the sunlit count carried over by F)
        E = p["Ow"] / r1 if r1 > 0 else p["nw"] * p["f"] / (1 - p["f"])
        out["plates"].append({"plate": p["plate"], "shadow": round(p["area"], 3), "weighted": round(p["f"], 3), "n": p["n"], "O": round(p["O"], 1),
                              "nw": round(p["nw"], 1), "Ow": round(p["Ow"], 1), "E": round(E, 1), "r": round(r1, 3),
                              "lo": round(lo if lo > GRID[0] else 0.0, 3), "hi": None if hi >= GRID[-1] else round(hi, 3)})
    print(pd.DataFrame(out["plates"]).to_string(index=False))
    print(f"common r {r_hat:.3f}, counting-only 68% {c_lo:.3f}-{c_hi:.3f}")

    # 1. the fit on its own model: each plate's transients placed at random over its face (the no-effect density times
    # 1 - (1 - r) s at each sub-point, the fitted r), fitted back
    Lc = np.zeros((TCLOSE, len(GRID)))
    for p in P:
        sub = p["sub"].ravel(); w = np.repeat(ACC, sub.size // ACC.size) * (1 - (1 - r_hat) * sub); w = w / w.sum()
        zero, one = sub <= 0, sub >= 1; part = ~zero & ~one
        N = int(round(p["nw"] * p["k"])); ellp = np.log(1 - (1 - GRID)[None, :] * sub[part][:, None])
        for t0 in range(0, TCLOSE, 500):
            cnt = rng.multinomial(N, np.r_[w[zero].sum(), w[one].sum(), w[part]], size=min(500, TCLOSE - t0)).astype(float)
            Lc[t0:t0 + len(cnt)] += cnt[:, 1:2] * np.log(GRID)[None, :] + cnt[:, 2:] @ ellp - N * p["NF"][None, :]
    rb, lo_b, hi_b = best(Lc, True)
    out["closure"] = {"median": round(float(np.median(rb)), 3), "coverage68": round(float(np.mean((lo_b <= r_hat) & (hi_b >= r_hat))), 3),
                      **{k: v for k, v in summary(rb).items() if k in ("q16", "q84")}}
    print("closure on the model:", out["closure"])

    # 2. toys from the data: the shadow laid on other plates with a similar total
    pick = lambda p: p["similar"]
    d0 = toys(P, pick, 1.0, T)
    out["toys"] = {"trials": T, "null": summary(d0), "p_as_high": round(float(np.mean(d0 >= r_hat)), 3),
                   "p_as_low": round(float(np.mean(d0 <= r_hat)), 3), "hist_null": np.histogram(d0, np.linspace(0, 3, 61))[0].tolist()}
    print("toys from the data, no effect:", out["toys"]["null"], "share as high as the data", out["toys"]["p_as_high"])
    if not belt:
        return out
    B, i68, i90 = neyman(P, pick, r_hat)
    # a 95% upper limit on the share of transients that vanish in shadow (1 - r): the lowest true r whose toys reach the
    # fitted r at least 5% of the time, from the belt's 95th percentile; and the limit expected with no effect, from the
    # no-effect toys in place of the data
    hi95 = np.maximum.accumulate(B[:, 5])
    up = lambda r_: 1 - np.interp(r_, hi95, B[:, 0])
    out["toys"].update({"limit": {"share95": round(float(up(r_hat)), 3), "expected": [round(float(np.quantile(up(d0), q)), 3) for q in (0.025, 0.16, 0.5, 0.84, 0.975)]},
                        "interval68": [round(v, 3) for v in i68], "interval90": [round(v, 3) for v in i90],
                        "belt": [[round(float(x), 3) for x in row] for row in B]})
    print("Neyman interval on r: 68%", i68, "90%", i90)
    # a shortfall as large as VASCO's (their pooled count over expected, on their construction in the same frame): how
    # often the fit would come back as high as it does on the real plates
    V_ = json.loads((ROOT / "data/shadow_pooled.json").read_text())[S.frame]["vasco"]
    r_v = round(V_["observed"] / V_["expected"], 2)
    d_v = toys(P, pick, r_v, T)
    out["injected"] = {"r_true": r_v, "p_as_high": round(float(np.mean(d_v >= r_hat)), 3),
                       "found": round(float(np.mean(d_v <= np.quantile(d0, 0.05))), 3),    # flagged at the 95% level
                       "hist": np.histogram(d_v, np.linspace(0, 3, 61))[0].tolist()}
    print("injected", out["injected"])
    if not full:
        return out

    # 3. cross-checks: other host choices; the shadow moved over its own plate; each plate left out
    out["checks"] = {}
    for key, pk in (("any_plate", lambda p: np.arange(len(H))), ("near_on_sky", lambda p: p["near"])):
        d = toys(P, pk, 1.0, T); _, a68, a90 = neyman(P, pk, r_hat)
        out["checks"][key] = {"null": summary(d), "p_as_high": round(float(np.mean(d >= r_hat)), 3),
                              "interval68": [round(v, 3) for v in a68], "interval90": [round(v, 3) for v in a90]}
        print(f"hosts {key}:", out["checks"][key])
    # same plate, transients fixed: the plate's own shadow turned, flipped, and shifted at random (in steps of the
    # sub-point lattice), its coverage within a factor 1.5 of the real one's; the shadow at each transient's place, and the
    # coverage, from the model at the moved places
    PX, PY = np.meshgrid(PLAT, PLAT); EX, EY = np.meshgrid(LAT, LAT); idx = lambda v: np.rint((v - LAT[0]) * 30).astype(int)
    same = []
    for p in P:
        E = SH.dark(EX, EY, p["row"], S)                                     # the shadow over a 12 x 12 deg square
        cand = np.array([(d, 0, 0) for d in range(8)] + list(zip(rng.integers(0, 8, 6000), *rng.integers(-len(PLAT) // 2, len(PLAT) // 2 + 1, (2, 6000)))))
        Fc = np.empty(len(cand))
        for i, (d, dx, dy) in enumerate(cand):
            X, Y = turn(PX - dx / 30, PY - dy / 30, d); Fc[i] = (E[idx(Y), idx(X)] * ACCs).sum() / ACCs.sum()
        ok = (Fc > p["f"] / 1.5) & (Fc < min(p["f"] * 1.5, 0.999))
        xy = [turn(p["xi"] - dx / 30, p["eta"] - dy / 30, d) for d, dx, dy in cand[ok]]
        so = SH.dark(np.concatenate([x for x, _ in xy]), np.concatenate([y for _, y in xy]), p["row"], S).reshape(ok.sum(), -1)
        same.append([(so[i], p["u"], p["k"], F_) for i, F_ in enumerate(Fc[ok])])
        # this plate on its own, with its shadow moved over it
        p["moved_r"] = np.array([maximise(so[i][so[i] > 0], p["k"] * p["u"][so[i] > 0], np.array([p["k"] * p["nw"]]), np.array([F_]))
                                 for i, F_ in enumerate(Fc[ok])])
    bs = toy_fits(options(same), [np.arange(len(o)) for o in same], 1.0, T)
    out["checks"]["same_plate"] = {"null": summary(bs), "p_as_high": round(float(np.mean(bs >= r_hat)), 3)}
    print("shadow moved over its own plate:", out["checks"]["same_plate"])
    out["checks"]["leave_one_out"] = {p["plate"]: round(float(data([i for i in range(len(P)) if i != j])), 3) for j, p in enumerate(P)}
    print("each plate left out:", out["checks"]["leave_one_out"])
    # the plate whose leaving out moves r most, on its own: how often a single plate gives a ratio as far out, with its
    # shadow laid on similar plates, on any plate, or moved over the plate itself
    loo = out["checks"]["leave_one_out"]; j = max(range(len(P)), key=lambda i: abs(loo[P[i]["plate"]] - r_hat)); p = P[j]
    r1 = float(data([j])); side = (lambda d: d <= r1) if r1 < 1 else (lambda d: d >= r1)
    fs, fw, olo, ohi, ok_, oF, base = P[0]["opts"]
    one = lambda pick: _toy_fits((base[j] + rng.choice(pick, T))[:, None].astype(np.int64), fs, fw, olo, ohi, ok_, oF, 1.0, int(rng.integers(2 ** 62)))
    out["checks"]["single"] = {"plate": p["plate"], "r": round(r1, 3), "p_similar": round(float(np.mean(side(one(p["similar"])))), 3),
                               "p_any": round(float(np.mean(side(one(np.arange(len(H)))))), 3),
                               "p_moved": round(float(np.mean(side(p["moved_r"]))), 3), "moved": int(len(p["moved_r"]))}
    print("the plate that moves r most, on its own:", out["checks"]["single"])
    return out


def joint_bands(cat, pl, band, S=SH.MAIN, size=T):
    """Signal purity in every band of the ratings at once. In band k a share pi_k of transients vanish in shadow (glints),
    the rest do not (flaws, missed stars), so its ratio is r_k = 1 - pi_k; pi_k is free in each band, as nothing says
    beforehand how purity varies with rating. Each band's transients on each plate partly covered follow the likelihood
    of measure, with that band's own density map (ratings vary across a plate). The test of no signal in any band, fixed
    before looking: q = sum over bands of 2 (ln L(pi hat) - ln L(0)), with each pi_k held at 0 or above (signal only
    removes transients). Its null: in each trial every shadow is laid on one host plate for all ten bands at once (a plate
    no shadow reaches, with a similar number of transients), so the bands keep how they move together."""
    inside = (cat.xi.abs() < HALF).values & (cat.eta.abs() < HALF).values
    cat = cat[inside].assign(k=band[inside])
    hosts = hosts_of(cat, crossed_plates()); H = list(hosts); NB = 10
    ACC = []
    for k in range(NB):
        m_ = [cell_counts(g.xi.values[g.k.values == k], g.eta.values[g.k.values == k]) / (g.k.values == k).sum() for g in hosts.values() if (g.k.values == k).any()]
        a = np.mean(m_, 0); ACC.append(a / a.mean())
    Hn = np.array([len(hosts[h]) for h in H]); hb = np.r_[0, np.cumsum(Hn)]
    hxi, heta, hk = (np.concatenate([hosts[h][c_].values for h in H]) for c_ in ("xi", "eta", "k"))
    R = np.linspace(0.0, 1.0, 401)[1:]                                                  # r = 1 - purity, from 1 (none) down
    Rw = np.exp(np.linspace(np.log(0.02), np.log(8), 601))                             # unconstrained r, for the trend
    one = np.eye(NB); hid_all = np.repeat(np.arange(len(H)), Hn)
    order = np.lexsort((hk, hid_all))                                          # transients by host, then band (stable)

    def band_terms(hs, key, nkb, F, grid):
        """Every host's curve in every band (hosts x bands x grid), from its transients' shares hs."""
        ks, ss = key[order], hs[order]; m = ss > 0
        st = np.r_[0, np.cumsum(np.bincount(ks[m], minlength=len(H) * NB))][:-1]
        Gp = segment_sums(ell(ss[m], np.ones(m.sum()), grid), st, len(H) * NB).reshape(len(H), NB, len(grid))
        return Gp - nkb[:, :, None] * np.log(1 - (1 - grid[None, None, :]) * F[None, :, None])

    def curves(s, kb, F, grid):
        """Per band, the log-likelihood curve (bands x grid) of transients with shares s and bands kb, no-effect shares F."""
        m = s > 0
        return one[kb[m]].T @ ell(s[m], np.ones(m.sum()), grid) - np.bincount(kb, minlength=NB)[:, None] * np.log(1 - (1 - grid[None, :]) * F[:, None])
    P = []
    for q in sorted(crossed_plates() & set(pl.index)):
        p = pl.loc[q]; sc = SH.cell_dark(p, S)
        if not sc.any() or (sc > 1 - 1e-9).all():
            continue
        c = cat[cat.plate == q]; s = SH.dark(c.xi.values, c.eta.values, p, S)
        F = np.array([float(sc @ a / a.sum()) for a in ACC]); hs = SH.dark(hxi, heta, p, S)
        s_ = np.flatnonzero((Hn > len(c) / 1.5) & (Hn < len(c) * 1.5))
        # each host's curves per band: its transients grouped by host then band, summed in segments
        key = hid_all * NB + hk; nkb = np.bincount(key, minlength=len(H) * NB).reshape(len(H), NB)
        G, Gw = (band_terms(hs, key, nkb, F, g_) for g_ in (R, Rw))
        P.append(dict(L=curves(s, c.k.values, F, R), Lw=curves(s, c.k.values, F, Rw), n=np.bincount(c.k.values, minlength=NB),
                      G=G, Gw=Gw, NF=np.log(1 - (1 - R[None, :]) * F[:, None]), hs=hs, sim=s_ if len(s_) >= 20 else np.argsort(np.abs(Hn - len(c)))[:20]))

    # benchmark D, VASCO's labels at face value: the signal spread over the bands as their labelled-real examples are, so
    # band k's purity is mu t_k / u_k (t: the labelled-real examples' share in each band; u: each band's share of the
    # transients on these plates), one strength mu for every band, up to the value at which a band's purity reaches 1
    lab = sample.labelled(); lr = lab[(lab.label == 1) & (lab.row >= 0)]
    n = np.array([p["n"] for p in P])
    t_ = np.bincount(sample.decile(lr.prob.values), minlength=NB) / len(lr); u_ = n.sum(0) / n.sum()
    MU = np.linspace(0, float(np.min(u_[t_ > 0] / t_[t_ > 0])), 201)
    RI = np.abs(R[None, None, :] - (1 - MU[:, None] * t_[None, :] / u_[None, :])[..., None]).argmin(2)          # mu x bands -> r index

    def q_stat(L):
        """L: trials x bands x R. Per band, 2 (ln L at the best r <= 1 - ln L at r = 1); the best r; and benchmark D's
        2 (ln L at the best mu - ln L at mu = 0), with the best mu."""
        LD = L[:, np.arange(NB)[None, :], RI].sum(2)                                                 # trials x mu
        return 2 * (L.max(2) - L[:, :, -1]), R[L.argmax(2)], 2 * (LD.max(1) - LD[:, 0]), MU[LD.argmax(1)]
    qb, rb, qd, mud = q_stat(sum(p["L"] for p in P)[None]); q_obs = float(qb.sum())

    def trials(purity, size, grid="R"):
        """Curves (trials x bands x grid) on trials from host plates, each band's transients removed with chance
        purity[k] times their share in shadow."""
        L = 0.0
        for p in P:
            h = rng.choice(p["sim"], size)
            if not np.any(purity):
                L = L + p["G" if grid == "R" else "Gw"][h]; continue
            Lp = np.zeros((size, NB, len(R)))
            for hh in np.unique(h):
                i = np.flatnonzero(h == hh); s, kb = p["hs"][hb[hh]:hb[hh + 1]], hk[hb[hh]:hb[hh + 1]]
                kept = rng.random((len(i), len(s))) < (1 - np.asarray(purity)[kb] * s)[None, :]; m = s > 0
                e = ell(s[m], np.ones(m.sum()), R)                                                  # transients x R
                for b in np.unique(kb[m]):
                    Lp[i, b] += kept[:, m][:, kb[m] == b] @ e[kb[m] == b]
                Lp[i] -= (kept.astype(float) @ one[kb])[:, :, None] * p["NF"][None]
            L = L + Lp
        return L

    def q_trials(purity, size):
        Q = []
        for t0 in range(0, size, 1000):
            q_ = q_stat(trials(purity, min(1000, size - t0))); Q.append(np.c_[q_[0].sum(1), q_[2]])
        return np.concatenate(Q)
    QQ = q_trials(np.zeros(NB), size); Q, QD = QQ[:, 0], QQ[:, 1]; q95 = float(np.quantile(Q, 0.95))
    # a trend with rating, as a signal among genuine transients would give: the slope of log r on each band's mean rating
    # (r unconstrained), against the same slope on trials from host plates, every band on the same hosts
    mr = np.array([cat.w[cat.k == k].mean() for k in range(NB)])
    line = lambda r_: np.polyfit(mr, np.log(r_), 1)                                   # slope, intercept
    s_obs = line(Rw[sum(p["Lw"] for p in P).argmax(1)])[0]; S_, Lt = [], []
    for t0 in range(0, size, 1000):
        for r_ in Rw[trials(np.zeros(NB), min(1000, size - t0), "Rw").argmax(2)]:
            b_, a_ = line(r_); S_.append(b_); Lt.append(a_ + b_ * mr)
    S_ = np.array(S_)
    # the uncertainty on a line fitted this way, slope and intercept together: each trial's line at every band's mean
    # rating, less the trials' median line (the band to draw about the fitted line, 68% and 90%)
    dev = np.array(Lt) - np.median(Lt, 0)
    # the test's reach: the top band's fitted purity in the top band alone, and VASCO's pooled shortfall in every band
    V_ = json.loads((ROOT / "data/shadow_pooled.json").read_text())[S.frame]["vasco"]; v_ = round(1 - V_["observed"] / V_["expected"], 2)
    inj = {"top_band": [0.0] * (NB - 1) + [round(float(1 - rb[0][-1]), 2)], "every_band": [v_] * NB}
    reach = {k: {"purity": v, "found": round(float(np.mean(q_trials(v, 5000)[:, 0] >= q95)), 3)} for k, v in inj.items()}
    return {"q": round(q_obs, 2), "p": round(float((np.sum(Q >= q_obs) + 1) / (len(Q) + 1)), 3), "trials": int(len(Q)),
            "q_by_band": [round(float(v), 2) for v in qb[0]], "purity": [round(float(1 - v), 3) for v in rb[0]],
            "null_q": {"median": round(float(np.median(Q)), 2), "q95": round(q95, 2)}, "reach": reach,
            "labels_template": {"template": [round(float(v), 3) for v in t_], "mu": round(float(mud[0]), 3), "q": round(float(qd[0]), 2),
                                "p": round(float((np.sum(QD >= qd[0]) + 1) / (len(QD) + 1)), 3), "top_band_purity_at_mu": round(float(mud[0] * t_[-1] / u_[-1]), 3)},
            "trend": {"slope": round(float(s_obs), 3), "null_median": round(float(np.median(S_)), 3),
                      "null_90": [round(float(np.quantile(S_, 0.05)), 3), round(float(np.quantile(S_, 0.95)), 3)],
                      "p": round(float(np.mean(S_ <= s_obs)), 3), "mean_rating": [round(float(v), 3) for v in mr],
                      "line_dev": {q: [round(float(v), 3) for v in np.quantile(dev, float(q), 0)] for q in ("0.05", "0.16", "0.84", "0.95")}}}


def task(args):
    """One piece of the run: the full fit, one rating band, every band at once, or one systematic variant, each with its
    own seed."""
    global rng
    kind, k = args
    rng = np.random.default_rng([11, ("main", "band", "joint", "variant").index(kind), k])
    cat, pl = sample.ours(); band = sample.decile(cat.w)
    if kind == "main":
        return measure(cat, pl, weighted=True, full=True)
    if kind == "band":
        edges = np.r_[0.0, sample.decile_edges(), 1.0]; m = measure(cat[band == k], pl, weighted=False, full=False)
        return {"band": k, "lo": round(float(edges[k]), 3), "hi": round(float(edges[k + 1]), 3), "fit": m["fit"],
                "plates": [{q: p_[q] for q in ("plate", "n", "O", "E", "r")} for p_ in m["plates"]],
                "toys": {q: m["toys"][q] for q in ("null", "p_as_low", "p_as_high", "interval68", "interval90", "limit")}}
    if kind == "joint":
        return joint_bands(cat, pl, band)
    name = list(SH.variants(SH.MAIN))[k]; m = measure(cat, pl, weighted=True, full=False, S=SH.variants(SH.MAIN)[name], belt=False)
    return {"r": m["fit"]["r"], "plates": len(m["plates"]), "whole": m["whole"], "p_as_low": m["toys"]["p_as_low"],
            "null_median": m["toys"]["null"]["median"]}


if __name__ == "__main__":
    import numba
    numba.set_num_threads(min(6, numba.config.NUMBA_NUM_THREADS))            # leave the machine some room
    jobs = [(kind, k) for kind, n in (("main", 1), ("joint", 1), ("band", 10), ("variant", len(SH.variants(SH.MAIN)))) for k in range(n)]
    res = {}
    for j in jobs:                                                             # one at a time
        res[j] = task(j); print("done", j, flush=True)
    out = res[("main", 0)]
    out["bands"] = [res[("band", k)] for k in range(10)]
    out["bands_joint"] = res[("joint", 0)]
    out["systematics"] = {name: res[("variant", k)] for k, name in enumerate(SH.variants(SH.MAIN))}
    print("r", out["fit"]["r"], "toys", out["toys"]["interval90"], "joint p", out["bands_joint"]["p"], flush=True)
    (ROOT / "data/shadow_fit.json").write_text(json.dumps(out))
