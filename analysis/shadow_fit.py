"""Earth's shadow: maximum-likelihood fit of r = (transient density in shadow) / (density in sunlight), plate by plate,
using only the plates the shadow crosses.

Each plate is gridded into 0.1 deg cells. The count in cell c is Poisson with mean
    lambda_c = A_c * T_c * exp(b . X_c) * (1 - s_c + r s_c)
  A_c   available sky: the share of the cell that the selection could use, exp(-rho pi (5")^2) from the local density
        rho of Gaia and Pan-STARRS sources (the "no catalogued source within 5 arcsec" requirement)
  T_c   radial template: how transients rise towards the plate edge, measured on the sunlit cells of plates of the
        same scan type (25-micron, 15-micron resampled, or VASCO's catalogue), never on shadowed cells
  X_c   per plate: a normalisation and a linear gradient across the plate (x, y). A free radial shape per plate would
        mimic the shadow edge, which is a gentle arc across a 6 deg plate.
  s_c   share of the cell inside Earth's shadow (geocentric, 8.5 deg about the anti-sun point, plate mid-exposure)
Background parameters are profiled out at each r. Plates are independent, so their -2 ln L curves add (common r); the
interval is where the sum rises by 1. Plates wholly in shadow or sunlight carry no information on r and are left out.
Error calibration: fake shadows (8.5 deg discs crossing the plate in random directions) on each plate's own sunlit sky,
real shadow masked, should give r = 1; the width of the pulls ln(r)/sigma rescales the curves. Closure: removing 32% of transients in the real
shadow should return r = 0.68.
Sets: VASCO's catalogue on its nine shadow plates; our extraction on the southern plates, on VASCO's nine, and both.
python analysis/shadow_fit.py   -> data/shadow_fit.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False     # the Sun's position needs no fresh Earth-rotation tables; downloading them stalls the run

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from scan_io import centre, pixel_arcsec, plate_header, plate_to_sky   # noqa: E402
from shadow import SHADOW_DEG, antisun_axis, unit                          # noqa: E402

rng = np.random.default_rng(3)
HALF, CELL = 2.9, 0.1
g = np.arange(-HALF + CELL / 2, HALF, CELL)
GX, GY = [a.ravel() for a in np.meshgrid(g, g)]
SUB = np.array([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]) * CELL / 3       # 3 x 3 points per cell for s
RHO = np.maximum(np.abs(GX), np.abs(GY)) / HALF
RB = np.linspace(0, 1, 11)                                    # rings for the radial template
RING = np.clip(np.digitize(RHO, RB) - 1, 0, len(RB) - 2)
X = np.c_[np.ones_like(GX), GX / HALF, GY / HALF]             # per plate: normalisation and a linear gradient
RGRID = np.exp(np.linspace(np.log(0.25), np.log(4), 161))


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


def deproject(xi, eta, ra0, dec0):
    """Plate-frame offsets (deg) -> RA, Dec (deg), gnomonic, as in shadow_paste_null.py."""
    xi, eta, ra0, dec0 = map(np.radians, (xi, eta, ra0, dec0))
    rho = np.hypot(xi, eta); c = np.arctan(rho)
    dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def available(cat_plate, to_sky, ut):
    """Share of each cell the selection could use, from the local catalogue density rho: exp(-rho pi (5")^2), the chance
    that a random point has no Gaia (moved to the plate epoch) or Pan-STARRS source within 5". Most Gaia stars are also
    in Pan-STARRS, so the denser of the two is used. Catalogue positions are put on the cell grid through a quadratic
    map fitted to the cell centres' tangent-plane positions."""
    ra0, dec0 = to_sky(0.0, 0.0)
    tp = lambda ra, dec: tangent(ra, dec, float(ra0), float(dec0))
    cxi, ceta = tp(*to_sky(GX, GY))
    Q = lambda a, b: np.c_[np.ones_like(a), a, b, a * a, a * b, b * b]
    mx = np.linalg.lstsq(Q(cxi, ceta), GX, rcond=None)[0]; my = np.linalg.lstsq(Q(cxi, ceta), GY, rcond=None)[0]
    dt = (pd.Timestamp("2016-01-01") - pd.Timestamp(ut)).days / 365.25
    dens = []
    for key in ("gaia", "ps1"):
        cols = ["ra", "dec", "pmra", "pmdec"] if key == "gaia" else ["ra", "dec"]
        c = pd.read_parquet(ROOT / f"data/scans/{key}_{cat_plate}.parquet", columns=cols)
        ra, dec = c.ra.values, c.dec.values
        if key == "gaia":
            ra = ra - np.nan_to_num(c.pmra.values) / 3.6e6 * dt / np.cos(np.radians(dec))
            dec = dec - np.nan_to_num(c.pmdec.values) / 3.6e6 * dt
        xi, eta = tp(ra, dec)
        x, y = Q(xi, eta) @ mx, Q(xi, eta) @ my
        ix, iy = np.floor((x + HALF) / CELL).astype(int), np.floor((y + HALF) / CELL).astype(int)
        ok = (ix >= 0) & (ix < len(g)) & (iy >= 0) & (iy < len(g))
        dens.append(np.bincount(iy[ok] * len(g) + ix[ok], minlength=GX.size) / CELL ** 2)
    return np.exp(-np.maximum(*dens) * np.pi * (5 / 3600) ** 2)


def tangent(ra, dec, ra0, dec0):
    """RA, Dec -> gnomonic offsets (deg) about (ra0, dec0)."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cc
    return np.degrees(xi), np.degrees(eta)


def shadow_share(sub_xyz, ut):
    ax = antisun_axis(Time(pd.Timestamp(ut).to_datetime64()))
    ax = np.asarray(ax).reshape(3)
    ang = np.degrees(np.arccos(np.clip(sub_xyz @ ax, -1, 1)))
    return (ang < SHADOW_DEG).reshape(-1, len(SUB)).mean(1)


def our_plate(plate, ut, s92=False):
    h = plate_header(plate); xc, yc = centre(plate); sc = pixel_arcsec(plate) / 3600
    to_sky = lambda x, y: plate_to_sky(h, xc + np.asarray(x) / sc, yc + np.asarray(y) / sc)
    sra, sdec = to_sky((GX[:, None] + SUB[:, 0]).ravel(), (GY[:, None] + SUB[:, 1]).ravel())
    s = pd.read_parquet(ROOT / f"data/scans/selected_{plate}.parquet")
    ix = np.floor(((s.x - xc) * sc + HALF) / CELL).astype(int); iy = np.floor(((s.y - yc) * sc + HALF) / CELL).astype(int)
    ok = (ix >= 0) & (ix < len(g)) & (iy >= 0) & (iy < len(g))
    n = np.bincount(iy[ok] * len(g) + ix[ok], minlength=GX.size)
    use = np.ones(GX.size, bool)
    if s92:                                              # Pan-STARRS covers it only north of -29.5
        use = to_sky(GX, GY)[1] > -29.5
    fine = h["NAXIS1"] == 13824                                # 15-micron scan resampled to 25 microns (23040 * 3/5)
    return dict(plate=plate, n=n, A=available(plate, to_sky, ut), use=use, sub=xyz(sra, sdec), ut=ut,
                type="ours-15um" if fine else "ours-25um")


def vasco_plate(plate, p, cat):
    to_sky = lambda x, y: deproject(np.asarray(x), np.asarray(y), p.ra_c, p.dec_c)
    sra, sdec = to_sky((GX[:, None] + SUB[:, 0]).ravel(), (GY[:, None] + SUB[:, 1]).ravel())
    c = cat[cat.plate == plate]
    ix = np.floor((c.xi.values + HALF) / CELL).astype(int); iy = np.floor((c.eta.values + HALF) / CELL).astype(int)
    ok = (ix >= 0) & (ix < len(g)) & (iy >= 0) & (iy < len(g))
    n = np.bincount(iy[ok] * len(g) + ix[ok], minlength=GX.size)
    return dict(plate=plate, n=n, A=available(plate, to_sky, p.ut), use=np.ones(GX.size, bool), sub=xyz(sra, sdec), ut=p.ut,
                type="vasco")


def template(plates):
    """Radial template from the sunlit cells of the given plates (one scan type): transients per unit of available sky
    in each ring, relative to the mean."""
    num, den = np.zeros(len(RB) - 1), np.zeros(len(RB) - 1)
    for P in plates:
        m = P["use"] & (shadow_share(P["sub"], P["ut"]) == 0)
        num += np.bincount(RING[m], P["n"][m], minlength=len(num)); den += np.bincount(RING[m], P["A"][m], minlength=len(den))
    t = num / np.where(den > 0, den, np.nan)
    t = np.where(np.isfinite(t) & (t > 0), t, np.nanmin(t[t > 0]))
    return (t / np.nanmean(t))[RING]


def fit_curve(n, A, s, use, T):
    """-2 ln L(r) on RGRID with the background (normalisation, gradient) profiled; T is the fixed radial template."""
    m = use & (A > 0)
    n, A, s, Xm = n[m], A[m] * T[m], s[m], X[m]
    beta = np.zeros(X.shape[1]); beta[0] = np.log(max(n.sum(), 1) / A.sum())
    out = np.empty(len(RGRID))
    for k, r in enumerate(RGRID):
        off = np.log(A * (1 - s + r * s))
        for _ in range(30):
            mu = np.exp(Xm @ beta + off)
            step = np.linalg.solve(Xm.T @ (Xm * mu[:, None]) + 1e-9 * np.eye(len(beta)), Xm.T @ (n - mu))
            beta += step
            if np.max(np.abs(step)) < 1e-7:
                break
        mu = np.exp(Xm @ beta + off)
        out[k] = 2 * np.sum(mu - n * np.log(np.where(mu > 0, mu, 1e-300)))
    return out - out.min()


def summarise(curve, scale=1.0):
    c = curve / scale ** 2
    k = int(np.argmin(c)); lr = np.log(RGRID)
    if 0 < k < len(c) - 1:                               # parabola through the minimum and its neighbours
        a, b, _ = np.polyfit(lr[k - 1:k + 2], c[k - 1:k + 2], 2); best = -b / (2 * a)
    else:
        best = lr[k]
    inside = lr[c <= 1]
    return dict(r=round(float(np.exp(best)), 3), lo=round(float(np.exp(inside.min())), 3), hi=round(float(np.exp(inside.max())), 3))


def informative(P, s):
    a = (P["A"] * P["use"]); f = (a * s).sum() / a.sum()
    return 0.02 < f < 0.97, float(f)


def fake_disc(P, k):
    """A fake shadow: an 8.5 deg disc whose edge crosses the plate, centred 6-14 deg from the plate centre in a random
    direction."""
    c = P["sub"].reshape(-1, len(SUB), 3).mean(1).mean(0); c /= np.linalg.norm(c)
    e1 = np.cross(c, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(c, e1)
    th, d = np.radians(rng.uniform(0, 360)), np.radians(rng.uniform(6, 14))
    ax = np.cos(d) * c + np.sin(d) * (np.cos(th) * e1 + np.sin(th) * e2)
    return (np.degrees(np.arccos(np.clip(P["sub"] @ ax, -1, 1))) < SHADOW_DEG).reshape(-1, len(SUB)).mean(1)


def pulls(P, T, nfake=40):
    """Controls: fake shadows on the plate's own sunlit sky (real shadow masked); ln(r_fake) / sigma for each."""
    use = P["use"] & ~(shadow_share(P["sub"], P["ut"]) > 0)
    out = []
    for k in range(nfake):
        sf = fake_disc(P, k)
        a = P["A"] * use; f = (a * sf).sum() / max(a.sum(), 1e-9)
        if not 0.05 < f < 0.95:
            continue
        sm = summarise(fit_curve(P["n"], P["A"], sf, use, T))
        lr = np.log(sm["r"]); sig = np.log(sm["hi"] / sm["r"]) if lr < 0 else np.log(sm["r"] / sm["lo"])
        if sig > 0:
            out.append(lr / sig)
    return out


def run_set(plates, label):
    rows, curves, all_pulls = [], [], []
    T = {ty: template([P for P in plates if P["type"] == ty]) for ty in {P["type"] for P in plates}}
    for P in plates:
        s = shadow_share(P["sub"], P["ut"])
        ok, f = informative(P, s)
        if not ok:
            continue
        c = fit_curve(P["n"], P["A"], s, P["use"], T[P["type"]])
        curves.append(c); all_pulls += pulls(P, T[P["type"]])
        rows.append(dict(plate=P["plate"], shadow=round(f, 3), n=int((P["n"] * P["use"]).sum()), **summarise(c)))
    tot = np.sum(curves, axis=0)
    pw = float(np.std(all_pulls)) if len(all_pulls) > 5 else 1.0
    k = max(pw, 1.0)
    res = dict(label=label, plates=rows, stat=summarise(tot), pull_width=round(pw, 2), n_pulls=len(all_pulls),
               calibrated=summarise(tot, k), curve=(tot / k ** 2).round(3).tolist(), rgrid=RGRID.round(4).tolist())
    cal = res["calibrated"]
    sig = (np.log(cal["hi"] / cal["r"]) if cal["r"] < 1 else np.log(cal["r"] / cal["lo"]))
    res["sigma_from_1"] = round(float(np.log(cal["r"]) / sig), 2)
    c068 = np.interp(np.log(0.68), np.log(RGRID), tot / k ** 2) - (tot / k ** 2).min()
    res["sigma_from_068"] = round(float(np.sqrt(max(c068, 0))), 2)
    print(f"\n{label}: {len(rows)} plates; r = {cal['r']:.3f} [{cal['lo']:.3f}, {cal['hi']:.3f}] (pull width {pw:.2f} "
          f"from {len(all_pulls)} time slides; uncalibrated [{res['stat']['lo']:.3f}, {res['stat']['hi']:.3f}]); "
          f"{res['sigma_from_1']:+.1f} sigma from 1, {res['sigma_from_068']:.1f} sigma from 0.68", flush=True)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)
    return res


def closure(plates, frac=0.32, ntoy=100):
    """Remove each transient in the real shadow with probability frac; the fit should return r = 1 - frac (times the
    real r, which the ratio of means corrects for)."""
    base, toys = [], []
    sets = [(P, shadow_share(P["sub"], P["ut"])) for P in plates]
    sets = [(P, s) for P, s in sets if informative(P, s)[0]]
    T = {ty: template([P for P in plates if P["type"] == ty]) for ty in {P["type"] for P in plates}}
    base = summarise(sum(fit_curve(P["n"], P["A"], s, P["use"], T[P["type"]]) for P, s in sets))["r"]
    for _ in range(ntoy):
        tot = 0
        for P, s in sets:
            n = rng.binomial(P["n"], 1 - frac * s)
            tot = tot + fit_curve(n, P["A"], s, P["use"], T[P["type"]])
        toys.append(summarise(tot)["r"])
    return dict(base=base, mean_toy=round(float(np.mean(toys)), 3), expected=round(base * (1 - frac), 3),
                sd_toy=round(float(np.std(toys)), 3))


if __name__ == "__main__":
    log = pd.read_csv(ROOT / "data/shadow_plates_log.csv").astype({"plate_number": str})
    south = pd.read_csv(ROOT / "data/southern_shadow_scans.csv").astype({"plate_number": str}).merge(
        log[["plate_number", "evening_date", "ra", "dec", "ut"]], on=["plate_number", "evening_date", "ra", "dec"])
    north = pd.read_csv(ROOT / "data/north_shadow_scans.csv")
    pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
    have = lambda q: all((ROOT / f"data/scans/{k}_{q}.parquet").exists() for k in ("gaia", "ps1"))
    S = [our_plate(p.plate, p.ut, p.plate_number == "S 92") for p in south.itertuples()
         if (ROOT / f"data/scans/selected_{p.plate}.parquet").exists()]
    N = [our_plate(q, pl.loc[q, "ut"]) for q in north.plate if (ROOT / f"data/scans/selected_{q}.parquet").exists()]
    cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
    strict = set(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date.astype(str).str[:10])
    cat = cat[cat.obs_date.astype(str).str[:10].isin(strict)]
    V = [vasco_plate(q, pl.loc[q], cat) for q in north.plate if have(q)]
    out = {}
    if V:
        out["vasco"] = run_set(V, "VASCO's catalogue, VASCO's nine plates")
    out["south"] = run_set(S, "Our transients, southern plates")
    if N:
        out["north"] = run_set(N, "Our transients, VASCO's nine plates")
    out["combined"] = run_set(S + N, "Our transients, all plates")
    out["closure"] = closure(S + N)
    print("\nclosure:", out["closure"])
    (ROOT / "data/shadow_fit.json").write_text(json.dumps(out))
