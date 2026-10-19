"""Do the transients in the strongest tracks look like images formed by the telescope? An image formed by the optics has
the shape of the stars around it on the same plate. A flaw in the emulsion need not (Hambly & Blair 2024 found
transients narrower than stars; Busko 2026 argues by eye from coma that transients on other plates are real images).

For each of the strongest tracks (track score 1e-4 or better, alignments.py), on a DSS cutout of the red POSS-I plate at
least 10' across: detect the stars (photutils), leaving out saturated ones and every catalogue transient; fit every star
and every transient in the track with the same elliptical Gaussian; and compare each transient's width with the stars of
similar peak brightness (within a factor of two) in the same cutout. Each transient and those stars are also cut out
centred on their fitted positions, less the local sky, and scaled to the same light within 4 pixels of the centre; the
stars' average is the star model, and pixels of a transient's stamp nearer another source are left out.

The spread of each image against the star model: every pixel's brightness at its distance from the fitted centre,
smoothed over distance (a Gaussian kernel half a pixel wide), scaled to the value at the centre, and divided by the star
model's profile scaled the same way. This is 1 at every distance for an image shaped like the stars, rises above 1 for
a wider one and falls below for a narrower one. For each track: the average over its transients, a band from resampling
its transients and redrawing their pixels within the sky noise, and the band the same number of real stars give (each
star against the model without it). A track is called wider or narrower when its average spread is beyond what those
stars give in 19 draws out of 20, either way (a flag for a closer look, not a detection). As a control, the same for
catalogue transients picked at random, each in its own cutout.

Also how far away a source would have to be for the telescope to blur it as wide as the least wide of the tracks
called wider: focused on the stars, the telescope spreads a point at distance L over a disc about APERTURE / L across,
taken here to add to the stars' width in quadrature.

python analysis/tail_shapes.py -> data/tail_shapes.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import astropy.units as u
from photutils.detection import DAOStarFinder
from scipy.ndimage import median_filter, shift as subshift

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geometry import deproject                                                   # noqa: E402
from stars import STAMP, fit                                                     # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = ROOT / "data/cutouts/tracks"; (CUT / "control").mkdir(exist_ok=True)
NCONTROL = 300                                                         # catalogue transients picked at random as a control
RG = np.arange(0, 2.31, 0.1)                                           # distances from the centre for the fade-out, pixels
H = 0.5                                                                # kernel width, pixels
NDRAW = 200                                                            # redraws of each transient within its sky noise
APERTURE = 1.22                                                        # metres, the Schmidt's 48-inch corrector plate
LEVEL = 0.05                                                           # either way, for a track to be called wider or narrower

cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet"); cat = cat[cat.plate != ""]
pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")


def dss(c, size, q, name=None):
    """A DSS cutout of the red plate, cached under the track's plate and size, or under its own name (one per random
    transient, so that two on the same plate do not share a cutout)."""
    f = CUT / f"{name or f'{q}_poss1_red_{size:g}'}.fits"
    if not f.exists():
        f.write_bytes(requests.get(f"https://archive.stsci.edu/cgi-bin/dss_search?v=poss1_red&r={c.ra.deg:.5f}&d={c.dec.deg:.5f}"
                                   f"&e=J2000&h={size}&w={size}&f=fits&c=none", timeout=300).content)
    h = fits.getheader(f); assert h["REGION"].strip() == q, (q, h["REGION"])
    return h, fits.getdata(f).astype(float)


def field(q, c, size, targets_sky, name=None):
    """Fit the targets and the stars in one cutout; return the targets' fits, the stars' fits, the image, and the positions
    of every source found and every catalogue transient (to keep neighbours out of a transient's stamp)."""
    h, R = dss(c, size, q, name); w = WCS(h)
    im = R - median_filter(R, size=31)
    _, med, sd = sigma_clipped_stats(im)
    found = DAOStarFinder(fwhm=2.5, threshold=5 * sd, sharplo=0.0, sharphi=3.0, roundlo=-3.0, roundhi=3.0)(im - med)   # no shape cuts on the stars
    tx, ty = w.world_to_pixel(targets_sky)
    # every catalogue transient in the cutout is left out of the stars
    g = cat[cat.plate == q]; ra_, dec_ = deproject(g.xi.values, g.eta.values, pl.loc[q].ra_c, pl.loc[q].dec_c)
    ax, ay = w.world_to_pixel(SkyCoord(ra_ * u.deg, dec_ * u.deg))
    sat = 0.8 * R.max() - np.median(R)
    stars = []
    xc, yc = ("x_centroid", "y_centroid") if found is not None and "x_centroid" in found.colnames else ("xcentroid", "ycentroid")
    for r in (found or []):
        x, y = float(r[xc]), float(r[yc])
        if np.min(np.hypot(ax - x, ay - y)) < 4 or im[int(round(y)), int(round(x))] > sat:
            continue
        f = fit(im, x, y)
        if f and f["fwhm"] < 10:
            stars.append({**f, "x": x, "y": y})
    out = []
    for x, y in zip(tx, ty):
        f = fit(im, x, y)
        out.append({**f, "x": float(x), "y": float(y)} if f else None)
    others = np.r_[np.c_[[float(r[xc]) for r in (found or [])], [float(r[yc]) for r in (found or [])]].reshape(-1, 2), np.c_[ax, ay]]
    return out, stars, im, others


RR = np.hypot(*np.mgrid[-STAMP:STAMP + 1, -STAMP:STAMP + 1])          # distance from the stamp's centre, pixels


def stamp(im, x, y):
    """A stamp centred on (x, y) to a fraction of a pixel, less the local sky (the median beyond 5 pixels)."""
    i, j = int(round(y)), int(round(x)); r = STAMP + 2
    if i < r or j < r or i >= im.shape[0] - r or j >= im.shape[1] - r:
        return None
    s = subshift(im[i - r:i + r + 1, j - r:j + r + 1], (i - y, j - x), order=3, mode="nearest")[2:-2, 2:-2]
    return s - np.median(s[RR >= 5])


def light(s):
    """The light within 4 pixels of the centre."""
    return None if s is None else float(s[RR <= 4].sum())


def images(im, t, stars, others):
    """The transient and the star model (the average of the stars of similar peak brightness), both scaled to unit light
    within 4 pixels; pixels of the transient's stamp nearer another source than the transient are left out (NaN)."""
    st = [stamp(im, s["xc"], s["yc"]) for s in stars if 0.5 <= s["peak"] / t["peak"] <= 2]
    st = [d / light(d) for d in st if d is not None and light(d) > 0]; d_t = stamp(im, t["xc"], t["yc"])
    if d_t is None or light(d_t) <= 0 or len(st) < 5:
        return None, None
    f = light(d_t); d_t = d_t / f
    yy, xx = np.mgrid[-STAMP:STAMP + 1, -STAMP:STAMP + 1]
    near = others[(np.hypot(*(others - [t["xc"], t["yc"]]).T) > 2) & (np.hypot(*(others - [t["xc"], t["yc"]]).T) < 2 * STAMP + 2)]
    for ox, oy in near:
        d_t[np.hypot(xx + t["xc"] - ox, yy + t["yc"] - oy) < RR] = np.nan
    return d_t, np.mean(st, 0)


def pixels(im, x, y):
    """The pixels within 6 of (x, y), unshifted: their distances from (x, y) and their values less the local sky (the
    median beyond 5 pixels), scaled to unit light within 4 pixels; the sky noise on that scale; and their positions."""
    i, j = int(round(y)), int(round(x))
    if i < STAMP or j < STAMP or i >= im.shape[0] - STAMP or j >= im.shape[1] - STAMP:
        return None
    yy, xx = np.mgrid[i - STAMP:i + STAMP + 1, j - STAMP:j + STAMP + 1]; r = np.hypot(xx - x, yy - y)
    v = im[i - STAMP:i + STAMP + 1, j - STAMP:j + STAMP + 1].astype(float); sky = v[r >= 5]
    v = v - np.median(sky); f = v[r <= 4].sum()
    return None if f <= 0 else (r.ravel(), v.ravel() / f, 1.4826 * np.median(np.abs(sky - np.median(sky))) / f, xx.ravel(), yy.ravel())


def kern(r, v):
    """Kernel sums of v over distance at each radius in RG: (sum of weight x value, sum of weights)."""
    ok = np.isfinite(v); w = np.exp(-0.5 * ((RG[:, None] - r[ok][None]) / H) ** 2)
    return w @ v[ok], w.sum(1)


def fall(num, den, m_num, m_den):
    """An image's fade-out, scaled to its centre, over the star model's scaled the same way."""
    p, m = num / den, m_num / m_den
    return (p / p[0]) / (m / m[0])


def curves(im, t, stars, others, rng):
    """The transient's fade-out against its star model: as measured, and NDRAW times with its pixels redrawn within the sky
    noise; and each star of similar peak brightness against the model without it. Pixels nearer another source than the
    transient are left out."""
    S = [pixels(im, s["xc"], s["yc"]) for s in stars if 0.5 <= s["peak"] / t["peak"] <= 2]; S = [q_ for q_ in S if q_]
    P = pixels(im, t["xc"], t["yc"])
    if P is None or len(S) < 5:
        return None
    r, v, sd, X, Y = P
    near = others[(np.hypot(*(others - [t["xc"], t["yc"]]).T) > 2) & (np.hypot(*(others - [t["xc"], t["yc"]]).T) < 2 * STAMP + 2)]
    for ox, oy in near:
        v[np.hypot(X - ox, Y - oy) < r] = np.nan
    ks = [kern(q_[0], q_[1]) for q_ in S]; Nm, Dm = np.sum([k[0] for k in ks], 0), np.sum([k[1] for k in ks], 0)
    ok = np.isfinite(v); w = np.exp(-0.5 * ((RG[:, None] - r[ok][None]) / H) ** 2); den = w.sum(1)
    draws = np.array([fall(w @ x, den, Nm, Dm) for x in v[ok] + rng.normal(0, sd, (NDRAW, ok.sum()))])
    return fall(w @ v[ok], den, Nm, Dm), draws, np.array([fall(n_, d_, Nm - n_, Dm - d_) for n_, d_ in ks])


def profile(items, rng, B=20000):
    """A track's average fade-out (geometric mean over its transients), its 68% band from resampling the transients and
    their noise redraws, and the 90% band of the same from one real star per transient. The track is called wider or
    narrower when its average fade-out over distance (the mean of its log) is beyond what those stars give in 19 draws
    out of 20, either way."""
    lg = lambda a: np.log(np.clip(a, 0.2, None))
    D = np.array([lg(d) for _, d, _ in items]); n = len(items)
    boot = np.exp(D[rng.integers(n, size=(B, n)), rng.integers(NDRAW, size=(B, n))].mean(1))
    null = np.mean([lg(sc[rng.integers(len(sc), size=B)]) for _, _, sc in items], 0)
    mid = np.mean([lg(c) for c, _, _ in items], 0); obs, nul = mid[1:].mean(), null[:, 1:].mean(1)
    p_w, p_n = float(np.mean(nul >= obs)), float(np.mean(nul <= obs))
    rd = lambda a: np.round(a, 3).tolist()
    return {"mid": rd(np.exp(mid)), "lo": rd(np.percentile(boot, 16, 0)), "hi": rd(np.percentile(boot, 84, 0)),
            "stars_lo": rd(np.exp(np.percentile(null, 5, 0))), "stars_hi": rd(np.exp(np.percentile(null, 95, 0))),
            "p_wider": p_w, "p_narrower": p_n, "call": "wider" if p_w < LEVEL else "narrower" if p_n < LEVEL else "like the stars"}


def compare(t, stars):
    """A transient against the stars of similar peak brightness in its field."""
    S = [s for s in stars if 0.5 <= s["peak"] / t["peak"] <= 2]
    if len(S) < 5:
        return None
    fw = np.array([s["fwhm"] for s in S]); el = np.array([s["elong"] for s in S])
    m = float(np.median(fw))
    return {**t, "n_stars": len(S), "fwhm_stars": m, "fwhm_pct": float(np.mean(fw < t["fwhm"]) * 100), "ratio": t["fwhm"] / m,
            "stars_p10": float(np.percentile(fw, 10) / m), "stars_p90": float(np.percentile(fw, 90) / m), "star_ratios": (fw / m).round(3).tolist(),
            "elong_stars": float(np.median(el)), "angle_stars": float(np.degrees(np.angle(np.mean(np.exp(2j * np.radians([s["angle"] for s in S]))))) / 2 % 180)}


if __name__ == "__main__":
    T = [a for a in json.loads((ROOT / "data/alignments.json").read_text())["track"]["candidates"] if a["score"] <= 1e-4]
    tracks = []; rp = np.random.default_rng(7)
    for a in T:
        q = a["plate"]; x, y = np.array(a["x"]) / 3600, np.array(a["y"]) / 3600
        ra_, dec_ = deproject(x, y, pl.loc[q].ra_c, pl.loc[q].dec_c); sk = SkyCoord(ra_ * u.deg, dec_ * u.deg)
        size = float(max(10, np.ceil(1.6 * a["d_arcmin"])))
        fits_, stars, im, others = field(q, SkyCoord(ra_.mean() * u.deg, dec_.mean() * u.deg), size, sk)
        rows = [compare(f, stars) if f else None for f in fits_]
        for r in rows:
            if r:
                d_t, m_ = images(im, r, stars, others)
                r.update(img=None if d_t is None else d_t.round(5).tolist(), model=None if m_ is None else m_.round(5).tolist())
        cv = [c for c in (curves(im, r, stars, others, rp) for r in rows if r) if c is not None]
        pr = profile(cv, rp)
        tracks.append({"plate": q, "n_stars": len(stars), "transients": rows, **{k: pr.pop(k) for k in ("call", "p_wider", "p_narrower")}, "profile": pr})
        print(q, tracks[-1]["call"], f"(p wider {tracks[-1]['p_wider']:.4f}, narrower {tracks[-1]['p_narrower']:.4f});", f"{len(stars)} stars;", " | ".join("no fit" if r is None else
              f"FWHM {r['fwhm']:.2f}px vs stars {r['fwhm_stars']:.2f} (percentile {r['fwhm_pct']:.0f}), elong {r['elong']:.2f} vs {r['elong_stars']:.2f}" for r in rows))
    # control: catalogue transients at random, each in a 10' cutout of its own plate
    rng = np.random.default_rng(11); ctrl = []; cv = []
    for k in rng.choice(len(cat), NCONTROL, replace=False):
        r = cat.iloc[k]; q = r.plate
        ra_, dec_ = deproject(np.array([r.xi]), np.array([r.eta]), pl.loc[q].ra_c, pl.loc[q].dec_c); sk = SkyCoord(ra_ * u.deg, dec_ * u.deg)
        try:
            fits_, stars, im, others = field(q, sk[0], 10.0, sk, f"control/{q}_{k}")
        except (AssertionError, OSError, ValueError):
            continue
        if fits_[0]:
            c_ = compare(fits_[0], stars)
            if c_:
                d_t, m_ = images(im, c_, stars, others)
                c2 = curves(im, c_, stars, others, rp)
                if d_t is not None and c2 is not None:
                    ctrl.append({"plate": q, **c_, "img": d_t.round(5).tolist(), "model": m_.round(5).tolist()}); cv.append(c2)
    pct = np.array([c_["fwhm_pct"] for c_ in ctrl])
    km = []
    for t in tracks:
        if t["call"] == "wider":
            ok = [r for r in t["transients"] if r]
            g = np.exp(np.mean([np.log(r["ratio"]) for r in ok])); ws = np.median([r["fwhm_stars"] for r in ok]) * 1.7
            km.append(APERTURE / np.radians(ws * np.sqrt(g ** 2 - 1) / 3600) / 1000)
            print(f"{t['plate']}: average width {g:.2f}x stars of {ws:.1f}\"; blur that wide from a source within {km[-1]:.0f} km")
    print(f"control: {len(ctrl)} catalogue transients; median width percentile among similar stars {np.median(pct):.0f}; "
          f"share narrower than 80% of their stars {np.mean(pct < 20):.2f}")
    (ROOT / "data/tail_shapes.json").write_text(json.dumps({"tracks": tracks, "control": ctrl, "control_profile": profile(cv, rp), "r_arcsec": np.round(RG * 1.7, 2).tolist(),
                                                                 "kernel_arcsec": H * 1.7, "pixel_arcsec": 1.7, "level": LEVEL,
                                                                 "defocus_km": float(max(km))}))
