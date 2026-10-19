"""Do the transients in the strongest tracks look like images formed by the telescope? A probability for each transient,
and for each track, that its image is consistent with the stars imaged on the same plate.

For each transient, in a DSS cutout of its red POSS-I plate (tail_shapes.py fetches and fits the fields): the stars of
similar peak brightness (within a factor of two; three if there are too few) are cut out, centred on their fitted
positions and scaled to unit flux; the transient is centred on its own fitted position in the same way. A principal-component model of them (a linear autoencoder with K components, fitted to
that field's stars) describes what a star of that brightness looks like there, including any elongation or asymmetry.
The transient is put through the same model, and the misfit is the chi-squared of what is left over against the sky
noise. Each star gets the same treatment with itself left out of the model, so the transient's p-value is the share of
real stars that fit the model worse than it does. The p-values of a track's transients are combined by Fisher's method.
The same is done for transients picked at random across the catalogue, as a control.

python analysis/tail_psf.py -> data/tail_psf.json"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import shift as subshift
from scipy.stats import chi2, combine_pvalues

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tail_shapes as TS                                               # noqa: E402

ROOT = TS.ROOT
R_ = 6                                                                 # stamp half-width, pixels
K = 3                                                                  # principal components in the star model
MIN_STARS = 15
NCONTROL = 200                                                         # transients picked at random across the catalogue


def stamp(im, x, y):
    """A stamp centred on (x, y) to a fraction of a pixel, and its flux."""
    i, j = int(round(y)), int(round(x)); r = R_ + 2
    if i < r or j < r or i >= im.shape[0] - r or j >= im.shape[1] - r:
        return None, None
    s = subshift(im[i - r:i + r + 1, j - r:j + r + 1], (i - y, j - x), order=3, mode="nearest")[2:-2, 2:-2]
    f = s.sum()
    return (s / f, f) if f > 0 else (None, None)


def misfit(d, flux, M, sky):
    """Chi-squared of a unit-flux stamp against the star model M (mean, components), noise from the sky."""
    mean, comps = M
    r = d.ravel() - mean
    rec = mean + comps.T @ (comps @ r)
    return float(np.sum(((d.ravel() - rec) * flux / sky) ** 2)), rec.reshape(d.shape)


def model(S):
    X = np.array([s.ravel() for s in S]); mean = X.mean(0)
    comps = np.linalg.svd(X - mean, full_matrices=False)[2][:K]
    return mean, comps


def score(im, sky, t, stars):
    """The transient's p-value among the stars of similar brightness, each star scored with itself left out."""
    for fac in (2, 3):
        S = [s for s in stars if 1 / fac <= s["peak"] / t["peak"] <= fac]
        if len(S) >= MIN_STARS:
            break
    st = [stamp(im, s["xc"], s["yc"]) for s in S]; st = [(d, f) for d, f in st if d is not None]
    d_t, f_t = stamp(im, t["xc"], t["yc"])                            # both centred on their own fitted positions
    if d_t is None or len(st) < MIN_STARS:
        return None
    M = model([d for d, _ in st]); c_t, rec = misfit(d_t, f_t, M, sky)
    loo = [misfit(d, f, model([e for k, (e, _) in enumerate(st) if k != i]), sky)[0] for i, (d, f) in enumerate(st)]
    typical = st[int(np.argsort(loo)[len(loo) // 2])][0]                  # a star with a middling fit, for the page
    return {"chi2": c_t, "p": float((1 + np.sum(np.array(loo) >= c_t)) / (len(loo) + 1)), "n_stars": len(st),
            "stars_chi2": np.round(loo, 1).tolist(), "img": d_t, "rec": rec, "mean": M[0].reshape(d_t.shape), "star": typical}


def fields_of(q, c, size, sky_targets, name=None):
    fits_, stars, im, _ = TS.field(q, c, size, sky_targets, name)
    sky = 1.4826 * np.median(np.abs(im - np.median(im)))
    return fits_, stars, im, sky


if __name__ == "__main__":
    T = [a for a in json.loads((ROOT / "data/alignments.json").read_text())["track"]["candidates"] if a["score"] <= 1e-4]
    tracks, examples = [], []
    for a in T:
        q = a["plate"]; x, y = np.array(a["x"]) / 3600, np.array(a["y"]) / 3600
        ra_, dec_ = TS.SH.deproject(x, y, TS.pl.loc[q].ra_c, TS.pl.loc[q].dec_c); sk = SkyCoord(ra_ * u.deg, dec_ * u.deg)
        size = float(max(10, np.ceil(1.6 * a["d_arcmin"])))
        fits_, stars, im, sky = fields_of(q, SkyCoord(ra_.mean() * u.deg, dec_.mean() * u.deg), size, sk)
        rows = [score(im, sky, f, stars) if f else None for f in fits_]
        ps = [r["p"] for r in rows if r]
        comb = float(combine_pvalues(ps, method="fisher")[1]) if ps else None
        tracks.append({"plate": q, "p": [round(r["p"], 4) if r else None for r in rows], "chi2": [round(r["chi2"], 1) if r else None for r in rows],
                       "n_stars": [r["n_stars"] if r else None for r in rows], "combined_p": comb})
        print(f"{q}: transient p-values {[round(p, 3) for p in ps]}; combined (Fisher) {comb:.2g}")
        # the transient that fits worst, with its model and residual, for the page
        worst = min((r for r in rows if r), key=lambda r: r["p"])
        examples.append({"plate": q, "p": worst["p"], **{k: np.round(worst[k] * 100, 3).tolist() for k in ("img", "rec", "mean", "star")}})
    rng = np.random.default_rng(11); ctrl = []
    for k in rng.choice(len(TS.cat), NCONTROL, replace=False):
        r = TS.cat.iloc[k]; q = r.plate
        ra_, dec_ = TS.SH.deproject(np.array([r.xi]), np.array([r.eta]), TS.pl.loc[q].ra_c, TS.pl.loc[q].dec_c); sk = SkyCoord(ra_ * u.deg, dec_ * u.deg)
        try:
            fits_, stars, im, sky = fields_of(q, sk[0], 10.0, sk, f"control/{q}_{k}")       # its own cutout
        except (AssertionError, OSError, ValueError):
            continue
        if fits_[0]:
            s_ = score(im, sky, fits_[0], stars)
            if s_:
                ctrl.append({"plate": q, "p": round(s_["p"], 4), "chi2": round(s_["chi2"], 1)})
    pc = np.array([c_["p"] for c_ in ctrl])
    # each track against random sets of the same number of control transients, combined the same way
    rs = np.random.default_rng(3)
    for t in tracks:
        ps = [p for p in t["p"] if p is not None]
        sims = np.array([combine_pvalues(rs.choice(pc, len(ps)), method="fisher")[1] for _ in range(20000)])
        t["vs_transients"] = float(np.mean(sims <= t["combined_p"]))
        print(f"{t['plate']}: combined p {t['combined_p']:.3g}; random sets of {len(ps)} transients do as badly {t['vs_transients']:.4f}")
    print(f"control: {len(ctrl)} transients; share with p < 0.05: {np.mean(pc < 0.05):.2f}; < 0.01: {np.mean(pc < 0.01):.2f}; median p {np.median(pc):.2f}")
    (ROOT / "data/tail_psf.json").write_text(json.dumps({"tracks": tracks, "control": ctrl, "examples": examples, "k": K}))
