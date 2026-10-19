"""The dotted trail on plate XE492, followed across the plate.

The strongest track on XE492 (alignments.py) sits on a dotted trail. Using the straight line through the track's
transients, DSS cutouts of the red plate 60' across are taken every 0.75 deg along it (kept only where the cutout comes
from XE492 itself), the sky is removed, and a strip is resampled along the line, 10' either side. The trail is then
followed outward from its clear middle in half-arcminute steps: in each step, the offset across the strip where the
summed brightness peaks, searching near the curve fitted so far. A straight line and a quadratic are fitted to the
offsets; anything moving through space (a meteor, or an object in orbit) crosses the plate in a straight line here.
The dots along the trail are then found, and their spacing measured over its length and in six sections.

python analysis/xe492_trail.py -> data/xe492_trail.npz (the strip near the trail) and data/xe492_trail.json (the traced
offsets, both fits, and the trail's extent)"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import map_coordinates, median_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import shadow as SH                                                    # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = ROOT / "data/cutouts/xe492"; CUT.mkdir(parents=True, exist_ok=True)
PIX = 1.7                                                              # arcsec, DSS pixel
pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet"); ra_c, dec_c = pl.loc["XE492"].ra_c, pl.loc["XE492"].dec_c
best = {}
for c in json.loads((ROOT / "data/alignments.json").read_text())["track"]["candidates"]:     # sorted best first
    best.setdefault(c["plate"], c)
P = np.c_[best["XE492"]["x"], best["XE492"]["y"]] / 3600; c0 = P.mean(0)
d = np.linalg.svd(P - c0)[2][0]; n = np.array([-d[1], d[0]])

# the line's extent inside the plate's searched area, and the strip along it
s_all = np.arange(-10, 10, 0.01); inside = np.all(np.abs(c0 + np.outer(s_all, d)) < 3.2, 1)
s_lo, s_hi = s_all[inside].min(), s_all[inside].max()
step = PIX / 3600; S = np.arange(s_lo, s_hi, step); O = np.arange(-600, 600.1, PIX) / 3600
strip = np.full((len(O), len(S)), np.nan); used = np.full(len(S), np.inf)
for sc in np.arange(s_lo + 0.4, s_hi, 0.75):
    xi, eta = c0 + sc * d; ra, dec = SH.deproject(np.array([xi]), np.array([eta]), ra_c, dec_c)
    f = CUT / f"s{sc:+.2f}.fits"
    if not f.exists():
        f.write_bytes(requests.get(f"https://archive.stsci.edu/cgi-bin/dss_search?v=poss1_red&r={ra[0]:.5f}&d={dec[0]:.5f}"
                                   f"&e=J2000&h=60&w=60&f=fits&c=none", timeout=600).content)
    h = fits.getheader(f)
    if h["REGION"].strip() != "XE492":                                 # DSS serves the neighbouring plate near the edge
        continue
    im = fits.getdata(f).astype(float); im = im - median_filter(im, size=31); w = WCS(h)
    sel = np.abs(S - sc) < used
    SS, OO = np.meshgrid(S[sel], O)
    rr, dd = SH.deproject((c0[0] + SS * d[0] + OO * n[0]).ravel(), (c0[1] + SS * d[1] + OO * n[1]).ravel(), ra_c, dec_c)
    px, py = w.world_to_pixel(SkyCoord(rr * u.deg, dd * u.deg))
    ok = (px > 0) & (px < im.shape[1] - 1) & (py > 0) & (py < im.shape[0] - 1)
    v = np.full(px.shape, np.nan); v[ok] = map_coordinates(im, [py[ok], px[ok]], order=1); v = v.reshape(SS.shape)
    keep = np.isfinite(v).mean(0) > 0.8; idx = np.flatnonzero(sel)[keep]
    strip[:, idx] = v[:, keep]; used[idx] = np.abs(S[idx] - sc)
sd = np.nanmedian(np.abs(strip - np.nanmedian(strip))) * 1.4826

# follow the trail outward from the transients
chunk = int(30 / PIX)


def ridge(k, guess, win):
    blk = strip[:, k:k + chunk]
    if blk.shape[1] < chunk or np.isnan(blk).mean() > 0.3:
        return None, None
    prof = np.nanmedian(blk, 1); m = np.abs(O * 3600 - guess) <= win
    j = np.flatnonzero(m)[np.nanargmax(prof[m])]
    return O[j] * 3600, prof[j] / sd * np.sqrt(chunk) / 1.25


pts = []
for direction in (+1, -1):
    k = int(np.argmin(np.abs(S))); guess, miss, run = 0.0, 0, []
    while 0 <= k < len(S) - chunk and miss < 6:
        o, snr = ridge(k, guess, 30 if len(run) > 3 else 60)
        if o is not None and snr > 6:
            run.append((S[k] + chunk * step / 2, o)); miss = 0
            arr = np.array(run[-12:])
            guess = np.polyval(np.polyfit(arr[:, 0], arr[:, 1], min(2, len(arr) - 1)), S[k] + direction * chunk * step) if len(run) >= 4 else o
        else:
            miss += 1
        k += direction * chunk
    pts += run
pts = np.array(sorted(set(pts)))
s_pts, o_pts = pts[:, 0], pts[:, 1]
lin = np.polyfit(s_pts, o_pts, 1); quad = np.polyfit(s_pts, o_pts, 2)
rms = lambda cf: float(np.std(o_pts - np.polyval(cf, s_pts)))
radius = 1 / abs(2 * np.polyfit(s_pts * 3600, o_pts, 2)[0]) / 3600
covered_to = max(float(f.stem[1:]) for f in CUT.glob("s*.fits") if fits.getheader(f)["REGION"].strip() == "XE492") + 0.5   # end of XE492's cutouts
out = {"start_arcmin": float((s_pts.min() - s_lo) * 60), "end_arcmin": float((s_pts.max() - s_lo) * 60), "length_arcmin": float(np.ptp(s_pts) * 60),
       "ends_at_coverage": bool(s_pts.max() > covered_to - 1 / 60), "coverage_end_arcmin": float((covered_to - s_lo) * 60), "line_length_deg": float(s_hi - s_lo),
       "rms_line": rms(lin), "rms_quad": rms(quad), "radius_deg": float(radius), "sag_arcsec": float(np.ptp(np.polyval(quad, s_pts))),
       "s": np.round((s_pts - s_lo) * 60, 3).tolist(), "offset": np.round(o_pts, 1).tolist(), "quad": quad.tolist(), "s_lo": float(s_lo)}
# the dots along the fitted curve: brightness sampled along it (+/-1 pixel across), peaks above 3 sigma at least two
# pixels apart, and their spacing over the whole trail and in six sections; a steadily flashing source would give the
# same spacing all the way along
on = (S >= s_pts.min()) & (S <= s_pts.max()); rows = np.interp(np.polyval(quad, S[on]), O * 3600, np.arange(len(O)))
cols = np.flatnonzero(on)
prof = np.nanmean([map_coordinates(np.nan_to_num(strip), [rows + k, cols], order=1) for k in (-1, 0, 1)], 0)
a_ = (S[on] - s_lo) * 60
pk = [i for i in range(1, len(prof) - 1) if prof[i] > 3 * sd and prof[i] >= prof[i - 1] and prof[i] >= prof[i + 1]]
pk = [i for k, i in enumerate(pk) if k == 0 or i - pk[k - 1] >= 2]
gaps = np.diff(a_[pk]) * 60
sections = []
for lo, hi in zip(np.linspace(a_.min(), a_.max(), 7)[:-1], np.linspace(a_.min(), a_.max(), 7)[1:]):
    g = np.diff([a_[i] for i in pk if lo <= a_[i] < hi]) * 60
    x_ = prof[(a_ >= lo) & (a_ < hi)]; x_ = x_ - x_.mean(); ac = np.correlate(x_, x_, "full")[len(x_) - 1:]; ac /= ac[0]
    lag = np.arange(len(ac)) * PIX; win = (lag > 6) & (lag < 40)
    sections.append({"from": round(float(lo), 1), "to": round(float(hi), 1), "dots": int(len(g) + 1), "spacing": round(float(np.median(g)), 1),
                     "repeat": round(float(lag[win][np.argmax(ac[win])]), 1), "correlation": round(float(ac[win].max()), 2)})
out.update({"dots": len(pk), "spacing_median": float(np.median(gaps)), "spacing_16_84": [float(np.percentile(gaps, 16)), float(np.percentile(gaps, 84))],
            "sections": sections})
print(f"{len(pk)} dots; spacing median {np.median(gaps):.1f}\" (16-84% {np.percentile(gaps, 16):.1f}-{np.percentile(gaps, 84):.1f}\"); sections:",
      "; ".join(f"{x['spacing']}\" (repeat {x['repeat']}\", r {x['correlation']})" for x in sections))
print(f"trail from {out['start_arcmin']:.1f}' to {out['end_arcmin']:.1f}' along the line ({out['length_arcmin']:.1f}'), "
      f"{'running to where our scans of XE492 end' if out['ends_at_coverage'] else 'ending on the plate'}; "
      f"rms about a straight line {out['rms_line']:.1f}\", about a quadratic {out['rms_quad']:.1f}\"; radius of curvature {radius:.1f} deg")
(ROOT / "data/xe492_trail.json").write_text(json.dumps(out))
# the strip near the trail, for the page: from 4' before the trail to its end, 2' either side of the line
i0, i1 = int((s_pts.min() - 4 / 60 - S[0]) / step), int((s_pts.max() + 0.5 / 60 - S[0]) / step); j = np.abs(O) <= 120 / 3600
np.savez_compressed(ROOT / "data/xe492_trail.npz", strip=strip[j][:, max(i0, 0):i1].astype(np.float32), sd=sd,
                    along0=float((S[max(i0, 0)] - s_lo) * 60), step_arcmin=PIX / 60, across=O[j] * 3600)
