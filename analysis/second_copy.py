"""Are the strongest tracks, and the features around them, on the original plates? A second, independent scan.

The DSS scans used on the page are of one set of glass copies of the POSS-I plates (STScI). The SuperCOSMOS Sky Survey
(Royal Observatory Edinburgh) scanned its own set of copies, at 0.67" per pixel. For each of the strongest tracks
(alignments.py, track score 1e-4 or better), the SuperCOSMOS scan of the same exposure (matched by date and time) is
fetched, the DSS scan is resampled onto its pixels, and the two are compared within 2' of the track: the correlation of
their brightness after a flat sky level is removed, and the brightness of each transient on both. Something on only one
set of copies (copy damage, a scanner artefact) would appear on one scan only.

python analysis/second_copy.py -> data/second_copy.json"""
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import map_coordinates

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geometry import deproject                                                   # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = ROOT / "data/cutouts/second_copy"; CUT.mkdir(parents=True, exist_ok=True)
SIAP = "http://niuserre.roe.ac.uk:8080/ssa/SSS_SIAP"
pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")


def sss_wcs(h):
    w = fits.Header({k: h[k] for k in ("CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CTYPE1", "CTYPE2", "CD1_1", "CD1_2", "CD2_1", "CD2_2", "EQUINOX")})
    w["NAXIS"] = 2
    return WCS(w)


def flat(im):
    _, med, sd = sigma_clipped_stats(im[np.isfinite(im)])
    return (im - med) / sd


if __name__ == "__main__":
    T = [a for a in json.loads((ROOT / "data/alignments.json").read_text())["track"]["candidates"] if a["score"] <= 1e-4]
    out = []
    for a in T:
        q = a["plate"]; x, y = np.array(a["x"]) / 3600, np.array(a["y"]) / 3600
        ra, dec = deproject(x, y, pl.loc[q].ra_c, pl.loc[q].dec_c); sk = SkyCoord(ra * u.deg, dec * u.deg)
        c = SkyCoord(ra.mean() * u.deg, dec.mean() * u.deg); ut = Time(pd.Timestamp(pl.loc[q].ut).to_pydatetime())
        xml = subprocess.run(["curl", "-s", "-m", "120", f"{SIAP}?POS={c.ra.deg:.5f},{c.dec.deg:.5f}&SIZE=0.25"], capture_output=True).stdout
        rows = [r for r in parse_single_table(io.BytesIO(xml)).to_table() if str(r["bandpassID"]) == "POSSI_E"]
        r = min(rows, key=lambda r: abs(float(r["modified-Julian-date"]) - ut.mjd))          # the same exposure
        dt = abs(float(r["modified-Julian-date"]) - ut.mjd) * 24 * 60
        f = CUT / f"{q}_field{int(r['field'])}.fits.gz"
        if not f.exists():
            subprocess.run(["curl", "-s", "-m", "600", "-o", str(f), str(r["url"])], check=True)
        hs = fits.open(f)[0]; S_ = hs.data.astype(float); ws = sss_wcs(hs.header)
        dss = max(ROOT.glob(f"data/cutouts/tracks/{q}_poss1_red_*.fits"), key=lambda p: float(p.stem.split("_")[-1]))
        D_ = fits.getdata(dss).astype(float); wd = WCS(fits.getheader(dss))
        yy, xx = np.mgrid[0:S_.shape[0], 0:S_.shape[1]]; dx, dy = wd.world_to_pixel(ws.pixel_to_world(xx, yy))
        Don = map_coordinates(D_, [dy, dx], order=1, mode="constant", cval=np.nan)
        px, py = ws.world_to_pixel(sk)
        near = np.min(np.hypot(xx[..., None] - px, yy[..., None] - py), -1) < 120 / 0.67          # within 2' of the track
        A, B = flat(Don), flat(S_); m = near & np.isfinite(A)
        corr = float(np.corrcoef(A[m], B[m])[0, 1])
        peak = lambda im: [float(np.nanmax(im[int(round(j)) - 2:int(round(j)) + 3, int(round(i)) - 2:int(round(i)) + 3])) for i, j in zip(px, py)]
        out.append({"plate": q, "field": int(r["field"]), "minutes_apart": round(dt, 1), "correlation": round(corr, 3),
                    "peaks_dss": np.round(peak(A), 1).tolist(), "peaks_sss": np.round(peak(B), 1).tolist()})
        print(f"{q}: SuperCOSMOS field {int(r['field'])} ({dt:.0f} min from our time); correlation within 2' {corr:.2f}; "
              f"transient peaks (sky sd) DSS {np.round(peak(A), 0).tolist()} SuperCOSMOS {np.round(peak(B), 0).tolist()}")
    (ROOT / "data/second_copy.json").write_text(json.dumps(out))
