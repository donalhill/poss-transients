"""Measure the true mid-exposure time of POSS-I red (XE) plates from asteroid trails.

For each plate: known numbered asteroids on the plate (SkyBoT), JPL Horizons astrometric ephemeris from Palomar (675)
at 1-min steps, a DSS cutout from the same plate, a local Gaia-based astrometric recalibration, then the trail centroid
projected onto the ephemeris track gives the mid-exposure UT. Compared with the header DATE-OBS, read as hours and
minutes as in VASCO's catalogue (06:82 becomes 07:22; the exposure in fact starts at 06:49, the field being hundredths of
an hour: see date_check.resolve), and with the time implied by OBSHA.

Usage: python analysis/asteroid_clock.py XE422 XE369 ...
"""
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.imcce import Skybot
from astroquery.jplhorizons import Horizons
from astroquery.vizier import Vizier
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, SourceCatalog

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parents[1]
HDR = ROOT / "data/dss_headers"
CUT = ROOT / "data/cutouts"
CUT.mkdir(parents=True, exist_ok=True); HDR.mkdir(parents=True, exist_ok=True)
DSS = "https://archive.stsci.edu/missions/dss/FitsArchive/XEsurvey/"
PALOMAR = EarthLocation.from_geodetic(lon=-116.863 * u.deg, lat=33.356 * u.deg, height=1706 * u.m)
CUT_ARCMIN = 12


def plate_header(plate):
    f = HDR / f"{plate.lower()}.hdr"
    if not f.exists():
        listing = requests.get(DSS, timeout=60).text
        fn = sorted(set(re.findall(rf"{plate.lower()}_[0-9a-z]{{4}}_00_0\d\.fits", listing)))[0]
        f.write_bytes(requests.get(DSS + fn, headers={"Range": "bytes=0-57599"}, timeout=60).content)
    h = f.read_bytes().decode("ascii", "replace")
    cards = {c[:8].strip(): c[10:].split("/")[0].strip().strip("'").strip() for c in (h[i:i + 80] for i in range(0, len(h), 80))}
    m = re.match(r"(\d{4}-\d\d-\d\d)T(\d\d):(\d\d):(\d\d)", cards["DATE-OBS"])
    t = pd.Timestamp(m.group(1)) + pd.Timedelta(hours=int(m.group(2)), minutes=int(m.group(3)), seconds=int(m.group(4)))
    ra, dec, ha, exp = (float(cards[k]) for k in ("PLATERA", "PLATEDEC", "OBSHA", "EXPOSURE"))
    lst_ha = (ra / 15 + ha) % 24
    lst_hdr = Time(t).sidereal_time("mean", longitude=PALOMAR.lon).hour
    t_ha = t + pd.Timedelta(minutes=((lst_ha - lst_hdr + 12) % 24 - 12) * 60 / 1.0027379)
    return dict(plate=plate, plateid=cards["PLATEID"], date_obs_raw=cards["DATE-OBS"], t_hdr=t, t_ha=t_ha,
                exp=exp, ha=ha, ra=ra, dec=dec)


def tangent(ra, dec, ra0, dec0):
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    c = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    return (np.degrees(np.cos(dec) * np.sin(ra - ra0) / c),
            np.degrees((np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / c))


def candidates(p, nmax=8):
    t = Time(p["t_hdr"].to_pydatetime(), scale="utc")
    try:
        df = Skybot.cone_search(SkyCoord(p["ra"] * u.deg, p["dec"] * u.deg), 2.9 * u.deg, t, location="675").to_pandas()
    except Exception:
        return pd.DataFrame()
    if "Number" not in df.columns:
        return pd.DataFrame()
    df = df[pd.to_numeric(df["Number"], errors="coerce").notna()]
    xi, eta = tangent(df.RA.values, df.DEC.values, p["ra"], p["dec"])
    df = df[(np.abs(xi) < 2.7) & (np.abs(eta) < 2.7) & df.V.between(12.0, 16.5) & (df.posunc < 2)]
    return df.sort_values("V").head(nmax)


def ephemeris(number, t0, hours=4):
    start = (t0 - pd.Timedelta(hours=hours)).strftime("%Y-%m-%d %H:%M")
    stop = (t0 + pd.Timedelta(hours=hours)).strftime("%Y-%m-%d %H:%M")
    for attempt in range(3):
        try:
            e = Horizons(id=str(int(number)), id_type="smallbody", location="675",
                         epochs={"start": start, "stop": stop, "step": "1m"}).ephemerides(quantities="1", extra_precision=True)
            break
        except Exception:
            time.sleep(3)
    else:
        raise RuntimeError(f"Horizons failed for {number}")
    return pd.DataFrame({"t": pd.to_datetime(Time(e["datetime_jd"], format="jd").iso), "ra": e["RA"], "dec": e["DEC"]})


def cutout(name, ra, dec):
    f = CUT / f"{name}.fits"
    if not f.exists():
        url = (f"https://archive.stsci.edu/cgi-bin/dss_search?v=poss1_red&r={ra:.6f}&d={dec:.6f}"
               f"&e=J2000&h={CUT_ARCMIN}&w={CUT_ARCMIN}&f=fits&c=none")
        for attempt in range(5):
            try:
                f.write_bytes(requests.get(url, timeout=120).content)
                break
            except requests.RequestException:
                time.sleep(5 * (attempt + 1))
    with fits.open(f) as hd:
        return hd[0].data.astype(float), hd[0].header


def gaia(ra, dec, epoch_year):
    v = Vizier(columns=["RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "Gmag"], column_filters={"Gmag": "<18.5"}, row_limit=-1)
    g = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=CUT_ARCMIN * 0.75 * u.arcmin, catalog="I/355/gaiadr3")[0].to_pandas()
    dt = epoch_year - 2016.0
    pmra, pmde = g.pmRA.fillna(0).values, g.pmDE.fillna(0).values
    g["ra"] = g.RA_ICRS + pmra * dt / 3.6e6 / np.cos(np.radians(g.DE_ICRS))
    g["dec"] = g.DE_ICRS + pmde * dt / 3.6e6
    return g


def measure(p, ast):
    eph = ephemeris(ast.Number, p["t_hdr"])
    e0 = eph.iloc[(eph.t - p["t_hdr"]).abs().argmin()]
    img, h = cutout(f"{p['plate']}_{int(ast.Number)}", e0.ra, e0.dec)
    if h.get("PLATEID", "").strip() != p["plateid"]:
        return None, "cutout from another plate"
    w = WCS(h)
    bkg = Background2D(img, (32, 32), filter_size=(3, 3), bkg_estimator=MedianBackground())
    data = img - bkg.background
    seg = detect_sources(data, 3 * bkg.background_rms, npixels=6)
    if seg is None:
        return None, "no sources"
    cat = SourceCatalog(data, seg).to_table(["x_centroid", "y_centroid", "elongation", "orientation", "semimajor_axis",
                                             "semiminor_axis", "segment_flux"]).to_pandas()
    cat = cat.rename(columns={"x_centroid": "xcentroid", "y_centroid": "ycentroid",
                              "semimajor_axis": "semimajor_sigma", "semiminor_axis": "semiminor_sigma"})
    # local astrometric recalibration: affine map from header-WCS pixel prediction of Gaia stars to measured centroids
    epoch = Time(p["t_hdr"].to_pydatetime()).jyear
    g = gaia(e0.ra, e0.dec, epoch)
    gx, gy = w.world_to_pixel(SkyCoord(g.ra.values * u.deg, g.dec.values * u.deg))
    stars = cat[cat.elongation < 1.4]
    A = None
    for tol in (5.0, 2.5, 1.5):
        if A is not None:
            gx2, gy2 = (np.c_[gx, gy, np.ones_like(gx)] @ A).T
        else:
            gx2, gy2 = gx, gy
        d = np.hypot(stars.xcentroid.values[:, None] - gx2[None], stars.ycentroid.values[:, None] - gy2[None])
        j = d.argmin(1)
        ok = d[np.arange(len(stars)), j] < tol
        if ok.sum() < 8:
            return None, f"too few Gaia matches ({ok.sum()})"
        X = np.c_[gx[j[ok]], gy[j[ok]], np.ones(ok.sum())]
        A, *_ = np.linalg.lstsq(X, stars[["xcentroid", "ycentroid"]].values[ok], rcond=None)
    res = stars[["xcentroid", "ycentroid"]].values[ok] - X @ A
    rms = np.sqrt(np.mean(np.sum(res**2, 1))) * 1.7
    # ephemeris track in calibrated pixel coordinates
    tx, ty = w.world_to_pixel(SkyCoord(eph.ra.values * u.deg, eph.dec.values * u.deg))
    tx, ty = (np.c_[tx, ty, np.ones_like(tx)] @ A).T
    # remove detections matching Gaia stars, keep those near the track
    gxa, gya = (np.c_[gx, gy, np.ones_like(gx)] @ A).T
    dstar = np.hypot(cat.xcentroid.values[:, None] - gxa[None], cat.ycentroid.values[:, None] - gya[None]).min(1)
    dtrk = np.hypot(cat.xcentroid.values[:, None] - tx[None], cat.ycentroid.values[:, None] - ty[None])
    k = dtrk.argmin(1)
    near = (dtrk[np.arange(len(cat)), k] < 3.0) & (dstar > 3.0)
    if near.sum() == 0:
        return None, f"no unmatched source on track (astrom rms {rms:.2f}\")"
    cand = cat[near].assign(k=k[near], off=dtrk[np.arange(len(cat)), k][near]).sort_values("segment_flux", ascending=False).iloc[0]
    t_mid = eph.t.iloc[int(cand.k)]
    rate = np.hypot(np.diff(tx), np.diff(ty)).mean() * 1.7  # arcsec per minute
    trail = np.sqrt(max(12 * (cand.semimajor_sigma**2 - cand.semiminor_sigma**2), 0)) * 1.7  # arcsec
    return dict(plate=p["plate"], number=int(ast.Number), name=ast.Name, V=float(ast.V), astrom_rms=rms, n_ref=int(ok.sum()),
                off_track=cand.off * 1.7, elong=cand.elongation, rate_as_min=rate, trail_as=trail,
                trail_min=trail / rate if rate > 0 else np.nan,
                t_mid=t_mid, mid_minus_hdr=(t_mid - p["t_hdr"]).total_seconds() / 60,
                mid_minus_ha=(t_mid - p["t_ha"]).total_seconds() / 60, exp=p["exp"]), "ok"


if __name__ == "__main__":
    out = ROOT / "data/asteroid_clock.csv"
    rows = []
    for plate in sys.argv[1:]:
        p = plate_header(plate)
        print(f"\n{plate} {p['plateid']} DATE-OBS {p['date_obs_raw']} -> {p['t_hdr']}  exp {p['exp']:.0f} min  "
              f"HA-implied {p['t_ha'].strftime('%H:%M')}")
        for _, a in candidates(p).iterrows():
            try:
                r, msg = measure(p, a)
            except Exception as ex:
                r, msg = None, f"error: {ex}"
            if r is None:
                print(f"   {int(a.Number):>7} {a.Name:<14} V={a.V:5.1f}  skipped: {msg}")
                continue
            rows.append(r)
            print(f"   {r['number']:>7} {r['name']:<14} V={r['V']:5.1f}  mid-hdr {r['mid_minus_hdr']:+6.0f} min  "
                  f"mid-HA {r['mid_minus_ha']:+6.0f}  trail {r['trail_min']:4.0f} min  off {r['off_track']:.1f}\"  "
                  f"elong {r['elong']:.1f}  astrom {r['astrom_rms']:.2f}\" ({r['n_ref']})")
    df = pd.DataFrame(rows)
    if out.exists() and out.stat().st_size > 0:
        df = pd.concat([pd.read_csv(out, parse_dates=["t_mid"]), df]).drop_duplicates(["plate", "number"], keep="last")
    df.to_csv(out, index=False)
