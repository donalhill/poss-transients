"""Where each known asteroid's trail lies along its track on the plate, giving the plate's middle of exposure.

For each plate of VASCO's catalogue and each asteroid dates/asteroid_list.py keeps on it (the brightest, V 12 to 18.5):
a box of BOX x BOX pixels from the plate's own scan in STScI's archive around the asteroid's position at the plate's
nominal middle (read with HTTP requests of byte ranges, a row each, so only the box is transferred), a local astrometric
recalibration against Gaia DR3 at the plate's epoch (an affine map from the header's WCS to the measured centroids of
the stars), then the detection on the track (within 3 pixels, and not a Gaia star) with the most flux: its centroid
projected onto the track's local segment gives the middle of the exposure, not rounded to the minute. The track is JPL
Horizons' astrometric position from Palomar (675) every minute over +-5 hours of the plate's nominal middle, wide enough
to hold the log's time where it differs from the header's. Each trail's offset from its track, elongation, and length
(from its second moments) are kept for the cleanliness rules of dates/plate_times.py. The tracks are fetched first, a
few at a time; scan headers, boxes, Gaia queries, and tracks are cached in data/cutouts/, and the measurements run
across processes.
usage: python dates/asteroid_clock.py [WORKERS] -> data/asteroid_clock.csv"""
import re
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parents[1]
CUT = ROOT / "data/cutouts"
BOX, PIX = 240, 1.7                                                           # box side (pixels); scan pixel (arcsec)
DSS = "https://archive.stsci.edu/missions/dss/FitsArchive/Palomar_Observatory/XEsurvey/"


def get(url, **kw):
    for attempt in range(5):
        try:
            r = requests.get(url, timeout=30, **kw); r.raise_for_status(); return r
        except requests.RequestException:
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"failed: {url}")


def scan_header(fn):
    """The plate scan's FITS header (cached), and where its data begin."""
    f = CUT / f"{fn}.hdr"
    if not f.exists():
        f.write_bytes(get(DSS + fn + ".fits", headers={"Range": "bytes=0-115199"}).content)
    raw = f.read_bytes(); end = next(i for i in range(0, len(raw), 80) if raw[i:i + 8] == b"END     ")
    return fits.Header.fromstring(raw[:end + 80].decode("ascii")), (end // 2880 + 1) * 2880


def box(fn, name, x0, y0):
    """The scan's pixels [y0, y0 + BOX) x [x0, x0 + BOX), from one request of byte ranges (a row each); cached."""
    f = CUT / f"{name}.npy"
    if not f.exists():
        h, off = scan_header(fn); nx, bp = h["NAXIS1"], abs(h["BITPIX"]) // 8
        ranges = [(off + ((y0 + i) * nx + x0) * bp, off + ((y0 + i) * nx + x0 + BOX) * bp - 1) for i in range(BOX)]
        parts = {}
        for k in range(0, BOX, 100):                                         # at most 100 ranges a request (the server's 8 kB
            r = get(DSS + fn + ".fits", headers={"Range": "bytes=" + ",".join(f"{a}-{b}" for a, b in ranges[k:k + 100])})   # header limit)
            for chunk in r.content.split(b"--" + r.headers["Content-Type"].split("boundary=")[1].encode())[1:-1]:
                head, body = chunk.split(b"\r\n\r\n", 1)
                a = int(re.search(rb"(?i)content-range: bytes (\d+)-", head).group(1)); parts[a] = body[:BOX * bp]
        img = np.frombuffer(b"".join(parts[a] for a, _ in ranges), dtype=f">i{bp}").reshape(BOX, BOX)
        np.save(f, (img * h.get("BSCALE", 1.0) + h.get("BZERO", 0.0)).astype(np.float32))
    return np.load(f).astype(float)


def gaia(name, ra, dec, epoch_year):
    """Gaia DR3 stars (G < 18.5) around the cutout, moved to the plate's epoch by their proper motions; cached."""
    f = CUT / f"{name}_gaia.parquet"
    if not f.exists():
        from astroquery.vizier import Vizier
        v = Vizier(columns=["RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "Gmag"], column_filters={"Gmag": "<18.5"}, row_limit=-1)
        for attempt in range(5):
            try:
                v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=BOX * PIX / 60 * 0.75 * u.arcmin, catalog="I/355/gaiadr3")[0] \
                    .to_pandas().to_parquet(f)
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
    g = pd.read_parquet(f); dt = epoch_year - 2016.0
    g["ra"] = g.RA_ICRS + g.pmRA.fillna(0).values * dt / 3.6e6 / np.cos(np.radians(g.DE_ICRS))
    g["dec"] = g.DE_ICRS + g.pmDE.fillna(0).values * dt / 3.6e6
    return g


def track(plate, number, mid):
    """JPL Horizons' astrometric track of asteroid `number` from Palomar, every minute over +-5 hours of `mid`; cached."""
    f = CUT / f"{plate}_{int(number)}_track.parquet"
    if not f.exists():
        r = get("https://ssd.jpl.nasa.gov/api/horizons.api", params={
            "format": "json", "COMMAND": f"'{int(number)};'", "OBJ_DATA": "NO", "MAKE_EPHEM": "YES", "EPHEM_TYPE": "OBSERVER",
            "CENTER": "'675'", "START_TIME": f"'{(mid - pd.Timedelta(hours=5)):%Y-%m-%d %H:%M}'",
            "STOP_TIME": f"'{(mid + pd.Timedelta(hours=5)):%Y-%m-%d %H:%M}'", "STEP_SIZE": "'1 m'", "QUANTITIES": "'1'",
            "ANG_FORMAT": "DEG", "EXTRA_PREC": "YES", "CSV_FORMAT": "YES", "TIME_TYPE": "UT"})
        txt = r.json()["result"]; rows = txt[txt.index("$$SOE") + 5:txt.index("$$EOE")].strip().splitlines()
        parts = [[c.strip() for c in row.split(",")] for row in rows]
        pd.DataFrame({"t": pd.to_datetime([q[0] for q in parts], format="%Y-%b-%d %H:%M"),
                      "ra": [float(q[3]) for q in parts], "dec": [float(q[4]) for q in parts]}).to_parquet(f)
    return pd.read_parquet(f)


def along(cx, cy, tx, ty, t):
    """The time at which the track (tx, ty at times t) passes closest to (cx, cy), on its local segment."""
    k = int(np.argmin(np.hypot(tx - cx, ty - cy)))
    best = (np.inf, t.iloc[k])
    for a in (k - 1, k):
        if 0 <= a < len(tx) - 1:
            dx, dy = tx[a + 1] - tx[a], ty[a + 1] - ty[a]; L2 = dx * dx + dy * dy
            f = np.clip(((cx - tx[a]) * dx + (cy - ty[a]) * dy) / L2, 0, 1) if L2 > 0 else 0.0
            d = np.hypot(tx[a] + f * dx - cx, ty[a] + f * dy - cy)
            if d < best[0]:
                best = (d, t.iloc[a] + (t.iloc[a + 1] - t.iloc[a]) * f)
    return best[1]


def measure(job):
    """One asteroid on one plate: (plate, scan file, nominal middle, number, V, track) -> the trail's record, or the reason
    none was measured. Pixel positions are the scan's own (its header's plate solution), less the box's corner."""
    from photutils.background import Background2D, MedianBackground
    from photutils.segmentation import SourceCatalog, detect_sources
    plate, fn, mid, number, V, trk = job
    base = dict(plate=plate, number=int(number), V=float(V))
    try:
        e0 = trk.iloc[(trk.t - mid).abs().argmin()]
        name = f"{plate}_{int(number)}"
        h, _ = scan_header(fn); W = WCS(h)
        cx, cy = W.world_to_pixel(SkyCoord(e0.ra * u.deg, e0.dec * u.deg))
        x0, y0 = int(round(float(cx))) - BOX // 2, int(round(float(cy))) - BOX // 2
        if not (0 <= x0 and 0 <= y0 and x0 + BOX <= h["NAXIS1"] and y0 + BOX <= h["NAXIS2"]):
            return {**base, "status": "off the scan"}
        img = box(fn, name, x0, y0)
        pix = lambda ra_, dec_: tuple(np.asarray(v) - o for v, o in zip(W.world_to_pixel(SkyCoord(ra_ * u.deg, dec_ * u.deg)), (x0, y0)))
        bkg = Background2D(img, (32, 32), filter_size=(3, 3), bkg_estimator=MedianBackground())
        data = img - bkg.background
        seg = detect_sources(data, 3 * bkg.background_rms, n_pixels=6)
        if seg is None:
            return {**base, "status": "no sources"}
        cat = SourceCatalog(data, seg).to_table(["x_centroid", "y_centroid", "elongation", "semimajor_axis", "semiminor_axis",
                                                 "segment_flux"]).to_pandas()
        g = gaia(name, e0.ra, e0.dec, Time(mid.to_pydatetime()).jyear)
        gx, gy = pix(g.ra.values, g.dec.values)
        stars = cat[cat.elongation < 1.4]
        A = None
        for tol in (5.0, 2.5, 1.5):                                         # the affine map, refined at tighter tolerances
            gx2, gy2 = (np.c_[gx, gy, np.ones_like(gx)] @ A).T if A is not None else (gx, gy)
            d = np.hypot(stars.x_centroid.values[:, None] - gx2[None], stars.y_centroid.values[:, None] - gy2[None])
            j = d.argmin(1); ok = d[np.arange(len(stars)), j] < tol
            if ok.sum() < 8:
                return {**base, "status": f"too few Gaia matches ({ok.sum()})"}
            X = np.c_[gx[j[ok]], gy[j[ok]], np.ones(ok.sum())]
            A, *_ = np.linalg.lstsq(X, stars[["x_centroid", "y_centroid"]].values[ok], rcond=None)
        rms = float(np.sqrt(np.mean(np.sum((stars[["x_centroid", "y_centroid"]].values[ok] - X @ A) ** 2, 1))) * PIX)
        tx, ty = pix(trk.ra.values, trk.dec.values)
        tx, ty = (np.c_[tx, ty, np.ones_like(tx)] @ A).T
        gxa, gya = (np.c_[gx, gy, np.ones_like(gx)] @ A).T
        dstar = np.hypot(cat.x_centroid.values[:, None] - gxa[None], cat.y_centroid.values[:, None] - gya[None]).min(1)
        dtrk = np.hypot(cat.x_centroid.values[:, None] - tx[None], cat.y_centroid.values[:, None] - ty[None])
        k = dtrk.argmin(1); off = dtrk[np.arange(len(cat)), k]
        near = (off < 3.0) & (dstar > 3.0)
        if not near.any():
            return {**base, "status": "nothing on the track", "astrom_rms": rms}
        i = np.flatnonzero(near)[np.argmax(cat.segment_flux.values[near])]
        c = cat.iloc[i]
        rate = float(np.hypot(np.diff(tx), np.diff(ty)).mean() * PIX)       # arcsec per minute
        trail = float(np.sqrt(max(12 * (c.semimajor_axis ** 2 - c.semiminor_axis ** 2), 0)) * PIX)
        return {**base, "status": "ok", "astrom_rms": rms, "n_ref": int(ok.sum()), "off_track": float(off[i] * PIX),
                "elong": float(c.elongation), "rate_as_min": rate, "trail_as": trail, "trail_min": trail / rate if rate > 0 else np.nan,
                "t_mid": along(c.x_centroid, c.y_centroid, tx, ty, trk.t)}
    except Exception as ex:                                                   # a failed download or fit: recorded, not fatal
        return {**base, "status": f"error: {type(ex).__name__}: {ex}"[:200]}


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    CUT.mkdir(parents=True, exist_ok=True)
    N = pd.read_csv(ROOT / "data/plate_times_nominal.csv", index_col=0, parse_dates=["start"])
    hdr = pd.read_parquet(ROOT / "data/xe_headers.parquet").drop_duplicates("plate").set_index("plate").file
    L = pd.read_parquet(ROOT / "data/asteroids_on_plates.parquet")
    L = L[L.plate.isin(hdr.index)]
    mid = N.start + pd.to_timedelta(N.exp / 2, unit="min")
    t0 = time.time()
    def safe_track(r):
        try:
            return track(r[0], r[1], mid[r[0]])
        except Exception as ex:                                              # recorded below as the trail's status
            return f"error: no Horizons track: {type(ex).__name__}: {ex}"[:200]
    with ThreadPoolExecutor(4) as pool:                                      # the tracks first, a few requests at a time
        tracks = dict(zip(zip(L.plate, L.number), pool.map(safe_track, zip(L.plate, L.number))))
    bad = {k: v for k, v in tracks.items() if isinstance(v, str)}
    print(f"{len(tracks) - len(bad):,} tracks from Horizons ({len(bad)} failed), {time.time() - t0:.0f} s", flush=True)
    jobs = [(r.plate, str(hdr[r.plate]), mid[r.plate], r.number, r.V, tracks[(r.plate, r.number)]) for r in L.itertuples()
            if (r.plate, r.number) not in bad]
    rows = [dict(plate=p_, number=int(n_), V=float(v_), status=bad[(p_, n_)]) for p_, n_, v_ in zip(L.plate, L.number, L.V)
            if (p_, n_) in bad]
    print(f"{len(jobs):,} trails to measure on {len({j[0] for j in jobs})} plates ({len(set(L.plate) - set(hdr.index))} plates "
          f"have no scan of their own to cut from)", flush=True)
    with ProcessPoolExecutor(workers) as pool:
        for i, r in enumerate(pool.map(measure, jobs, chunksize=4)):
            rows.append(r)
            if (i + 1) % 200 == 0:
                print(f"{i + 1:,} measured, {time.time() - t0:.0f} s", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "data/asteroid_clock.csv", index=False)
    print(out.status.str.split(":").str[0].value_counts().to_string(), f"\n{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
