"""Which known asteroids lie on each plate of VASCO's catalogue, computed from their orbits: the Minor Planet Center's
osculating elements (MPCORB) integrated back to the plate epochs with ASSIST (the DE440 planets, the 16 massive asteroids
of sb441-n16, and general relativity), seen from Palomar (675) with light time. Every numbered asteroid bright enough to
matter (H <= HMAX) is integrated once through every plate's nominal middle of exposure (data/plate_times_nominal.csv,
from the header), in batches across processes, orbits sharing an epoch integrated together, backwards in time. An
asteroid is kept for a plate when it lies within the plate's central 5.4 x 5.4 deg with V between VMIN and VMAX; the
NMAX brightest are kept.
This only chooses the asteroids: the positions need be right to an arcminute, not better, so the integrator's tolerance
is EPS (against its default of 1e-9: at most 9", median 0", on a test of 300 orbits, twice as fast). The tracks the
clock measures against come from JPL Horizons (dates/asteroid_clock.py): over 70 years, MPCORB's elements integrated
here differ from JPL's orbits by up to 2", minutes along a slow asteroid's track.
Inputs (data/external/ephemerides/, fetched by the Makefile): MPCORB.DAT.gz, de440.bsp, sb441-n16.bsp. Runs under a
python with ASSIST (conda-forge's assist and rebound; see the README).
usage: python dates/asteroid_list.py [WORKERS] -> data/asteroids_on_plates.parquet"""
import gzip
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EPH = ROOT / "data/external/ephemerides"
HMAX, VMIN, VMAX, NMAX, BOX, EPS = 17.5, 12.0, 18.5, 12, 2.7, 1e-4
K = 0.01720209895                       # Gauss, AU^1.5 / day
EPS = np.radians(84381.448 / 3600)      # J2000 obliquity, as MPC
C_AU_D = 173.144632674                  # speed of light, AU / day
LON, LAT, HGT = np.radians(-116.8650), np.radians(33.3563), 1.706   # Palomar (675), km
RE_AU = 6378.137 / 1.495978707e8


def unpack_epoch(s):
    """MPC packed date (e.g. K25BL) -> JD (TT)."""
    yr = {"I": 1800, "J": 1900, "K": 2000}[s[0]] + int(s[1:3])
    dig = lambda c: int(c) if c.isdigit() else ord(c) - 55
    return pd.Timestamp(year=yr, month=dig(s[3]), day=dig(s[4])).to_julian_date()


def read_mpcorb():
    """Numbered asteroids with H <= HMAX: elements, and the number from the readable designation."""
    cols = [(0, 7), (8, 13), (14, 19), (20, 25), (26, 35), (37, 46), (48, 57), (59, 68), (70, 79), (92, 103), (166, 194)]
    names = ["desig", "H", "G", "epoch", "M", "peri", "node", "incl", "e", "a", "name"]
    path = EPH / "MPCORB.DAT.gz"
    with gzip.open(path, "rt", errors="replace") as f:
        skip = next(i for i, line in enumerate(f) if line.startswith("-----")) + 1
    D = pd.read_fwf(path, colspecs=cols, names=names, skiprows=skip, compression="infer", dtype=str, encoding_errors="replace")
    for c in ["H", "G", "M", "peri", "node", "incl", "e", "a"]:
        D[c] = pd.to_numeric(D[c], errors="coerce")
    D["number"] = pd.to_numeric(D.name.str.extract(r"^\((\d+)\)", expand=False), errors="coerce")
    D = D[(D.H <= HMAX) & (D.e < 1) & D.a.notna() & D.number.notna()].copy()
    D["G"] = D.G.fillna(0.15); D["number"] = D.number.astype(int)
    D["jd"] = D.epoch.map({e: unpack_epoch(e) for e in D.epoch.unique()})
    return D.reset_index(drop=True)


def elements_to_state(D):
    """Heliocentric ecliptic osculating elements -> heliocentric equatorial (ICRF) position and velocity, AU and AU/day."""
    a, e = D.a.values, D.e.values
    i, O, w, M = (np.radians(D[c].values) for c in ["incl", "node", "peri", "M"])
    E = M.copy()
    for _ in range(50):
        E -= (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    n = K / a ** 1.5
    xp, yp = a * (np.cos(E) - e), a * np.sqrt(1 - e ** 2) * np.sin(E)
    Ed = n / (1 - e * np.cos(E))
    vxp, vyp = -a * np.sin(E) * Ed, a * np.sqrt(1 - e ** 2) * np.cos(E) * Ed
    cO, sO, cw, sw, ci, si = np.cos(O), np.sin(O), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
    P = np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si])
    Q = np.array([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si])
    r = P * xp + Q * yp; v = P * vxp + Q * vyp
    rot = np.array([[1, 0, 0], [0, np.cos(EPS), -np.sin(EPS)], [0, np.sin(EPS), np.cos(EPS)]])
    return (rot @ r).T, (rot @ v).T


def site(jd_utc):
    """Palomar's offset from the Earth's centre (AU, equator of date; the 1950s precession moves the 6,400 km vector by
    under 1%, a sub-arcsecond effect), for an array of UTC Julian dates."""
    jd_utc = np.atleast_1d(jd_utc)
    T = (jd_utc - 2451545.0) / 36525
    lst = np.radians((280.46061837 + 360.98564736629 * (jd_utc - 2451545.0) + 0.000387933 * T ** 2) % 360) + LON
    f = 1 / 298.257223563; C = 1 / np.sqrt(np.cos(LAT) ** 2 + (1 - f) ** 2 * np.sin(LAT) ** 2); Sf = (1 - f) ** 2 * C
    rc = (C + HGT / 6378.137) * np.cos(LAT) * RE_AU; rs = (Sf + HGT / 6378.137) * np.sin(LAT) * RE_AU
    return np.c_[rc * np.cos(lst), rc * np.sin(lst), np.full(len(jd_utc), rs)]


def delta_t(ts):
    """TT - UT in the 1950s (s), linear in the year."""
    return 29.15 + 0.38 * (ts.dt.year + ts.dt.day_of_year / 365.25 - 1950)


def sky(p, pv, obs):
    """Astrometric RA, Dec (deg) and distance (AU) of barycentric positions p (velocities pv) from observer obs, with one
    light-time iteration (under 0.1 day)."""
    d = p - obs; tau = np.linalg.norm(d, axis=-1, keepdims=True) / C_AU_D
    d = d - pv * tau; dist = np.linalg.norm(d, axis=-1)
    return np.degrees(np.arctan2(d[..., 1], d[..., 0])) % 360, np.degrees(np.arcsin(d[..., 2] / dist)), dist


def box(ra, dec, ra0, dec0):
    """Gnomonic offsets (deg) about each plate centre."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    c = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    return (np.degrees(np.cos(dec) * np.sin(ra - ra0) / c),
            np.degrees((np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / c), c)


def ephem():
    import assist
    return assist.Ephem(str(EPH / "de440.bsp"), str(EPH / "sb441-n16.bsp"))


def listing_batch(args):
    """One batch of orbits through every plate time: the asteroids in each plate's box, with V and their state there."""
    D, times = args
    import assist
    import rebound
    E = ephem(); t0 = time.time()
    r, v = elements_to_state(D)
    states = np.zeros((len(D), 6, len(times)))
    order = np.argsort(-times.tdb.values)                                    # backwards in time from the epochs
    for jd, idx in D.groupby("jd").indices.items():
        sim = rebound.Simulation(); sun = E.get_particle("sun", jd - E.jd_ref)
        for j in idx:
            sim.add(x=r[j, 0] + sun.x, y=r[j, 1] + sun.y, z=r[j, 2] + sun.z, vx=v[j, 0] + sun.vx, vy=v[j, 1] + sun.vy, vz=v[j, 2] + sun.vz)
        sim.t = jd - E.jd_ref; sim.ri_ias15.min_dt = 1e-3; sim.ri_ias15.epsilon = EPS
        ex = assist.Extras(sim, E)                                           # noqa: F841  (keeps the forces attached)
        buf, bv = np.zeros((len(idx), 3)), np.zeros((len(idx), 3))
        for q in order:
            sim.integrate(times.tdb.values[q] - E.jd_ref)
            sim.serialize_particle_data(xyz=buf, vxvyvz=bv)
            states[idx, :3, q] = buf; states[idx, 3:, q] = bv
    rows = []
    for q, tm in enumerate(times.itertuples()):
        tt = tm.tdb - E.jd_ref; e = E.get_particle("earth", tt); s = E.get_particle("sun", tt)
        obs = np.array([e.x, e.y, e.z]) + site(tm.jd_utc)[0]
        ra, dec, dist = sky(states[:, :3, q], states[:, 3:, q], obs)
        xi, eta, c = box(ra, dec, tm.ra_c, tm.dec_c)
        m = (c > 0) & (np.abs(xi) < BOX) & (np.abs(eta) < BOX)
        if not m.any():
            continue
        hs = states[m, :3, q] - np.array([s.x, s.y, s.z]); rh = np.linalg.norm(hs, axis=1)
        d = states[m, :3, q] - obs
        ph = np.arccos(np.clip(np.sum(hs * d, 1) / (rh * dist[m]), -1, 1)); t2 = np.tan(ph / 2); G = D.G.values[m]
        phi = (1 - G) * np.exp(-3.33 * t2 ** 0.63) + G * np.exp(-1.87 * t2 ** 1.22)
        V = D.H.values[m] + 5 * np.log10(rh * dist[m]) - 2.5 * np.log10(phi)
        k = V <= VMAX
        rows.append(pd.DataFrame({"plate": tm.plate, "when": tm.when, "number": D.number.values[m][k], "name": D.name.values[m][k],
                                  "H": D.H.values[m][k], "V": V[k], "ra": ra[m][k], "dec": dec[m][k], "tdb": tm.tdb, "jd_utc": tm.jd_utc,
                                  **{c_: states[m, j, q][k] for j, c_ in enumerate(["x", "y", "z", "vx", "vy", "vz"])}}))
    return (pd.concat(rows) if rows else pd.DataFrame()), len(D), time.time() - t0


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    N = pd.read_csv(ROOT / "data/plate_times_nominal.csv", index_col=0, parse_dates=["start", "log_start"])
    mid = N.start + pd.to_timedelta(N.exp / 2, unit="min")
    times = pd.DataFrame({"plate": N.index, "when": "header", "jd_utc": mid.map(pd.Timestamp.to_julian_date).values,
                          "tdb": (mid.map(pd.Timestamp.to_julian_date) + delta_t(mid) / 86400).values,
                          "ra_c": N.ra_c.values, "dec_c": N.dec_c.values})
    D = read_mpcorb(); t0 = time.time()
    batches = [(D.iloc[i:i + 1000], times) for i in range(0, len(D), 1000)]
    print(f"{len(D):,} numbered asteroids with H <= {HMAX}; {len(times)} plate times; {len(batches)} batches", flush=True)
    out = []
    with ProcessPoolExecutor(workers) as pool:
        for j, (L, n, dt) in enumerate(pool.map(listing_batch, batches)):
            out.append(L)
            if j % 10 == 0:
                print(f"batch {j}: {n} orbits, {len(L)} in boxes, {dt:.0f} s; {time.time() - t0:.0f} s so far", flush=True)
    L = pd.concat(out, ignore_index=True)
    L = L[L.V >= VMIN].sort_values(["plate", "when", "V"]).groupby(["plate", "when"]).head(NMAX)
    L[["plate", "number", "name", "H", "V", "ra", "dec", "jd_utc"]].to_parquet(ROOT / "data/asteroids_on_plates.parquet")
    print(f"listing: {len(L):,} asteroids on {L.plate.nunique()} plates, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
