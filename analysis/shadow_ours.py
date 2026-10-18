"""Earth's-shadow test with our own transient extraction (frozen_selection.py, thresholds fixed on XE320), on the 22
southern plates that VASCO's catalogue leaves out and on VASCO's nine shadow plates, re-measured the same way.

Per plate: candidates inside the shadow (geocentric, 8.5 deg about the anti-sun point at geostationary height, at the
plate's mid-exposure time: Palomar's log for the southern plates, the corrected dates for VASCO's) against the number
expected from the plate's shadowed area, weighted by how candidates are spread across a plate (pooled over all these
plates, as a function of distance from the plate centre). Statistic: R = observed / expected, pooled (1 = no deficit).
Null: shadow footprints pasted onto plates of the set, in a random one of the 8 orientations of the square, using only
each host's own unshadowed sky (whole plates if wholly in shadow).
All plates are used, as VASCO's catalogue keeps plates on the Milky Way; the result is also given without the plates where
more than half of random positions have a Gaia or Pan-STARRS source within 5" (crowding.py). S 92 is used only north of
dec -29.5, where Pan-STARRS covers it.
python analysis/shadow_ours.py   -> data/shadow_ours.csv, data/shadow_ours_grid.npz (for shadow_phi_combined.py)"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
from scan_io import centre, pixel_arcsec, plate_header, plate_to_sky   # noqa: E402
from shadow import SHADOW_DEG, antisun_geocentric                        # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(7)
log = pd.read_csv(ROOT / "data/shadow_plates_log.csv").astype({"plate_number": str})
south = pd.read_csv(ROOT / "data/southern_shadow_scans.csv").astype({"plate_number": str}).merge(
    log[["plate_number", "evening_date", "ra", "dec", "ut"]], on=["plate_number", "evening_date", "ra", "dec"])   # numbers repeat
north = pd.read_csv(ROOT / "data/north_shadow_scans.csv").astype({"plate_number": str})
north["ut"] = pd.read_parquet(ROOT / "data/plates_corrected.parquet").ut.reindex(north.plate).astype(str).values
plates = pd.concat([south.assign(set="south"), north.assign(set="north")], ignore_index=True)
HALF, STEP = 2.9, 0.05
g = np.arange(-HALF + STEP / 2, HALF, STEP)                     # symmetric, so the 8 orientations map cells to cells
GX, GY = np.meshgrid(g, g)
NG = len(g)
RB = np.linspace(0, HALF, 11)                                   # bins of distance from the centre (Chebyshev)
cell_r = np.maximum(np.abs(GX), np.abs(GY))


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


crowd = pd.read_csv(ROOT / "data/crowding.csv").drop_duplicates("plate", keep="last").set_index("plate").veto   # crowding.py


rows, masks, usable, cells = [], {}, {}, {}
for p in plates.itertuples():
    f = ROOT / f"data/scans/selected_{p.plate}.parquet"
    if not f.exists():
        print(f"{p.plate} ({p.plate_number}): no selection yet, skipped"); continue
    h = plate_header(p.plate)
    xc, yc = centre(p.plate); scale = pixel_arcsec(p.plate) / 3600
    ra, dec = plate_to_sky(h, xc + GX.ravel() / scale, yc + GY.ravel() / scale)
    veto = crowd.get(p.plate, np.nan)
    shadow = antisun_geocentric(ra, dec, Time(np.full(ra.size, pd.Timestamp(p.ut).to_datetime64()))) < SHADOW_DEG
    use = np.ones(ra.size, bool) if p.plate_number != "S 92" else dec > -29.5
    s = pd.read_parquet(f)
    ix = np.clip(np.round(((s.x - xc) * scale - g[0]) / STEP).astype(int), 0, NG - 1)
    iy = np.clip(np.round(((s.y - yc) * scale - g[0]) / STEP).astype(int), 0, NG - 1)
    keep = (s.r < HALF).values & use.reshape(NG, NG)[iy, ix]
    rows.append(dict(plate=p.plate, plate_number=p.plate_number, set=p.set, dec=p.dec, veto=round(veto, 3),
                     shadow_area=round(float(shadow[use].mean()), 3), n=int(keep.sum())))
    masks[p.plate] = shadow.reshape(NG, NG); usable[p.plate] = use.reshape(NG, NG)
    cells[p.plate] = (iy[keep], ix[keep])
    print(f"{p.plate} ({p.plate_number}): {keep.sum()} candidates, shadow covers {shadow[use].mean():.0%}, "
          f"catalogue veto at random positions {veto:.0%}", flush=True)

t = pd.DataFrame(rows).set_index("plate")
D8 = [lambda m: m, np.rot90, lambda m: np.rot90(m, 2), lambda m: np.rot90(m, 3), np.fliplr, np.flipud, np.transpose,
      lambda m: np.rot90(m, 2).T]


def test(K, label):
    """R = observed / expected in the shadow, pooled over plates K, against pasted footprints."""
    if not K:
        print(f"\n{label}: no plates processed yet"); return None, None
    cnt = sum(np.histogram(cell_r[cells[q]], RB)[0] for q in K)      # spread of candidates with distance from centre
    area = sum(np.histogram(cell_r[usable[q]], RB)[0] for q in K)
    W = (cnt / area)[np.clip(np.digitize(cell_r, RB) - 1, 0, len(RB) - 2)]

    def obs_exp(q, mask, use=None):
        use = usable[q] if use is None else use
        k = use[cells[q]]
        c = (cells[q][0][k], cells[q][1][k])
        return mask[c].sum(), k.sum() * (W * mask * use).sum() / (W * use).sum()

    # hosts lend only sky where nothing could have vanished (a plate wholly in shadow is kept whole)
    host = {q: usable[q] if masks[q][usable[q]].mean() >= 0.97 else usable[q] & ~masks[q] for q in K}

    oe = np.array([obs_exp(q, masks[q]) for q in K])
    R = oe[:, 0].sum() / oe[:, 1].sum()
    null = np.empty(10000)
    for n in range(null.size):
        o = e = 0.0
        for q in K:
            a, b = obs_exp(q, D8[rng.integers(8)](masks[K[rng.integers(len(K))]]), host[q])
            o += a; e += b
        null[n] = o / e
    if label == "All plates":                                   # for shadow_phi_combined.py
        np.savez(ROOT / "data/shadow_ours_grid.npz", plates=np.array(K), sets=t.set[K].values, n=t.n[K].values.astype(int), W=W,
                 masks=np.array([masks[q] for q in K]), usable=np.array([usable[q] for q in K]),
                 cells=np.array([np.ravel_multi_index(cells[q], (NG, NG)) for q in K], dtype=object))
    print(f"\n{label}: {len(K)} plates, observed {oe[:, 0].sum():.0f}, expected {oe[:, 1].sum():.1f}, R = {R:.3f}")
    print(f"  pasted-footprint null: median {np.median(null):.3f}, 16-84% {np.quantile(null, .16):.3f}-{np.quantile(null, .84):.3f}; "
          f"share at or below R {np.mean(null <= R):.3f}")
    return oe, null


oe, null = test(list(t.index), "All plates")
t["observed"], t["expected"] = oe[:, 0], oe[:, 1].round(1)
test(list(t.index[t.set == "south"]), "Southern plates")
test(list(t.index[t.set == "north"]), "VASCO's nine plates, our transients")
test(list(t.index[t.veto <= 0.5]), "Without the crowded plates (veto > 50%)")
print(t.to_string())
t.to_csv(ROOT / "data/shadow_ours.csv")
