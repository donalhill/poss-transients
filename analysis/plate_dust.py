"""Interstellar dust towards each plate: E(B-V) at the plate centre from the Schlegel, Finkbeiner & Davis (1998) map,
queried from IRSA's dust service; and E(B-V) averaged over the plate's searched area (6 x 6 deg, a 25 x 25 grid), from the
same map read locally (the SFD files, fetched once into data/external/dustmaps/sfd). Dust is patchy on scales smaller than
a plate, so the plate average is what filters a plate's transients. With the plate's number of transients.
python analysis/plate_dust.py -> data/plate_dust.csv"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa.irsa_dust import IrsaDust

sys.path.insert(0, str(Path(__file__).parent))
import sample                                                                  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SFD = ROOT / "data/external/dustmaps/sfd"
SFD_URL = "https://github.com/kbarbary/sfddata/raw/7a5fe7fadf086561ba4748756e59a4c51d0ec632/"                     # a mirror of the SFD files used by sfdmap


def plate_mean(pl):
    """E(B-V) averaged over each plate's central 6 x 6 deg, on a 25 x 25 grid of points."""
    import numpy as np
    import requests
    from dustmaps.sfd import SFDQuery
    from geometry import deproject
    SFD.mkdir(parents=True, exist_ok=True)
    for f in ("SFD_dust_4096_ngp.fits", "SFD_dust_4096_sgp.fits"):
        if not (SFD / f).exists():
            (SFD / f).write_bytes(requests.get(SFD_URL + f, timeout=600).content)
    q = SFDQuery(map_dir=str(SFD))
    g = np.linspace(-3, 3, 25); gx, gy = [a.ravel() for a in np.meshgrid(g, g)]
    out = []
    for _, p in pl.iterrows():
        ra, dec = deproject(gx, gy, np.full(gx.size, p.ra_c), np.full(gx.size, p.dec_c))
        out.append(float(np.mean(q(SkyCoord(ra * u.deg, dec * u.deg)))))
    return out

if __name__ == "__main__":
    cat, pl = sample.ours()
    n = cat.groupby("plate").size()
    out = ROOT / "data/plate_dust.csv"
    known = pd.read_csv(out).set_index("plate").ebv if out.exists() else pd.Series(dtype=float)   # the service is slow; keep what we have
    def ebv(q):
        p = pl.loc[q]
        for attempt in range(3):
            try:
                return float(IrsaDust.get_query_table(SkyCoord(p.ra_c * u.deg, p.dec_c * u.deg), section="ebv")["ext SFD ref"][0])
            except Exception:
                pass
        return float("nan")
    with ThreadPoolExecutor(8) as ex:
        todo = [q for q in n.index if q not in known.index or pd.isna(known.get(q))]
        d = pd.concat([known[known.index.isin(n.index)], pd.Series(dict(zip(todo, ex.map(ebv, todo))), dtype=float)]).rename("ebv").rename_axis("plate").to_frame()
    d["ebv_plate"] = plate_mean(pl.reindex(d.index))
    d["n"] = n
    d.to_csv(ROOT / "data/plate_dust.csv")
    print(d.describe())
