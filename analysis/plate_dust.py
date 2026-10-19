"""Interstellar dust towards each plate: E(B-V) at the plate centre from the Schlegel, Finkbeiner & Davis (1998) map,
queried from IRSA's dust service, with the plate's number of transients (every catalogue plate).
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
    d["n"] = n
    d.to_csv(ROOT / "data/plate_dust.csv")
    print(d.describe())
