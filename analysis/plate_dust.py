"""Interstellar dust towards each plate: E(B-V) at the plate centre from the Schlegel, Finkbeiner & Davis (1998) map,
queried from IRSA's dust service, with the plate's number of transients (VASCO's nights).
python analysis/plate_dust.py -> data/plate_dust.csv"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa.irsa_dust import IrsaDust

ROOT = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    pl = pd.read_parquet(ROOT / "data/plates_corrected.parquet")
    cat = pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
    strict = set(pd.read_csv(ROOT / "data/external/poss_i_observing_dates_strict.csv").obs_date.astype(str).str[:10])
    n = cat[cat.obs_date.astype(str).str[:10].isin(strict) & (cat.plate != "")].groupby("plate").size()
    def ebv(q):
        p = pl.loc[q]
        for attempt in range(3):
            try:
                return float(IrsaDust.get_query_table(SkyCoord(p.ra_c * u.deg, p.dec_c * u.deg), section="ebv")["ext SFD ref"][0])
            except Exception:
                pass
        return float("nan")
    with ThreadPoolExecutor(8) as ex:
        d = pd.DataFrame({"plate": list(n.index), "ebv": list(ex.map(ebv, n.index))}).set_index("plate")
    d["n"] = n
    d.to_csv(ROOT / "data/plate_dust.csv")
    print(d.describe())
