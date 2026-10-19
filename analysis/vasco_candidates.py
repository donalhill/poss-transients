"""VASCO's five shortlisted alignments (Villarroel et al. 2025, PASP, the transient tables of its candidates; positions in
data/vasco_2025_candidates.csv): which of their 24 transients are in the released catalogue, and where each candidate
lies on its plate (offset from the plate centre, in the frame of data/xe_plate_centres.csv).
python analysis/vasco_candidates.py -> data/vasco_candidates.json"""
import json
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    t = pd.read_csv(ROOT / "data/vasco_2025_candidates.csv")
    c = SkyCoord(t.ra_hms.values, t.dec_dms.values, unit=(u.hourangle, u.deg))
    cat = pd.read_csv(ROOT / "data/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv", usecols=["ra", "dec", "obs_date"])
    idx, sep, _ = c.match_to_catalog_sky(SkyCoord(cat.ra.values * u.deg, cat.dec.values * u.deg))
    t["nearest_arcsec"] = sep.arcsec.round(1)
    print(f"{(sep.arcsec < 5).sum()} of {len(t)} transients within 5\" of a released transient:")
    print(t[sep.arcsec < 5].to_string(index=False))
    cen = pd.read_csv(ROOT / "data/xe_plate_centres.csv").set_index("plate")
    hdr = pd.read_parquet(ROOT / "data/xe_headers.parquet").drop_duplicates("plate").set_index("plate")
    out = {"n_transients": len(t), "n_in_catalogue": int((sep.arcsec < 5).sum()), "candidates": []}
    for k, g in t.groupby("candidate"):
        cc = SkyCoord(g.ra_hms.values, g.dec_dms.values, unit=(u.hourangle, u.deg))
        ra, dec = np.mean(cc.ra.deg), np.mean(cc.dec.deg)
        r0, d0 = np.radians(cen.ra.values), np.radians(cen.dec.values); r1, d1 = np.radians(ra), np.radians(dec)
        cosc = np.sin(d0) * np.sin(d1) + np.cos(d0) * np.cos(d1) * np.cos(r1 - r0)
        xi = np.degrees(np.cos(d1) * np.sin(r1 - r0) / cosc); eta = np.degrees((np.cos(d0) * np.sin(d1) - np.sin(d0) * np.cos(d1) * np.cos(r1 - r0)) / cosc)
        on = (np.abs(xi) < 3.2) & (np.abs(eta) < 3.2) & (cosc > 0)
        same = [q for q in cen.index[on] if q in hdr.index and abs((hdr.loc[q].date - pd.Timestamp(g.date.iloc[0])).days) <= 1]
        for q in same:
            i = list(cen.index).index(q)
            out["candidates"].append({"candidate": int(k), "plate": q, "xi": round(float(xi[i]), 3), "eta": round(float(eta[i]), 3),
                                      "n": len(g), "nearest_arcsec": float(g.nearest_arcsec.min())})
            print(f"candidate {k} ({g.date.iloc[0]}, {len(g)} transients, nearest released transient {g.nearest_arcsec.min():.0f}\"): "
                  f"plate {q}, {abs(xi[i]):.2f} deg {'east' if xi[i] > 0 else 'west'} and {abs(eta[i]):.2f} deg {'north' if eta[i] > 0 else 'south'} of the centre")
    (ROOT / "data/vasco_candidates.json").write_text(json.dumps(out))
