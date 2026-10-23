"""Header facts for every POSS-I red (XE) plate scan in STScI's DSS archive: serial label, archive date and time,
exposure, and plate centre. The FITS headers (first 20 blocks of each scan) are cached in data/dss_headers/.
Writes data/xe_headers.parquet (one row per scan file) and data/xe_plate_centres.csv (one row per plate, first scan).
python analysis/plate_headers.py"""
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
HDR = ROOT / "data/dss_headers"
DSS = "https://archive.stsci.edu/missions/dss/FitsArchive/Palomar_Observatory/XEsurvey/"


def header(fn):
    f = HDR / f"full_{fn[:-5]}.hdr"
    if not f.exists():
        f.write_bytes(requests.get(DSS + fn, headers={"Range": "bytes=0-57599"}, timeout=120).content)
    return fits.Header.fromstring(f.read_bytes().decode("ascii", "ignore"))


def row(fn):
    h = header(fn)
    m = re.match(r"(\d{4}-\d\d-\d\d)T(\d\d):(\d\d)", h["DATE-OBS"])
    lab = h.get("PLTLABEL", "").strip()
    num = re.sub(r"\D", "", lab)
    ra = 15 * (h["PLTRAH"] + h["PLTRAM"] / 60 + h["PLTRAS"] / 3600)
    dec = h["PLTDECD"] + h["PLTDECM"] / 60 + h["PLTDECS"] / 3600
    return dict(file=fn[:-5], plate=h["REGION"].strip()[:5], plateid=h["PLATEID"].strip(), label=float(num) if num else None,
                labelraw=lab, date=pd.Timestamp(m.group(1)), hh=int(m.group(2)), mm=int(m.group(3)),
                exp=float(h["EXPOSURE"]), ra=ra, dec=-dec if h.get("PLTDECSN", "+").strip() == "-" else dec)


if __name__ == "__main__":
    HDR.mkdir(parents=True, exist_ok=True)
    files = sorted(set(re.findall(r"xe\d{3}_[0-9a-z]{4}_00_0\d\.fits", requests.get(DSS, timeout=60).text)))
    with ThreadPoolExecutor(8) as ex:
        t = pd.DataFrame(list(ex.map(row, files)))
    t.drop(columns=["ra", "dec"]).to_parquet(ROOT / "data/xe_headers.parquet")
    t[t.file.str.endswith("_00_00")].drop_duplicates("plate")[["plate", "ra", "dec"]].to_csv(ROOT / "data/xe_plate_centres.csv", index=False)
    print(f"{len(t)} scans of {t.plate.nunique()} plates")
