"""Download a full DSS plate scan from STScI in parallel byte ranges, resuming any partial file.
python analysis/fetch_scan.py xe320_06co_00_00"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
name = sys.argv[1]
url = f"https://archive.stsci.edu/missions/dss/FitsArchive/XEsurvey/{name}.fits"
out = ROOT / f"data/scans/{name}.fits"
size = int(requests.head(url, timeout=60).headers["Content-Length"])
have = out.stat().st_size if out.exists() else 0
existing = sorted(int(q.suffix[5:]) for q in out.parent.glob(out.name.replace(".fits", ".part*")))
if existing:
    have = existing[0]                                    # resume the earlier plan
NP = 8
step = (size - have + NP - 1) // NP
parts = [(have + k * step, min(size, have + (k + 1) * step) - 1) for k in range(NP) if have + k * step < size]


def get(rng):
    a, b = rng
    p = out.with_suffix(f".part{a}")
    for attempt in range(20):                              # the archive drops connections: resume where each left off
        got = p.stat().st_size if p.exists() else 0
        if got >= b - a + 1:
            return p
        try:
            with requests.get(url, headers={"Range": f"bytes={a + got}-{b}"}, stream=True, timeout=600) as r, open(p, "ab") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
        except requests.RequestException:
            continue
    raise SystemExit(f"part {a} incomplete")


with ThreadPoolExecutor(NP) as ex:
    files = list(ex.map(get, parts))
with open(out, "ab") as fh:
    for p in files:
        fh.write(p.read_bytes()); p.unlink()
print(out, out.stat().st_size, "of", size)
