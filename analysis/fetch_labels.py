"""Fetch the label corner (last 600 rows, first 2000 columns) of POSS-I XE plate scans and save PNG crops for reading.

Usage: python analysis/fetch_labels.py XE641 XE089 ...  (or a file of plate names via @file)"""
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
HDR = ROOT / "data/dss_headers"
OUT = ROOT / "data/labels"
OUT.mkdir(parents=True, exist_ok=True)
DSS = "https://archive.stsci.edu/missions/dss/FitsArchive/Palomar_Observatory/XEsurvey/"
ROWS, COLS = 600, 2000


def card(h, key):
    for i in range(0, len(h), 80):
        if h[i:i + 8].strip() == key:
            return h[i + 10:i + 80].split("/")[0].strip().strip("'").strip()


def fetch(plate):
    out = OUT / f"{plate}.npy"
    if out.exists():
        return plate, "cached"
    files = sorted(HDR.glob(f"full_{plate.lower()}_*_00_0*.hdr"))
    if not files:                                       # the scan's FITS header, from the archive (first 20 blocks)
        listing = requests.get(DSS, timeout=60).text
        fn = sorted(set(re.findall(rf"{plate.lower()}_[0-9a-z]{{4}}_00_0\d\.fits", listing)))[0]
        HDR.mkdir(parents=True, exist_ok=True)
        files = [HDR / f"full_{fn[:-5]}.hdr"]
        files[0].write_bytes(requests.get(DSS + fn, headers={"Range": "bytes=0-57599"}, timeout=60).content)
    f = files[0]
    h = f.read_bytes().decode("ascii", "replace")
    nx, ny = int(card(h, "NAXIS1")), int(card(h, "NAXIS2"))
    data0 = (h.find("END     ") // 2880 + 1) * 2880
    rows = list(range(ny - ROWS, ny))
    a = np.zeros((ROWS, COLS), dtype=np.float32)
    url = DSS + f.stem.removeprefix("full_") + ".fits"
    for k in range(0, ROWS, 150):  # multi-range requests: only the first COLS pixels of each row
        chunk = rows[k:k + 150]
        rng = ",".join(f"{data0 + 2 * nx * r}-{data0 + 2 * nx * r + 2 * COLS - 1}" for r in chunk)
        resp = requests.get(url, headers={"Range": f"bytes={rng}"}, timeout=300)
        boundary = resp.headers["Content-Type"].split("boundary=")[1].encode()
        parts = [p for p in resp.content.split(b"--" + boundary) if b"content-range" in p.lower()]
        for p in parts:
            hdr, body = p.split(b"\r\n\r\n", 1)
            first = int(hdr.split(b"bytes ")[1].split(b"-")[0])
            r = (first - data0) // (2 * nx)
            a[r - (ny - ROWS)] = np.frombuffer(body[:2 * COLS], dtype=">i2")
    np.save(out, a)
    return plate, "ok"


if __name__ == "__main__":
    args = sys.argv[1:]
    plates = Path(args[0][1:]).read_text().split() if args and args[0].startswith("@") else args
    with ThreadPoolExecutor(8) as ex:
        for plate, msg in ex.map(fetch, plates):
            print(plate, msg, flush=True)
    # contact sheets, 6 labels per page
    ok = [p for p in plates if (OUT / f"{p}.npy").exists()]
    for k in range(0, len(ok), 6):
        fig, axes = plt.subplots(6, 1, figsize=(14, 18))
        for ax, p in zip(axes, ok[k:k + 6]):
            a = np.load(OUT / f"{p}.npy")[100:600, 50:1100]
            lo, hi = np.percentile(a, [1, 99.7])
            ax.imshow(a, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)
            ax.set_title(p, fontsize=12, loc="left"); ax.set_xticks([]); ax.set_yticks([])
        for ax in axes[len(ok[k:k + 6]):]:
            ax.axis("off")
        fig.savefig(ROOT / f"figs/labels_sheet_{k // 6:02d}.png", dpi=70, bbox_inches="tight")
        plt.close(fig)
    print("sheets:", (len(ok) + 5) // 6)
