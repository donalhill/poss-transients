"""The example plate's scan as a small image for the page: the scan
averaged in 20 x 20 pixel blocks, stretched by asinh between its 1st and 99.8th percentiles inside the scan's edge (blocks below zero, where the scanned label
and wedges wrap round, shown as sky), drawn
with stars bright and the scan's y axis up (the chart's), with the scan's SHA-256 recorded beside the image.
usage: python report/plate_image.py SCAN.fits OUT.png"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from matplotlib import image

scan, out = Path(sys.argv[1]), Path(sys.argv[2])
B = 20
with fits.open(scan, memmap=True) as h:
    d = h[0].data
    ny, nx = (d.shape[0] // B) * B, (d.shape[1] // B) * B
    a = d[:ny, :nx].astype(np.float32).reshape(ny // B, B, nx // B, B).mean(axis=(1, 3))   # rows are the scan's y, columns its x
valid = a > 0                                                               # the scanned label and wedges wrap round below zero: not sky
inner = a[5:-5, 5:-5][valid[5:-5, 5:-5]]
lo, hi = np.percentile(inner, [1, 99.8])
a = np.where(valid, a, np.median(inner))
t = np.clip((a - lo) / (hi - lo), 0, 1)
t = np.arcsinh(t * 10) / np.arcsinh(10)                                     # soft: faint stars visible, bright ones not saturated
image.imsave(out, np.flipud(t), cmap="gray", vmin=0, vmax=1)                # row 0 at the bottom, as the chart's y axis
sha = hashlib.sha256(scan.read_bytes()).hexdigest()
out.with_suffix(".json").write_text(json.dumps({"scan": scan.name, "sha256": sha, "block_px": B, "shape": list(a.shape)}))
print(f"{out}: {a.shape[1]} x {a.shape[0]} from {scan.name} ({sha[:12]})")
