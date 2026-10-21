"""Fitting a star (or any compact source) on a scan: an elliptical Gaussian on a small stamp. Used by
plate_depth.py."""
import numpy as np
from scipy.optimize import curve_fit

STAMP = 6                                                              # pixels each side of the fitted source


def gauss(xy, a, x0, y0, sx, sy, th, c):
    x, y = xy; ct, st = np.cos(th), np.sin(th)
    xr, yr = (x - x0) * ct + (y - y0) * st, -(x - x0) * st + (y - y0) * ct
    return a * np.exp(-0.5 * ((xr / sx) ** 2 + (yr / sy) ** 2)) + c


def fit(im, x, y):
    """Elliptical Gaussian on a stamp around (x, y): peak, FWHM (geometric mean, pixels), elongation 1 - b/a, angle."""
    i, j = int(round(y)), int(round(x))
    if i < STAMP or j < STAMP or i >= im.shape[0] - STAMP or j >= im.shape[1] - STAMP:
        return None
    s = im[i - STAMP:i + STAMP + 1, j - STAMP:j + STAMP + 1]; yy, xx = np.mgrid[-STAMP:STAMP + 1, -STAMP:STAMP + 1]
    try:
        p, _ = curve_fit(gauss, (xx.ravel(), yy.ravel()), s.ravel(), p0=[s.max() - np.median(s), x - j, y - i, 1.2, 1.2, 0, np.median(s)],
                         bounds=([0, -2, -2, 0.4, 0.4, -np.pi, -np.inf], [np.inf, 2, 2, 6, 6, np.pi, np.inf]), maxfev=4000)
    except (RuntimeError, ValueError):
        return None
    a_, b_ = max(p[3], p[4]), min(p[3], p[4]); th = (p[5] + (np.pi / 2 if p[4] > p[3] else 0)) % np.pi
    return {"peak": float(p[0]), "fwhm": float(2.3548 * np.sqrt(p[3] * p[4])), "elong": float(1 - b_ / a_), "angle": float(np.degrees(th)),
            "xc": float(j + p[1]), "yc": float(i + p[2])}
