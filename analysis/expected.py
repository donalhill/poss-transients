"""Each plate's expected count of transients with no effect: the pointing model. A negative binomial in the plate's
count with offset log(area) and smooth terms in where the plate points: |b| (a spline), dust averaged over the plate's
central square (SFD, read locally: make data/external/dustmaps/sfd), ecliptic latitude, the plate's depth, and the row of
plates centred on the equator (eq_row, build_sample.py: on the same sky these plates hold far fewer of VASCO's
transients than the row to their north). The per-plate scatter is in the dispersion, fitted by profile likelihood.
Used by the nuclear test (analysis/nuclear.py) and the pooled shadow counts (analysis/shadow_stats.py)."""
from pathlib import Path

import astropy.units as u
import numpy as np
import patsy
import statsmodels.api as sm
from astropy.coordinates import SkyCoord
from scipy.optimize import minimize_scalar

import shadow as SH
from geometry import deproject

ROOT = Path(__file__).resolve().parents[1]


def TERMS(pl):
    """The model's terms; the equatorial row's only when the plates table has it (report/vasco_stats.py also fits the
    model without it, to measure the row's shortfall)."""
    return "bs(b, df=4) + lebv + ecl_lat + depth" + (" + eq_row" if "eq_row" in pl else "")


def design(pl):
    """Adds b, ecl_lat, lebv, and expo to the plates table and returns the design matrix and the offset."""
    sky = SkyCoord(pl.ra_c.values * u.deg, pl.dec_c.values * u.deg)
    pl["b"] = np.abs(sky.galactic.b.deg); pl["ecl_lat"] = np.abs(sky.barycentrictrueecliptic.lat.deg)
    from dustmaps.sfd import SFDQuery
    q = SFDQuery(map_dir=str(ROOT / "data/external/dustmaps/sfd"))
    dust = []
    for name, p in pl.iterrows():
        ra, dec = deproject(SH.GX, SH.GY, np.full(SH.GX.size, p.ra_c), np.full(SH.GX.size, p.dec_c))
        dust.append(float(np.mean(q(SkyCoord(ra * u.deg, dec * u.deg)))))
    pl["lebv"] = np.log10(dust)
    pl["expo"] = SH.GX.size * SH.CELL ** 2                                   # the central square's area, deg^2: the model
                                                                             # fits a density, the same area on every plate
    Z = lambda c: (pl[c] - pl[c].mean()) / pl[c].std()
    X = np.asarray(patsy.dmatrix(TERMS(pl), pl.assign(lebv=Z("lebv"), ecl_lat=Z("ecl_lat"),
                                                      depth=Z("depth")), NA_action="raise"))
    off = np.log(pl.expo.to_numpy())
    assert np.isfinite(off).all()
    return X, off


def fit(y, X, off, use):
    """The expected count of every plate from the model fitted on the plates in use, and the dispersion alpha."""
    # dispersion by profile likelihood
    at = lambda la: sm.GLM(y[use], X[use], offset=off[use], family=sm.families.NegativeBinomial(alpha=np.exp(la))).fit()
    la = minimize_scalar(lambda la: -at(la).llf, bounds=(np.log(1e-4), np.log(10)), method="bounded", options={"xatol": 1e-4}).x
    return at(la).predict(X, offset=off), float(np.exp(la))
