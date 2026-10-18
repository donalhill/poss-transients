"""Exploratory: the numbers on the page come from analysis/build_site.py, which sets each night's expected count from a
smooth per-plate curve. This script uses an earlier per-night expectation (mean of nights within 8 degrees).

VASCO's nuclear test, reproduced and then rerun after allowing for where the telescope pointed.

1. VASCO's random-night test on its +/-1 day window: archive dates (its 26 nights), corrected dates, and corrected dates with
   each night's count divided by the number expected from the mean distance of its plates from the Milky Way.
2. Look-elsewhere on the 2026 paper's own scan of single-day offsets -3..+3: the strongest of the seven, against the
   strongest of seven in test calendars shifted by any number of days or by whole years (same season)."""
import contextlib, io
from pathlib import Path
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
ROOT = Path(__file__).resolve().parents[1]
g={"__file__":str(ROOT / "analysis/nuclear_lee.py"),"__name__":"x"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(ROOT / "analysis/nuclear_lee.py").read().split("windows = {")[0], g)
nt=g["nt"]; day0=pd.Timestamp("1949-01-01"); nights=g["nights"]
us=np.unique(((pd.to_datetime(nt[nt.Country.str.upper()=="USA"].Date)-day0).dt.days).values)
Nc, vc = g["N"], g["values"]
cat=pd.read_parquet(ROOT / "data/transients_with_plates.parquet")
strict=pd.to_datetime(pd.read_csv(ROOT / "VASCO-ML/data/poss_i_observing_dates_strict.csv").obs_date.astype(str).str[:10])
d=pd.to_datetime(cat.obs_date.astype(str).str[:10]); d=d[d.isin(strict)]
cnt=d.value_counts().reindex(strict, fill_value=0); Na=((cnt.index-day0).days).values; va=cnt.values.astype(float)
pl=pd.read_parquet(ROOT / "data/plates_corrected.parquet"); pl["b"]=np.abs(SkyCoord(pl.ra_c.values*u.deg, pl.dec_c.values*u.deg).galactic.b.deg)
nb=pl.groupby(pl.ut.dt.normalize()).b.mean().reindex(nights); ok=nb.notna().values
bb=nb.values   # expected count: mean of all nights whose plates lie within 8 degrees of the same distance (as build_site.py)
exp_night=np.array([vc[ok & (np.abs(bb-b0)<=8)].mean() if o else np.nan for b0,o in zip(bb,ok)])
adj=vc/exp_night
rng=np.random.default_rng(42)
def perm(N,v,n=20000):
    f=np.isfinite(v); N,v=N[f],v[f]
    lab=np.isin(N,(us[:,None]+np.array([-1,0,1])).ravel()); k=lab.sum()
    obs=v[lab].mean()/v[~lab].mean(); tot=v.sum(); R=np.empty(n)
    for i in range(n):
        idx=rng.choice(len(v),k,replace=False); m=v[idx].mean(); R[i]=m/((tot-v[idx].sum())/(len(v)-k))
    return obs,k,np.mean(R>=obs),R
for nm,(N,v) in {"archive raw":(Na,va),"corrected raw":(Nc,vc),"corrected sky-adjusted":(Nc,adj)}.items():
    obs,k,p,R=perm(N,v); print(f"{nm:24s} nights {k:2d} ratio {obs:.2f}  random-night p {p:.3f} (1 in {1/p:.0f})")

# look-elsewhere on VASCO's own scan: strongest of seven single-day offsets (-3..+3), real vs shifted calendars
lags = np.arange(-3, 4)
anys = np.array([s for s in range(-1100, 1101) if abs(s) >= 30])
yrs = np.array([365*k+dd for k in (-5,-4,-3,-2,-1,1,2,3,4,5) for dd in range(-20,21)])
def lagratios(N, v, shift):
    f = np.isfinite(v); N, v = N[f], v[f]; out = []
    for k in lags:
        lab = np.isin(N, us + shift + k)
        out.append(v[lab].mean()/v[~lab].mean() if lab.sum() >= 2 else np.nan)
    return np.array(out)
for nm,(N,v) in {"archive raw":(Na,va),"corrected raw":(Nc,vc),"corrected sky-adjusted":(Nc,adj)}.items():
    real = lagratios(N, v, 0); rmax = np.nanmax(real); kbest = lags[np.nanargmax(real)]
    for sn, sh in (("any shift", anys), ("same season", yrs)):
        E = np.array([lagratios(N, v, s) for s in sh])
        single = np.nanmean(E[:, np.nanargmax(real)] >= rmax)
        glob = np.nanmean(np.nanmax(E, axis=1) >= rmax)
        print(f"{nm:24s} {sn:11s}: strongest lag {kbest:+d} at {rmax:.2f}x; that lag alone p={single:.3f}; strongest-of-7 p={glob:.3f}")
