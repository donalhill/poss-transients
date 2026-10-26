"""The numbers the page quotes on VASCO's catalogue itself (the transients each plate holds, before any test):
how far the counts spread beyond counting noise; the share of plates denser at their edges than their middles; how
counts fall towards the Milky Way and behind dust; the equatorial row's shortfall against the pointing model
(analysis/expected.py), and on the same sky where it overlaps the next row north (VASCO's control positions measuring
the area searched); how much of the spread in log counts the pointing model explains; and the survey's years.
Counts are of the whole catalogue, every transient once; the plates' pointing-model inputs come from a sample (build_sample.py).
Under "fig", the data of the section's charts from the same counts and the same fit (report/figures.py): every
plate's count, expected count, and distance from the Milky Way; the plates in eight groups by expected count; a smooth
curve of log count against distance from the Milky Way with the middle half of plates about it; and the density of
transients and of VASCO's control positions by distance from the plate's centre, and across the plate.
usage: python report/vasco_stats.py SAMPLE_DIR OUT.json   (SAMPLE_DIR: build_sample.py's)"""
import json
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "analysis"), str(ROOT / "vasco")]
import expected as EX                                                        # noqa: E402
import vasco_sample                                                          # noqa: E402
from geometry import deproject, tangent                                      # noqa: E402

d, out = Path(sys.argv[1]), Path(sys.argv[2])
cat = vasco_sample.ours()[0]                                                 # every transient, the whole plate
pl = pd.read_parquet(d / "plates.parquet")                                   # the plates, with the pointing model's inputs
N = cat.groupby("plate").size().reindex(pl.index).fillna(0).astype(int)
avg = float(N.mean()); noise = float(np.sqrt(avg)); lo, hi = avg - 2 * noise, avg + 2 * noise
res = {"transients": int(len(cat)), "plates": int(len(pl)), "min": int(N.min()), "max": int(N.max()), "avg": round(avg),
       "noise": round(noise), "band_lo": round(lo), "band_hi": round(hi), "in_band": int(((N >= lo) & (N <= hi)).sum()),
       "over1000": int((N > 1000).sum()), "years": [int(pl.ut.min().year), int(pl.ut.max().year)]}
# the share of plates whose band 2.0-2.75 deg from the centre is denser than their inner 1.5 deg
r = np.maximum(cat.xi.abs(), cat.eta.abs())
res["edge_share"] = round(float(np.mean([((g < 2.75) & (g >= 2.0)).sum() / (5.5 ** 2 - 4.0 ** 2) > (g < 1.5).sum() / 3.0 ** 2
                                         for _, g in r.groupby(cat.plate)])), 3)
# crowding: the typical plate near the Milky Way and far from it
b = np.abs(SkyCoord(pl.ra_c.values * u.deg, pl.dec_c.values * u.deg).galactic.b.deg)
has = N.values > 0
res["crowd_inner"] = int(np.median(N.values[has & (b < 5)])); res["crowd_outer"] = int(round(float(np.median(N.values[has & (b > 40)])), -1))
# dust, away from the Milky Way (|b| > 20 deg): the dustiest plates (E(B-V) >= 0.15) against the least dusty (< 0.05)
du = pd.read_csv(ROOT / "data/plate_dust.csv").set_index("plate").reindex(pl.index)
far = b > 20; dusty, clear = far & (du.ebv >= 0.15).values, far & (du.ebv < 0.05).values
res.update({"dust_n_dusty": int(dusty.sum()), "dust_median_dusty": int(np.median(N.values[dusty])),
            "dust_n_clear": int(clear.sum()), "dust_median_clear": int(np.median(N.values[clear]))})
# the pointing model fitted to every plate's count: the share of the spread in log counts it explains; and the plates
# centred on the equator (VASCO's southern edge) against the other rows, from the model without its equatorial term
fit_ = lambda P: EX.fit(N.values.astype(float), *EX.design(P), np.ones(len(P), bool))[0]
E = fit_(pl.copy()); E_norow = fit_(pl.drop(columns="eq_row"))
pos = N.values > 0
res["pointing_share"] = round(float(np.corrcoef(np.log(N.values[pos]), np.log(E[pos]))[0, 1] ** 2), 3)
row = (np.round(pl.dec_c / 6) * 6).astype(int).values; q = N.values / E_norow
res["eq_plates"] = int((row == 0).sum()); res["eq_share"] = round(float(np.median(q[row == 0]) / np.median(q[row > 0])), 3)
# the same sky on both rows, where the equatorial plates overlap the row to their north (dec +2.25 to +3): each row's
# transients per control position
dec_t = deproject(cat.xi.values, cat.eta.values, pl.ra_c.reindex(cat.plate).values, pl.dec_c.reindex(cat.plate).values)[1]
ctrl = pd.read_csv(ROOT / "data/external/XE_final_random_sources.csv", usecols=["Dec", "Name"])


def per_control(rows):
    S = set(pl.index[np.isin(row, rows)])
    nt = ((dec_t > 2.25) & (dec_t < 3.0) & cat.plate.isin(S).values).sum(); nc = ((ctrl.Dec > 2.25) & (ctrl.Dec < 3.0) & ctrl.Name.isin(S)).sum()
    return nt / nc


res["eq_same_sky"] = round(float(per_control([0]) / per_control([6])), 3)

# ---- the charts' data
import patsy                                                                 # noqa: E402
import statsmodels.api as sm                                                 # noqa: E402
r1 = lambda v: round(float(v), 1)
F = {"plates": [{"id": q, "N": int(N[q]), "E": r1(E[i]), "b": r1(b[i])} for i, q in enumerate(pl.index)],
     "mean": round(avg, 1), "lo": round(lo, 1), "hi": round(hi, 1)}
# calibration: the plates in eight groups by expected count, each group's mean expected and mean actual count, the latter's
# 68% range from resampling the group's plates
rng = np.random.default_rng(5); cal = []
for g in np.array_split(np.argsort(E), 8):
    n_ = N.values[g]; bm = n_[rng.integers(0, len(n_), (4000, len(n_)))].mean(1)
    cal.append({"E": r1(E[g].mean()), "N": r1(n_.mean()), "lo": r1(np.percentile(bm, 16)), "hi": r1(np.percentile(bm, 84)), "n": int(len(n_))})
F["calibration"] = cal
# log count against distance from the Milky Way: a straight-line fit on a smooth curve of the distance (the typical plate),
# and the middle half of plates about it from the residuals' robust spread
pc = pd.DataFrame({"absb": b[has], "n": N.values[has]}); top = float(np.ceil(pc.absb.max()))
Xc = patsy.dmatrix(f"bs(absb, df=4, lower_bound=0, upper_bound={top})", pc, return_type="dataframe")
fc = sm.OLS(np.log(pc.n.values), Xc).fit(); resid = np.log(pc.n.values) - fc.fittedvalues.values
sd = 1.4826 * float(np.median(np.abs(resid - np.median(resid))))
bg = np.arange(0.0, top + 0.01, 1.0); lg = np.asarray(fc.predict(patsy.build_design_matrices([Xc.design_info], {"absb": bg})[0]))
F["crowd_fit"] = {"b": bg.tolist(), "mid": np.exp(lg).round(1).tolist(), "q25": np.exp(lg - 0.6745 * sd).round(1).tolist(), "q75": np.exp(lg + 0.6745 * sd).round(1).tolist()}
# density by distance from the plate's centre (square rings, max(|xi|, |eta|), a quarter degree wide, out to 2.75 deg: the
# controls thin beyond, the area searched ending short of 3 deg on some plates) over the density
# 0.5 to 1 degree out, the transients' 68% range from resampling the plates; VASCO's control positions alongside
ctl = pd.read_csv(ROOT / "data/external/XE_final_random_sources.csv", usecols=["RA", "Dec", "Name"])
ctl = ctl[ctl.Name.isin(pl.index)]
cx, cy = tangent(ctl.RA.values, ctl.Dec.values, pl.ra_c.reindex(ctl.Name).values, pl.dec_c.reindex(ctl.Name).values)
ring = np.arange(0, 2.76, 0.25); area = (2 * ring[1:]) ** 2 - (2 * ring[:-1]) ** 2; ref = (ring[:-1] >= 0.5) & (ring[:-1] < 1.0)
cp = cat[cat.plate.isin(pl.index)]; rt = np.maximum(cp.xi.abs(), cp.eta.abs()).values
dt = np.histogram(rt, ring)[0] / area; dc = np.histogram(np.maximum(np.abs(cx), np.abs(cy)), ring)[0] / area
per = np.array([np.histogram(rt[cp.plate.values == q], ring)[0] for q in pl.index[N.values > 0]]) / area
boot = np.array([(lambda s: s / s[ref].mean())(per[rng.integers(0, len(per), len(per))].sum(0)) for _ in range(2000)])
F["radial"] = {"r": (ring[:-1] + 0.125).tolist(), "transients": (dt / dt[ref].mean()).round(3).tolist(), "control": (dc / dc[ref].mean()).round(3).tolist(),
               "lo": np.percentile(boot, 16, 0).round(3).tolist(), "hi": np.percentile(boot, 84, 0).round(3).tolist()}
# across the plate: 0.2-degree cells, each plate weighted equally, over the density within 1 degree of the centre (rows
# south to north, columns west to east)
G = np.arange(-3.0, 3.01, 0.2); mid = G[:-1] + 0.1; ctr = (np.abs(mid)[:, None] < 1) & (np.abs(mid)[None, :] < 1)
cmap = lambda x, y, w: (lambda h: (h / h[ctr].mean()).round(3))(np.histogram2d(y, x, [G, G], weights=w)[0])
F["map"] = {"edges": G.round(2).tolist(), "transients": cmap(cp.xi.values, cp.eta.values, 1 / cp.groupby("plate").plate.transform("size").values).tolist(),
            "control": cmap(cx, cy, None).tolist()}
res["fig"] = F
out.write_text(json.dumps(res, indent=1))
print({k: v for k, v in res.items() if k != "fig"})
