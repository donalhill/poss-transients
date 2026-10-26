"""Every plate's mid-exposure time, for the plates of VASCO's catalogue, in two stages.

nominal: start and exposure from the date records in data/. The header's time is read as hours and hundredths of an
hour (06:82 is 06:49; Palomar's log confirms), on the date corrected by labels, asteroid dating, the log, or the block
rule (plate_date_corrections.csv); centres from the scan headers. Palomar's log (palomar_log_vi25.csv, the record kept at
the telescope) is matched to each plate by position (its pointing within 0.1 deg of the plate's centre), not by number:
the log's plate numbers run one off the labels in places, and fields were photographed more than once. The log names a
night by its evening, so its UT date is the next day; its start is Pacific Standard Time (UT - 8). Each plate takes the
matched entry nearest in time to its own start, within a day, else the matched entry with its label's number (a plate
whose date the log disputes: its time of day). The exposure is the header's, except where the log gives
the plate's red exposure and differs: then the log's; a header exposure under 10 min with no log entry is read as a
slipped decimal (x10). Where the header and the log give times of day more than 10 minutes apart, both are kept as
candidates. -> data/plate_times_nominal.csv

calibrate: the time of day from the sky (dates/asteroid_clock.py's trails, data/asteroid_clock.csv). A trail is clean
when it lies within 1" of its track, is elongated (> 2), and lasts 0.6 to 1.5 times the exposure. Each clean trail gives
the offset of the exposure's middle from the nominal middle; a plate's offsets agree to within the per-trail scatter,
measured on the plates with three or more clean trails (a robust standard deviation about each plate's median); trails
more than 15 minutes, or three times that scatter, from their plate's median are set aside. Where the
header and the log disagree, the candidate the plate's trails support (within 15 minutes) is taken; without trails, the
candidate consistent with the rest of the night: lying wholly in astronomical night (the Sun 18 deg below the
horizon), overlapping no other plate's exposure that night, and keeping the night's plates in the order of their
numbers. Then each plate's offset tau from its nominal middle is the posterior of a Gaussian prior, centred on its year's
offset (the median offset of that year's plates with two or more clean trails, where the year has at least YEAR_MIN of
them, else the survey's: the clock's offset changes by up to about half a minute from year to year, and not with the
hour of the night) with the spread of those plates' offsets about their centres (robust), and its clean trails, each
with the per-trail scatter: the time used is the posterior mean, its 68% range the posterior's, and 200 draws give the
plate times' systematic. A plate whose two or more trails agree on an offset more than GROSS minutes from its year's has
a wrong written time, which the prior does not describe: its time comes from its trails alone.
-> data/plate_times_all.csv (plate, ut, ra_c, dec_c, exp, time_source, exp_source), data/plate_time_posterior.csv
   (tau_median, tau_q16, tau_q84, from_prior_centre, on_line, ut), data/plate_time_draws.parquet (200 draws, minutes,
   one column per plate), data/asteroid_clock_matches.csv (plate, desig, rate, dt_min, perp, V: every clean trail)
usage: python dates/plate_times.py nominal|calibrate"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun, search_around_sky
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"
YEAR_MIN = 9                                                                 # timed plates a year needs for its own centre
GROSS = 3.0                                                                  # min from the year's offset: a wrong written time
PALOMAR = EarthLocation.from_geodetic(lon=-116.865 * u.deg, lat=33.356 * u.deg, height=1706 * u.m)


def nominal():
    PC = pd.read_parquet(D / "plates_corrected.parquet"); plates = PC.index          # the plates of VASCO's catalogue
    h = pd.read_parquet(D / "xe_headers.parquet").drop_duplicates("plate").set_index("plate")
    c = pd.read_csv(D / "plate_date_corrections.csv", index_col=0)
    x = pd.read_csv(D / "xe_plate_centres.csv").drop_duplicates("plate").set_index("plate")
    T = h[["date", "hh", "mm", "exp", "label"]].join(c[["corrected_date", "source"]]).join(x).reindex(plates)
    # the few plates with no scan header: their centre from VASCO's plates, their start from dates/date_corrections.py
    nohdr = T.hh.isna()
    T.loc[nohdr, "exp"] = PC.exp[nohdr]
    T["ra"], T["dec"] = T.ra.fillna(PC.ra_c), T.dec.fillna(PC.dec_c)                 # (and any plate without a scan centre)
    log = pd.read_csv(D / "palomar_log_vi25.csv").dropna(subset=["ra_j2000", "red_start_local", "evening_date"])
    hh_ = log.red_start_local.str[:2].astype(int); mm_ = log.red_start_local.str[3:5].astype(int)
    log["start"] = pd.to_datetime(log.evening_date) + pd.to_timedelta(hh_ + 8 + np.where(hh_ < 12, 24, 0), unit="h") + pd.to_timedelta(mm_, unit="min")
    log["log_exp"] = pd.to_numeric(log.red_exposure_min, errors="coerce")
    lc = SkyCoord(log.ra_j2000.values, log.dec_j2000.values, unit=(u.hourangle, u.deg))
    start = pd.to_datetime(T.corrected_date) + pd.to_timedelta(T.hh + T.mm / 100, unit="h")
    start[nohdr] = PC.ut[nohdr] - pd.to_timedelta(PC.exp[nohdr] / 2, unit="min")
    exp, esrc, alt = T.exp.astype(float).copy(), pd.Series("header", index=T.index), pd.Series(pd.NaT, index=T.index)
    ip, il, _, _ = search_around_sky(SkyCoord(T.ra.values * u.deg, T.dec.values * u.deg), lc, 0.1 * u.deg)   # every pair within 0.1 deg
    pairs = pd.DataFrame({"q": T.index.values[ip], "l": il})
    for q, grp in pairs.groupby("q"):
        near = log.iloc[grp.l.values]
        if pd.isna(start[q]):
            continue
        dt_ = (near.start - start[q]).abs(); same = pd.to_numeric(near.plate_number, errors="coerce") == T.label[q]
        if dt_.min() <= pd.Timedelta(days=1, hours=1):
            l = near.loc[dt_.idxmin()]
        elif same.any():                                                     # none within a day: the entry with the label's
            l = near[same].iloc[0]                                           # number (the plate's date is disputed)
        else:
            continue
        a = start[q].normalize() + (l.start - l.start.normalize())
        if abs(((a - start[q]).total_seconds() / 60 + 720) % 1440 - 720) > 10:
            alt[q] = a
        if pd.notna(l.log_exp) and l.log_exp != exp[q]:
            exp[q] = l.log_exp; esrc[q] = "Palomar_log"
    slip = (esrc == "header") & (exp < 10); exp[slip] *= 10; esrc[slip] = "header_x10"
    out = pd.DataFrame({"start": start, "log_start": alt, "exp": exp, "exp_source": esrc, "ra_c": T.ra, "dec_c": T.dec,
                        "label": T.label, "date_source": T.source.fillna("unresolved").str.replace(" ", "_"), "header": ~nohdr}, index=T.index)
    out.index.name = "plate"
    out.to_csv(D / "plate_times_nominal.csv")
    print(f"{len(out)} plates; {int(out.log_start.notna().sum())} with header and log more than 10 minutes apart")


def clean(A, exp):
    """The clean trails: on the track, elongated, and as long as the exposure. The cuts are where the trails' scatter about
    their plate's median time rises (all trails measured, 2026-10-26: about 1.0 min within them; 2.4 min and more beyond
    0.7 arcsec from the track, 1.5 min and more at elongation below 2.5, 2 min and more outside 0.8 to 1.2 exposures)."""
    e = A.plate.map(exp)
    return A[(A.off_track < 0.7) & (A.elong > 2.5) & A.trail_min.between(0.8 * e, 1.2 * e)].copy()


def dark(t):
    """True where the Sun is more than 18 deg below Palomar's horizon at UT t (astronomical night)."""
    return get_sun(Time(t)).transform_to(AltAz(obstime=Time(t), location=PALOMAR)).alt.deg < -18


def night_rule(q, N, start):
    """For a plate whose header and log disagree and whose trails do not decide: the candidate that lies wholly in
    astronomical night, overlaps no other exposure of its night, and keeps the night's plates in the order of their
    numbers; the header's when both or neither do."""
    same = [p for p in N.index if p != q and pd.notna(start[p]) and abs((start[p] - start[q]).total_seconds()) < 14 * 3600]
    ok = {}
    for k, s in (("header", N.start[q]), ("log", N.log_start[q])):
        e = s + pd.Timedelta(minutes=N.exp[q]); good = bool(dark(s.to_pydatetime()) and dark(e.to_pydatetime()))
        for p in same:
            sp, ep = start[p], start[p] + pd.Timedelta(minutes=N.exp[p])
            good &= not (s < ep and sp < e)
            lab_q, lab_p = pd.to_numeric(N.label[q], errors="coerce"), pd.to_numeric(N.label[p], errors="coerce")
            if pd.notna(lab_q) and pd.notna(lab_p):
                good &= (lab_p < lab_q) == (sp < s)
        ok[k] = good
    return "log" if ok["log"] and not ok["header"] else "header"


def calibrate():
    N = pd.read_csv(D / "plate_times_nominal.csv", index_col=0, parse_dates=["start", "log_start"])
    A = pd.read_csv(D / "asteroid_clock.csv", parse_dates=["t_mid"])
    A = clean(A[A.plate.isin(N.index) & (A.status == "ok")], N.exp)
    mid = N.start + pd.to_timedelta(N.exp / 2, unit="min")
    A["dt_min"] = (A.t_mid - A.plate.map(mid)).dt.total_seconds() / 60
    # the plates where header and log disagree: the candidate the trails support, else the night's rule
    start, tsrc = N.start.copy(), "header+" + N.date_source
    for q in N.index[N.log_start.notna()]:
        d = A.dt_min[A.plate == q]; shift = (N.log_start[q] - N.start[q]).total_seconds() / 60
        if len(d) and abs(d.median() - shift) < 15 and abs(d.median()) >= 15:
            pick = "log"
        elif len(d) and abs(d.median()) < 15:
            pick = "header"
        else:
            pick = night_rule(q, N, start)
        if pick == "log":
            start[q] = N.log_start[q]; tsrc[q] += "+log_time"
            A.loc[A.plate == q, "dt_min"] -= shift
    mid = start + pd.to_timedelta(N.exp / 2, unit="min")
    # within each plate, trails more than 15 minutes from the plate's median are set aside (a faint star on the track
    # taken for the trail), then the per-trail scatter, and trails beyond 3 times it from their plate's median
    A = A[(A.dt_min - A.plate.map(A.groupby("plate").dt_min.median())).abs() < 15]
    g = A.groupby("plate").dt_min
    many = g.size()[g.size() >= 3].index
    dev = A[A.plate.isin(many)].dt_min - A[A.plate.isin(many)].plate.map(g.median())
    A = A[(A.dt_min - A.plate.map(g.median())).abs() < 3 * 1.4826 * float(np.median(np.abs(dev)))]
    # the per-trail scatter, the survey offset and its spread
    g = A.groupby("plate").dt_min
    many = g.size()[g.size() >= 3].index
    dev = A[A.plate.isin(many)].dt_min - A[A.plate.isin(many)].plate.map(g.median())
    sig_t = 1.4826 * float(np.median(np.abs(dev)))
    two = g.size()[g.size() >= 2].index
    S = float(g.median()[two].median()); year = mid.dt.year
    byyear = g.median()[two].groupby(year[two]).agg(["median", "size"])
    centre = year.map(byyear["median"].where(byyear["size"] >= YEAR_MIN)).fillna(S)
    sig_s = 1.4826 * float(np.median(np.abs(g.median()[two] - centre[two])))
    # each plate's posterior: Gaussian prior N(its year's centre, sig_s), each clean trail N(tau, sig_t)
    n = g.size().reindex(N.index).fillna(0).astype(int); sm = g.sum().reindex(N.index).fillna(0.0)
    prec = 1 / sig_s ** 2 + n / sig_t ** 2
    tau = (centre / sig_s ** 2 + sm / sig_t ** 2) / prec; sd = 1 / np.sqrt(prec)
    # a plate whose trails agree on an offset far from its year's (GROSS minutes, two or more trails) has a written time
    # that is wrong: the prior's narrow spread does not describe it, and its time is taken from its trails alone
    gross = (n >= 2) & ((g.median().reindex(N.index) - centre).abs() > GROSS)
    tau[gross] = (sm / n.clip(lower=1))[gross]; sd[gross] = (sig_t / np.sqrt(n.clip(lower=1)))[gross]
    rng = np.random.default_rng(2026)
    draws = pd.DataFrame(rng.standard_normal((200, len(N))) * sd.values[None, :], columns=N.index)
    own = (n > 0) & ((tau - centre).abs() > 0.5)
    tsrc = tsrc + np.where(own, "+own_asteroids", "+survey_clock")
    ut = mid + pd.to_timedelta(tau, unit="min")
    P = pd.DataFrame({"on_line": n, "tau_median": tau.round(2), "tau_q16": (tau - sd).round(2), "tau_q84": (tau + sd).round(2),
                      "from_prior_centre": (tau - centre).round(2), "ut": ut}, index=N.index)
    T = pd.DataFrame({"ut": ut, "ra_c": N.ra_c, "dec_c": N.dec_c, "exp": N.exp, "time_source": tsrc, "exp_source": N.exp_source}, index=N.index)
    for f_, df in (("plate_times_all.csv", T), ("plate_time_posterior.csv", P)):
        df.index.name = "plate"; df.to_csv(D / f_)
    draws.round(3).to_parquet(D / "plate_time_draws.parquet")
    A.assign(desig=A.number, rate=A.rate_as_min * 60 / 3600 * 60, perp=A.off_track)[["plate", "desig", "rate", "dt_min", "perp", "V"]] \
        .to_csv(D / "asteroid_clock_matches.csv", index=False)
    print(f"{len(A)} clean trails on {int((n > 0).sum())} plates; per-trail scatter {sig_t:.2f} min; survey offset {S:+.2f} min, by year "
          f"{byyear['median'].where(byyear['size'] >= YEAR_MIN).dropna().round(2).to_dict()} (spread {sig_s:.2f}); {int(own.sum())} plates timed "
          f"by their own trails; {int(tsrc.str.contains('log_time').sum())} take the log's time")


if __name__ == "__main__":
    {"nominal": nominal, "calibrate": calibrate}[sys.argv[1]]()
