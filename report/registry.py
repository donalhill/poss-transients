"""The numbers the page quotes: each key with the output it comes from, its place in that output, its format, and the
page section it belongs to. The copy refers to them as {n:key}; nothing numerical is typed into the copy. Writes
report/numbers.json: {"numbers": {key: text}, "sources": {key: "file#path"}}.
usage: python report/registry.py"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "release"
VS = R / "vasco_rated"                                                       # the tests' sample (build_sample.py)
FILES = {"fit": VS / "shadow_fit.json", "perm": VS / "permutation_null.json",  # the shadow test on the real plates; shuffles
         "rep": VS / "replicas.json",                                           # its null and limit from replicas of the shadows
         "total": VS / "shadow_total.json",                                     # the result with the plate times' systematic
         "geo": VS / "shadow_geometry.json",
         "nb": VS / "shadow_neighbours.json",                                   # the sky around the near-empty shadow plates
         "alfit": VS / "alignment_fit.json", "align": VS / "alignments.json",
         "alchecks": VS / "alignment_checks.json", "alexample": VS / "alignment_example.json",
         "nuc": VS / "nuclear.json", "nucsym": VS / "nuclear_window-1_0_1.json",
         "dates": ROOT / "report/dates_stats.json",
         "cat": ROOT / "report/vasco_stats.json",                               # VASCO's catalogue itself
         "repeat": VS / "vasco_repeat.json",                                       # VASCO's shadow count repeated step by step
         "flat": R / "vasco_flat/replicas.json",                               # the shadow test with a flat map
         "lowdec": R / "vasco_lowdec/replicas.json",                              # with the map from low-declination plates
         "raw": R / "vasco_raw/replicas.json",                                    # with each transient counted once
         "site": VS / "site_data.json",                                           # the illustration of the shadow test's model
         "ratings": VS / "vasco_ratings.json",                                     # VASCO-ML retrained
         "nucvlist": VS / "nuclear_vascolist.json"}                                 # the nuclear test with VASCO's list of tests

# key, file, path (keys separated by "/", list items by number), format, section of the page
REG = [
    ("dates.by_log", "dates", "times/date_source/Palomar_log", "int", 3),
    ("dates.by_label", "dates", "times/date_source/label", "int", 3),
    ("dates.by_asteroids", "dates", "times/date_source/asteroids", "int", 3),
    ("dates.early_pct", "dates", "early_pct", "plain", 3),
    ("dates.ast_late_trails", "dates", "asteroid_late_trails", "int", 3),
    ("dates.ast_late_plates", "dates", "asteroid_late_plates", "int", 3),
    ("dates.ast_ok_trails", "dates", "asteroid_ok_trails", "int", 3),
    ("dates.ast_ok_plates", "dates", "asteroid_ok_plates", "int", 3),
    ("dates.log_checked", "dates", "log_checked", "int", 3),
    ("dates.log_agree", "dates", "log_agree", "int", 3),
    ("dates.exceptions", "dates", "exceptions", "words", 3),
    ("dates.gross", "dates", "gross", "words", 3),
    ("dates.one_day", "dates", "one_day", "words", 3),
    ("times.own_asteroids", "dates", "times/calibrated_own_asteroids", "int", 3),
    ("times.survey", "dates", "times/calibrated_survey_offset", "int", 3),
    ("times.late", "dates", "times/survey_offset_min", "f1abs", 3),
    ("times.errors_n", "dates", "times/shared_record_errors/plates", "words", 3),
    ("times.errors_min", "dates", "times/shared_record_errors/min_abs_min", "int", 3),
    ("times.errors_max", "dates", "times/shared_record_errors/max_abs_min", "int", 3),
    ("times.precision", "dates", "times/time_halfwidth68_median_min", "f1", 3),
    ("times.apart", "dates", "times/header_log_apart", "words", 3),
    ("times.trail_sd", "dates", "times/trail_scatter_min", "f1", 3),
    ("times.late_lo", "dates", "times/year_offset_min/0", "f1abs", 3),
    ("times.late_hi", "dates", "times/year_offset_min/1", "f1abs", 3),
    ("sh.transients", "fit", "fit/n", "int", 5),
    ("sh.in_shadow", "fit", "fit/O", "f1", 5),
    ("sh.obs", "fit", "fit/O", "f1", 5),
    ("sh.exp", "fit", "fit/E", "f1", 5),
    ("sh.oe", "fit", "fit/r", "f2", 5),
    ("sh.z", "fit", "fit/Z", "f1", 5),
    ("sh.replicas", "rep", "replicas", "int", 5),
    ("sh.null_median", "rep", "matches/count_dec/null/OE/median", "f2", 5),
    ("sh.null_q16", "rep", "matches/count_dec/null/OE/q16", "f2", 5),
    ("sh.null_q84", "rep", "matches/count_dec/null/OE/q84", "f2", 5),
    ("sh.p_t", "total", "p", "times", 5),
    ("sh.perm_p_t", "perm", "p_near/OE", "times", 5),
    ("sh.limit95", "total", "limit95", "pct0", 5),
    ("sh.limit_count", "rep", "matches/count/limit95", "pct0", 5),
    ("sh.limit_pointing", "rep", "matches/pointing/limit95", "pct0", 5),
    ("sh.sys", "total", "oe_sd_time", "f2", 5),
    ("sh.flat.p_t", "flat", "matches/count_dec/p/OE", "times", 5),
    ("sh.raw.p_t", "raw", "matches/count_dec/p/OE", "times", 5),
    ("sh.raw.oe", "raw", "real/OE", "f2", 5),
    ("sh.ratio_fit", "fit", "fit/ratio_fit", "f2", 5),
    ("sh.share_fit", "fit", "fit/glint_share_fit", "pct0", 5),
    ("est.e1_r", "fit", "fit/E1", "f2", 5),
    ("est.e2_r", "fit", "fit/E2", "f2", 5),
    ("est.e1_null_n", "rep", "matches/count_dec/null/E1/median", "f2", 5),
    ("est.e1_null_mu", "rep", "matches/pointing/null/E1/median", "f2", 5),
    ("est.power_g0", "rep", "power_g0", "pct0", 5),
    ("est.z_power", "rep", "matches/count_dec/power/OE", "pct0", 5),
    ("est.e1_power", "rep", "matches/count_dec/power/E1", "pct0", 5),
    ("est.e2_power", "rep", "matches/count_dec/power/E2", "pct0", 5),
    ("est.h_r", "fit", "fit/hybrid", "f2", 5),
    ("est.nb_radius", "nb", "radius_deg", "intstr", 5),
    ("est.nb_median", "nb", "neighbour_median", "pct0", 5),
    ("est.nb_share", "nb", "neighbour_share_low", "pct0", 5),
    ("est.survey_share", "nb", "survey_share_low", "pct0", 5),
    ("est.h_power", "rep", "matches/count_dec/power/H", "pct0", 5),
    ("sh.fit_p_t", "rep", "matches/count_dec/p/r", "times", 5),
    ("model.art", "site", "model/artefacts", "int", 5),
    ("model.glints", "site", "model/glints", "int", 5),
    ("model.share", "site", "model/share", "pct0", 5),
    ("model.lost", "site", "model/lost", "int", 5),
    ("model.oe0", "site", "model/oe_no_glints", "f2", 5),
    ("sh.geo.umbra", "geo", "umbra_deg", "f2", 5),
    ("sh.geo.penumbra", "geo", "penumbra_deg", "f2", 5),
    ("sh.geo.shift_lo", "geo", "shift_lo", "f1", 5),
    ("sh.geo.shift_hi", "geo", "shift_hi", "f1", 5),
    ("sh.geo.sweep_lo", "geo", "sweep_lo", "f1", 5),
    ("sh.geo.sweep_hi", "geo", "sweep_hi", "f1", 5),
    ("sh.geo.plates", "geo", "plates", "int", 5),
    ("sh.geo.whole", "geo", "whole", "words", 5),
    ("sh.geo.thin", "geo", "thin", "words", 5),
    ("sh.geo.partial", "geo", "partial", "int", 5),
    ("sh.geo.worth", "geo", "plates_worth", "f1", 5),
    ("sh.geo.rich_plate", "geo", "rich/plate", "plain", 5),
    ("sh.geo.rich_n", "geo", "rich/n", "int", 5),
    ("sh.geo.poor_plate", "geo", "poor/plate", "plain", 5),
    ("sh.geo.poor_n", "geo", "poor/n", "int", 5),
    ("al.p_4_t", "alfit", "flashes/4/p", "times", 6),
    ("al.copies", "alfit", "copies", "int", 6),
    ("al.shift_arcsec", "align", "shift_arcsec", "int", 6),
    ("al.link_arcsec", "align", "link_arcsec", "int", 6),
    ("al.transients", "alchecks", "transients_searched", "int", 6),
    ("al.dir_tracks", "alchecks", "direction/tracks_per_angle", "int", 6),
    ("al.dir_angles", "alchecks", "direction/angles", "int", 6),
    ("al.dir_min", "alchecks", "direction/found_min", "pct0", 6),
    ("al.dir_max", "alchecks", "direction/found_max", "pct0", 6),
    ("al.absorb_tracks", "alchecks", "absorption/tracks", "int", 6),
    ("al.absorb4", "alchecks", "absorption/per_100_1e-4", "int", 6),
    ("al.absorb3", "alchecks", "absorption/per_100_1e-3", "int", 6),
    ("al.keep4", "alchecks", "absorption/excess_shown_1e-4", "pct0", 6),
    ("al.keep3", "alchecks", "absorption/excess_shown_1e-3", "pct0", 6),
    ("al.pairs20", "alchecks", "pairs/within_20_ratio", "f1", 6),
    ("al.close_kept3", "alchecks", "close_pairs/kept_1e-3", "f1", 6),
    ("al.close_split3", "alchecks", "close_pairs/split_1e-3", "f1", 6),
    ("al.sb_real3", "alchecks", "sideband/real_1e-3", "words", 6),
    ("al.sb_copies3", "alchecks", "sideband/copies_1e-3", "f1", 6),
    ("al.sb_real2", "alchecks", "sideband/real_1e-2", "int", 6),
    ("al.sb_copies2", "alchecks", "sideband/copies_1e-2", "f1", 6),
    ("al.ex_plate", "alexample", "plate", "plain", 6),
    ("nuc.sym.plates_near", "nucsym", "plates_near", "int", 7),
    ("nuc.sym.ratio", "nucsym", "ratio", "f2", 7),
    ("nuc.sym.p_calendar_t", "nucsym", "calendar_null/p", "times", 7),
    ("nuc.per_plate_raw", "nuc", "per_photograph", "f2", 7),
    ("nuc.raw.ratio", "nuc", "per_photograph", "f2", 7),
    ("nuc.raw.p_calendar_t", "nuc", "per_photograph_calendar_null/p", "times", 7),
    ("nuc.raw.cal_median", "nuc", "per_photograph_calendar_null/median", "f2", 7),
    ("nuc.raw.exp_real", "nuc", "calendar_pointing/expected_real", "f2", 7),
    ("nuc.raw.corr", "nuc", "calendar_pointing/corr_log_raw_expected", "f2", 7),
    ("nuc.raw.share_exp_t", "nuc", "calendar_pointing/share_expected_as_high", "times", 7),
    ("nuc.raw.spring_real", "nuc", "calendar_pointing/spring_real", "pct0", 7),
    ("nuc.raw.spring_median", "nuc", "calendar_pointing/spring_median", "pct0", 7),
    ("nuc.ctx.per_night", "nuc", "context/per_night", "f2", 7),
    ("nuc.ratio", "nuc", "ratio", "f2", 7),
    ("nuc.plates", "nuc", "plates", "int", 7),
    ("nuc.plates_near", "nuc", "plates_near", "int", 7),
    ("nuc.calendars", "nuc", "calendars", "int", 7),
    ("nuc.held_out", "nuc", "held_out", "f2", 7),
    ("nuc.p_calendar_t", "nuc", "calendar_null/p", "times", 7),
    ("nuc.p_random_t", "nuc", "random_nights_null/p", "times", 7),
    ("nuc.inj15", "nuc", "injected_x1.5/ratio", "f2", 7),
    ("nuc.ctx.nights_near", "nuc", "context/nights_near", "int", 7),
    ("nuc.ctx.ppn_near", "nuc", "context/plates_per_night_near", "f1", 7),
    ("nuc.ctx.ppn_other", "nuc", "context/plates_per_night_other", "f1", 7),
    ("nuc.ctx.empty_near", "nuc", "context/empty_nights_near", "pct0", 7),
    ("nuc.ctx.empty_other", "nuc", "context/empty_nights_other", "pct0", 7),
    ("nuc.ctx.spring_near", "nuc", "context/spring_share_near", "pct0", 7),
    ("nuc.ctx.spring_other", "nuc", "context/spring_share_other", "pct0", 7),
    ("data.year_first", "cat", "years/0", "year", 1),
    ("data.year_last", "cat", "years/1", "year", 1),
    ("cat.transients", "cat", "transients", "int", 4),
    ("cat.plates", "cat", "plates", "int", 4),
    ("cat.min", "cat", "min", "int", 4),
    ("cat.max", "cat", "max", "int", 4),
    ("cat.avg", "cat", "avg", "int", 4),
    ("cat.noise", "cat", "noise", "int", 4),
    ("cat.band_lo", "cat", "band_lo", "int", 4),
    ("cat.band_hi", "cat", "band_hi", "int", 4),
    ("cat.in_band", "cat", "in_band", "int", 4),
    ("cat.over1000", "cat", "over1000", "words", 4),
    ("cat.edge_share", "cat", "edge_share", "pct0", 4),
    ("crowd.inner", "cat", "crowd_inner", "int", 4),
    ("crowd.outer", "cat", "crowd_outer", "int", 4),
    ("dust.n_dusty", "cat", "dust_n_dusty", "int", 4),
    ("dust.median_dusty", "cat", "dust_median_dusty", "int", 4),
    ("dust.n_clear", "cat", "dust_n_clear", "int", 4),
    ("dust.median_clear", "cat", "dust_median_clear", "int", 4),
    ("cat.eq_plates", "cat", "eq_plates", "int", 4),
    ("cat.eq_share", "cat", "eq_share", "pct0", 4),
    ("cat.eq_same_sky", "cat", "eq_same_sky", "pct0", 4),
    ("cat.pointing_share", "cat", "pointing_share", "pct0", 4),
    ("rep.vasco_obs", "repeat", "steps/0/in_shadow", "int", 5),
    ("rep.vasco_exp", "repeat", "steps/0/expected", "int", 5),
    ("rep.vasco_z", "repeat", "steps/0/binomial_z", "f1", 5),
    ("rep.rated_z", "repeat", "steps/4/binomial_z", "f1", 5),
    ("rep.vasco_r", "repeat", "steps/0/ratio", "f2", 5),
    ("rep.palomar_r", "repeat", "steps/1/ratio", "f2", 5),
    ("rep.times_r", "repeat", "steps/2/ratio", "f2", 5),
    ("rep.moving_r", "repeat", "steps/3/ratio", "f2", 5),
    ("rep.rated_r", "repeat", "steps/4/ratio", "f2", 5),
    ("rep.vasco_n", "repeat", "steps/0/transients", "int", 5),
    ("rep.shuffle_p", "repeat", "shuffled/rated/p_as_low", "times", 5),
    ("rep.shuffle_p_raw", "repeat", "shuffled/raw/p_as_low", "times", 5),
    ("al.lib_tracks", "alfit", "flashes/3/tracks", "int", 6),
    ("nuc.vlist.plates_near", "nucvlist", "plates_near", "int", 7),
    ("nuc.vlist.ratio", "nucvlist", "ratio", "f2", 7),
    ("nuc.vlist.p_calendar_t", "nucvlist", "calendar_null/p", "times", 7),
]

# method constants quoted in the text, read from the line of code that sets them (a change to the code breaks the text
# loudly): key, file, regular expression whose first group is the value, format
CONSTS = [
    ("c.shadow_diam", "analysis/shadow.py", r'^SHADOW_DEG = ([\d.]+)', "twice1"),
    ("c.al_lmax", "analysis/alignments.py", r'^LMAX = ([\d.]+)', "arcmin"),
    ("c.al_minsep", "analysis/alignments.py", r'^MINSEP = ([\d.]+)', "intstr"),
    ("c.al_box", "analysis/alignments.py", r'^BOX = ([\d.]+)', "intstr"),
    ("c.al_side", "analysis/alignments.py", r'side = max\(([\d.]+), 1\.5 \* L\)', "arcmin"),
    ("c.al_chi2p", "analysis/alignments.py", r'chi2\.isf\(([\d.]+), np\.arange', "pct0"),
]
SECTION_OF_CONST = {"c.shadow_diam": 5, "c.al_lmax": 6, "c.al_minsep": 6, "c.al_box": 6, "c.al_side": 6, "c.al_chi2p": 6}


def get(x, path):
    for k in path.split("/"):
        x = x[int(k)] if isinstance(x, list) else x[k]
    return x


def fmt(v, f):
    """The value as the page prints it (British style: thousands with commas, percentages whole unless stated)."""
    if f == "int":
        return f"{int(round(float(v))):,}"
    if f == "year":
        return str(int(v))
    if f in ("f1", "f2", "f3"):
        return f"{float(v):.{f[1]}f}"
    if f == "f1abs":
        return f"{abs(float(v)):.1f}"
    if f == "pct0":
        return f"{round(100 * float(v)) + 0.0:.0f}%"                       # + 0.0: no negative zero
    if f == "times":                   # a chance in words, as the page's charts give it: "1 time in 6", "2 times in 3"
        v = float(v)
        if v < 0.25:
            return f"1 time in {round(1 / v):,}"
        k, n = min(((k, n) for n in range(2, 6) for k in range(1, n)), key=lambda t: abs(t[0] / t[1] - v))
        return f"{k} {'time' if k == 1 else 'times'} in {n}"
    if f == "words":                   # small counts in words
        n = int(v); return ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"][n] if n < 10 else f"{n:,}"
    if f == "twice1":
        return f"{2 * float(v):.1f}"
    if f == "intstr":
        return f"{float(v):g}"
    if f == "arcmin":                  # arcseconds as arcminutes
        return f"{float(v) / 60:g}"
    return str(v)


def build():
    cache, numbers, sources = {}, {}, {}
    keys = [r[0] for r in REG]
    assert len(keys) == len(set(keys)), "a key registered twice"
    for key, f, path, form, _ in REG:
        p = FILES[f]
        if p not in cache:
            if not p.exists():
                raise SystemExit(f"{key}: {p.relative_to(ROOT)} is not yet produced")
            cache[p] = json.loads(p.read_text())
        sources[key] = f"{p.relative_to(ROOT)}#{path}"
        try:
            numbers[key] = fmt(get(cache[p], path), form)
        except (KeyError, IndexError, TypeError):
            raise SystemExit(f"{key}: no {path} in {p.relative_to(ROOT)}")
    for key, f, rx, form in CONSTS:
        m = re.search(rx, (ROOT / f).read_text(), re.M)
        if not m:
            raise SystemExit(f"{key}: no match for its pattern in {f}")
        numbers[key] = fmt(m.group(1), form); sources[key] = f"{f}#{rx}"
    return numbers, sources


def derive(numbers, sources):
    """Numbers computed from others. The share of transients more expected on the plates near tests from where they
    pointed; the crowding ratio in words; the range of the shadow ratio with each plate left out in turn; the pooled
    measurements without the two plates that hold under a tenth of their expected total; VASCO-ML retrained, the range
    over its sets; and the plate times' systematic in words when below 0.01. The text's worded facts are held here too:
    no fitted excess of 3, 5, or 6 flashes; one plate wholly in shadow; no shifted calendar as high as the 1.5x injection."""
    nu = json.loads(FILES["nuc"].read_text())
    numbers["nuc.ctx.sky_pct"] = fmt(nu["context"]["sky_expects_near"] - 1, "pct0"); sources["nuc.ctx.sky_pct"] = "nuc#context/sky_expects_near - 1"
    al = json.loads(FILES["alfit"].read_text())["flashes"]
    assert all(al[k]["s_hat"] == 0 for k in "356") and al["4"]["s_hat"] > 0, "the alignment verdict's wording"
    assert json.loads(FILES["geo"].read_text())["whole"] == 1, "the shadow text's 'one plate has no part ... sunlit'"
    assert nu["injected_x1.5"]["ratio"] > max(nu["calendar_draws"]), "the nuclear text's 'none of the shifted calendars'"
    # the pooled measurements without the plates entering the shadow test that hold under a tenth of their expected total
    F = json.loads(FILES["fit"].read_text()); Rp = json.loads(FILES["rep"].read_text())
    assert all(m["power"][k] < 0.5 for m in Rp["matches"].values() for k in ("OE", "E1", "E2", "H")), "the text's 'misses more often than it finds'"
    # the plate wholly in shadow: its total against the pointing model's; the scatter of plate totals about the model
    wh = [p for p in F["pooled_by_plate"] if p["area_in_shadow"] >= 0.999]
    assert len(wh) == 1, "the text's 'the plate wholly in shadow'"
    numbers["est.whole_plate"] = wh[0]["plate"]; numbers["est.whole_ratio"] = fmt(wh[0]["n"] / wh[0]["expected"], "f1")
    numbers["est.disp"] = fmt(F["fit"]["dispersion"] ** 0.5, "pct0")
    sources["est.whole_plate"] = sources["est.whole_ratio"] = "fit#pooled_by_plate (area_in_shadow = 1; n / expected)"
    sources["est.disp"] = "fit#fit/dispersion (square root: the scatter of plate totals about the pointing model, as a share)"
    pb = F["pooled_by_plate"]; low = [p for p in pb if not p["set_aside"] and p["n"] < 0.1 * p["expected"]]
    rest = [p for p in pb if p not in low]
    assert len(low) == 2, "the estimator text's 'two plates'"
    names = sorted(p["plate"] for p in low)
    numbers["est.low_plates"] = " and ".join(names) if len(names) == 2 else ", ".join(names[:-1]) + ", and " + names[-1]
    numbers["est.low_share"] = fmt(sum(p["n"] for p in low) / sum(p["expected"] for p in low), "pct0")
    numbers["est.low_exp_share"] = fmt(sum(p["expected_in_shadow"] for p in low) / sum(p["expected_in_shadow"] for p in pb), "pct0")
    numbers["est.e2_without"] = fmt(sum(p["in_shadow"] for p in rest) / sum(p["expected_in_shadow"] for p in rest), "f2")
    N_, J_ = F["totals"]["N"] - sum(p["n"] for p in low), F["totals"]["J"] - len(low)
    numbers["est.e1_without"] = fmt((sum(p["in_shadow"] for p in rest) / N_) / (sum(p["area_in_shadow"] for p in rest) / J_), "f2")
    for k in ("est.low_plates", "est.low_share", "est.low_exp_share", "est.e2_without", "est.e1_without"):
        sources[k] = "fit#pooled_by_plate (plates with n < 0.1 expected left out)"
    # VASCO-ML retrained: the range over the sets of labels and inputs of each set's median result
    if FILES["ratings"].exists():
        V = json.loads(FILES["ratings"].read_text())["sets"].values()
        for k, m, f in (("ml.nuc", "nuclear", "f2"), ("ml.nuc_p", "nuclear_p", "times"), ("ml.sh", "shadow_r", "f2")):
            vals = [v[m]["median"] for v in V]
            numbers[f"{k}_min"], numbers[f"{k}_max"] = fmt(min(vals), f), fmt(max(vals), f)
            sources[f"{k}_min"] = sources[f"{k}_max"] = f"ratings#sets/*/{m}/median (min, max)"
        assert len(V) == 8, "the text's eight settings"
        numbers["ml.sets"] = fmt(len(V), "words"); numbers["ml.runs"] = fmt(sum(v["runs"] for v in V), "int")
        sources["ml.sets"] = sources["ml.runs"] = "ratings#sets (count; runs summed)"
    if json.loads(FILES["total"].read_text())["oe_sd_time"] < 0.005:
        numbers["sh.sys"] = "less than 0.01"; sources["sh.sys"] = "total#oe_sd_time (< 0.005)"
    # each plate left out in turn: the plate whose absence moves the ratio most, and the range without each of the others
    numbers["sh.plates"] = fmt(len(F["plates"]), "words"); sources["sh.plates"] = "fit#plates (count)"
    ld_ = json.loads(FILES["lowdec"].read_text())["matches"]["count_dec"]["p"]["OE"]
    fl, ld = numbers["sh.flat.p_t"], fmt(ld_, "times")                     # the two maps' chances, in words
    numbers["sh.maps_p_t"] = f"in both cases about {fl}" if fl == ld else f"about {fl} and {ld} respectively"
    sources["sh.maps_p_t"] = "flat#matches/count_dec/p/OE, lowdec#matches/count_dec/p/OE (worded together)"
    meds = sorted(m["null"]["OE"]["median"] for m in Rp["matches"].values())   # our ratio's median with no effect, over the choices of host
    numbers["sh.null_med_lo"], numbers["sh.null_med_hi"] = fmt(meds[0], "f2"), fmt(meds[-1], "f2")
    sources["sh.null_med_lo"] = sources["sh.null_med_hi"] = "rep#matches/*/null/OE/median (min, max)"
    ps = sorted(m["p"]["E2"] for m in Rp["matches"].values())                 # the pooled ratio's chance over the choices of host
    numbers["est.e2_p_lo"], numbers["est.e2_p_hi"] = fmt(ps[-1], "times"), fmt(ps[0], "times")   # the likelier first
    sources["est.e2_p_lo"] = sources["est.e2_p_hi"] = "rep#matches/*/p/E2 (min, max)"
    for k_, s_ in (("e1", "E1"), ("h", "H")):
        ps = sorted(m["p"][s_] for m in Rp["matches"].values())
        numbers[f"est.{k_}_p_lo"], numbers[f"est.{k_}_p_hi"] = fmt(ps[-1], "times"), fmt(ps[0], "times")
        sources[f"est.{k_}_p_lo"] = sources[f"est.{k_}_p_hi"] = f"rep#matches/*/p/{s_} (min, max)"
    M = json.loads(FILES["site"].read_text())["model"]
    numbers["model.oe"] = fmt(M["oe"], "f2"); sources["model.oe"] = "site#model/oe"
    lo = {q: v["OE"] for q, v in F["leave_one_out"].items()}; top = max(lo, key=lambda q: abs(lo[q] - F["fit"]["r"]))
    others = [v for q, v in lo.items() if q != top]
    numbers["sh.loo.top"], numbers["sh.loo.top_oe"] = top, fmt(lo[top], "f2")
    numbers["sh.loo.oe_min"], numbers["sh.loo.oe_max"] = fmt(min(others), "f2"), fmt(max(others), "f2")
    for k in ("sh.loo.top", "sh.loo.top_oe", "sh.loo.oe_min", "sh.loo.oe_max"):
        sources[k] = "fit#leave_one_out/*/OE (the plate moving it most; the range without each other plate)"
    v = json.loads(FILES["cat"].read_text())
    numbers["crowd.ratio"] = fmt(round(v["crowd_outer"] / v["crowd_inner"]), "words"); sources["crowd.ratio"] = "cat#crowd_outer / crowd_inner"


if __name__ == "__main__":
    numbers, sources = build()
    derive(numbers, sources)
    out = ROOT / "report/numbers.json"
    out.write_text(json.dumps({"numbers": numbers, "sources": sources, "sections": {**{r[0]: r[4] for r in REG}, **SECTION_OF_CONST}},
                              indent=1, ensure_ascii=False))
    print(f"{out.relative_to(ROOT)}: {len(numbers)} numbers")
