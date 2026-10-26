# POSS-I transients: VASCO's catalogue (Solano et al. 2022; the 107,875 transients of their shadow test, from the
# VASCO-ML repository), every plate on its corrected date and asteroid-calibrated time, through the three tests, and the
# page. Targets are files, so `make` reruns only what changed and `make -n` shows the chain.
#   make tests        the samples, the three tests, and their checks
#   make report       the page in docs/, its numbers held to the copy (report/check.py)

PY       ?= python3
PYML     ?= $(PY)
PYEPH    ?= $(PY)
D        := data
DUST     := data/external/dustmaps/sfd/SFD_dust_4096_ngp.fits
ANALYSIS := analysis/sample.py analysis/shadow.py analysis/geometry.py
SHADOW   := analysis/shadow_fit.py analysis/shadow_stats.py analysis/glint_fit.py analysis/expected.py $(DUST) $(ANALYSIS)
NUCLEAR  := analysis/nuclear.py analysis/expected.py $(ANALYSIS) $(DUST) $(D)/us_shot_times.csv $(D)/palomar_log_vi25.csv
ALIGN    := analysis/alignments.py $(ANALYSIS)
EXT      := $(D)/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv
ML       := $(D)/external/scored_catalog.zip
LOG      := $(D)/palomar_log_vi25.csv
EVIDENCE := $(D)/label_readings.csv $(D)/extra_labels.csv $(D)/date_check_summary.csv $(D)/xe_headers.parquet $(D)/xe_plate_centres.csv
TIME_DRAWS := $(D)/plate_time_draws.parquet

.PHONY: tests report plate-times
.SECONDARY:                                # keep every intermediate output (the samples' crossed.json among them)

# VASCO's public inputs, fetched once from fixed commits (checksums in vasco/fetch_external.py), Palomar's log, and the
# SFD dust map
$(EXT) $(ML):
	$(PY) vasco/fetch_external.py
$(LOG):
	$(PY) dates/fetch_palomar_log.py
$(DUST):
	$(PY) -c "from dustmaps.config import config; config['data_dir'] = 'data/external/dustmaps'; import dustmaps.sfd; dustmaps.sfd.fetch()"

# each transient on its plate (and every plate's centre and header time), and the plates' corrected dates (the times
# then come from data/plate_times_all.csv)
$(D)/transients_with_plates.parquet: vasco/plates.py analysis/geometry.py $(EXT)
	$(PY) vasco/plates.py
$(D)/plates.parquet: $(D)/transients_with_plates.parquet ;
$(D)/plates_corrected.parquet: dates/date_corrections.py dates/date_check.py $(D)/transients_with_plates.parquet $(LOG) $(EVIDENCE)
	$(PY) dates/date_corrections.py
# each plate's time of exposure, calibrated against its asteroids' trails (slow: the known asteroids on each plate from
# their orbits, under a python with ASSIST, PYEPH; JPL Horizons' tracks; boxes read from each plate's scan and Gaia):
# kept in the repository, built only when missing, or remade with `make plate-times`
EPHEM := $(D)/external/ephemerides
$(EPHEM)/de440.bsp:
	mkdir -p $(@D) && curl -sL -o $@ https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp
$(EPHEM)/sb441-n16.bsp:
	mkdir -p $(@D) && curl -sL -o $@ https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp
$(EPHEM)/MPCORB.DAT.gz:
	mkdir -p $(@D) && curl -sL -o $@ https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz
$(D)/plate_times_all.csv:
	$(MAKE) plate-times
plate-times: $(D)/plates_corrected.parquet $(EPHEM)/de440.bsp $(EPHEM)/sb441-n16.bsp $(EPHEM)/MPCORB.DAT.gz
	$(PY) dates/plate_times.py nominal
	$(PYEPH) dates/asteroid_list.py
	$(PY) dates/asteroid_clock.py 12
	$(PY) dates/plate_times.py calibrate
$(D)/plate_time_posterior.csv $(D)/plate_time_draws.parquet $(D)/asteroid_clock_matches.csv: $(D)/plate_times_all.csv ;
# each plate's depth (slow: DSS cutouts and Gaia) and dust (IRSA's dust service): kept in the repository, built only
# when missing
$(D)/plate_depth.csv:
	$(PY) vasco/plate_depth.py
$(D)/plate_dust.csv:
	$(PY) vasco/plate_dust.py
# VASCO-ML retrained with VASCO's own training code, four sets of labels with and without the per-night inputs (slow;
# under VASCO-ML's pinned versions: PYML, a python with vasco/requirements-vasco-ml.txt)
RATINGS := $(D)/external/vasco_ml/ratings/doc250_no_window_without_3.npy
$(RATINGS): vasco/ml_retrain.py
	$(PYML) vasco/ml_retrain.py

# the samples: transients by VASCO-ML's rating (the tests' sample, VS), each counted once (VR), with a flat map in place
# of the edge map (VF), and with the edge map from the host plates below dec 25 (VL)
VS := release/vasco_rated
VR := release/vasco_raw
VF := release/vasco_flat
VL := release/vasco_lowdec
VASCO_IN := build_sample.py vasco/vasco_sample.py analysis/shadow.py $(D)/transients_with_plates.parquet $(D)/plates_corrected.parquet $(D)/plate_times_all.csv $(D)/plate_depth.csv $(ML)
$(VS)/plates.parquet: $(VASCO_IN)
	$(PY) build_sample.py $(VS) > $(VS).build.log
$(VR)/plates.parquet: $(VASCO_IN)
	$(PY) build_sample.py $(VR) --weights none > $(VR).build.log
$(VF)/plates.parquet: $(VASCO_IN)
	$(PY) build_sample.py $(VF) --acceptance flat > $(VF).build.log
$(VL)/plates.parquet: $(VS)/plates.parquet diagnostics/acceptance_lowdec.py
	$(PY) diagnostics/acceptance_lowdec.py $(VS) $(VL) > $(VL).build.log

# the plates the shadow can reach, for any sample
release/%/crossed.json: release/%/plates.parquet analysis/crossed.py $(ANALYSIS)
	POSS_SAMPLE=$(@D) $(PY) analysis/crossed.py > $(@D)/crossed.log

# Earth's shadow: the statistics on the real plates, their null and the limit on glints from replicas of the real
# shadows on other plates, the cross-check with shuffled plate times, the plate times drawn from their posteriors, the
# result with the plate times' systematic included, the shadow's geometry, and the plate-level data for the charts; the
# checks of the edge map
release/%/shadow_fit.json: release/%/crossed.json $(SHADOW)
	POSS_SAMPLE=$(@D) $(PY) -u analysis/shadow_fit.py > $(@D)/shadow_fit.log
release/%/replicas.json: release/%/crossed.json analysis/replicas.py $(SHADOW)
	POSS_SAMPLE=$(@D) $(PY) -u analysis/replicas.py 20000 > $(@D)/replicas.log
%/permutation_null.json: %/crossed.json analysis/permutation_null.py $(SHADOW)
	POSS_SAMPLE=$* $(PY) -u analysis/permutation_null.py 10000 > $*/permutation_null.log
$(VS)/time_posterior/shadow_fit_only.json: $(VS)/shadow_fit.json $(TIME_DRAWS)
	mkdir -p $(@D) && for f in $(VS)/*.parquet $(VS)/*.json; do [ -e $(@D)/$$(basename $$f) ] || ln -s ../$$(basename $$f) $(@D)/; done
	POSS_FIT_ONLY=1 POSS_TIME_JITTER_FILE=$(TIME_DRAWS) POSS_TIME_JITTER_N=100 POSS_SAMPLE=$(@D) $(PY) -u analysis/shadow_fit.py > $(@D)/run.log
$(VS)/shadow_total.json: $(VS)/replicas.json $(VS)/time_posterior/shadow_fit_only.json diagnostics/shadow_total.py
	$(PY) diagnostics/shadow_total.py $(VS)
$(VS)/shadow_neighbours.json: $(VS)/shadow_fit.json diagnostics/shadow_neighbours.py $(SHADOW)
	POSS_SAMPLE=$(VS) $(PY) diagnostics/shadow_neighbours.py > /dev/null
$(VS)/shadow_geometry.json: $(VS)/shadow_fit.json analysis/shadow_geometry.py
	POSS_SAMPLE=$(VS) $(PY) analysis/shadow_geometry.py > /dev/null
$(VS)/site_data.json: $(VS)/shadow_fit.json analysis/site_data.py
	POSS_SAMPLE=$(VS) $(PY) analysis/site_data.py

# VASCO's shadow count repeated, then changed one step at a time to the pooled measurement against the control positions
$(VS)/vasco_repeat.json: $(VS)/permutation_null.json $(VR)/permutation_null.json diagnostics/vasco_repeat.py vasco/vasco_sample.py $(SHADOW) \
                         $(D)/external/XE_final_random_sources.csv $(D)/external/poss_i_observing_dates_strict.csv $(D)/plates.parquet
	POSS_SAMPLE=$(VS) $(PY) diagnostics/vasco_repeat.py $(VR) > $(VS)/vasco_repeat.log
# the rated nuclear and shadow tests on each set of retrained ratings
$(VS)/vasco_ratings.json: $(VS)/crossed.json diagnostics/vasco_ratings.py analysis/expected.py $(SHADOW) $(DUST) $(RATINGS)
	POSS_SAMPLE=$(VS) $(PY) -u diagnostics/vasco_ratings.py > $(VS)/vasco_ratings.log

# the nuclear tests: VASCO's window, the symmetric one (a day either side in UT), and VASCO's own list of US tests
# (without Operation Upshot-Knothole)
$(VS)/nuclear.json: $(VS)/crossed.json $(NUCLEAR)
	POSS_SAMPLE=$(VS) $(PY) analysis/nuclear.py > $(VS)/nuclear.log
$(VS)/nuclear_window-1_0_1.json: $(VS)/crossed.json $(NUCLEAR)
	POSS_NUC_WINDOW=-1,0,1 POSS_SAMPLE=$(VS) $(PY) analysis/nuclear.py > $(VS)/nuclear_window.log
$(VS)/nuclear_vascolist.json: $(VS)/crossed.json $(NUCLEAR) $(D)/external/nuclear_tests_1949_1957.csv
	POSS_NUC_LIST=vasco POSS_SAMPLE=$(VS) $(PY) analysis/nuclear.py > $(VS)/nuclear_vascolist.log

# the alignment search, against 1,000 random catalogues built from the data; the fit, one process per number of
# flashes, from tracks added to the catalogue; the checks on the random catalogues (direction, what they keep of a
# track, the sideband), and the example candidate for the method's chart
$(VS)/alignments.npz: $(VS)/plates.parquet $(ALIGN)
	POSS_SAMPLE=$(VS) $(PY) analysis/alignments.py 1000 > $(VS)/alignments.log
$(VS)/alignment_sideband.json: $(VS)/plates.parquet analysis/alignment_sideband.py $(ALIGN)
	POSS_SAMPLE=$(VS) $(PY) -u analysis/alignment_sideband.py > $(VS)/alignment_sideband.log
$(VS)/alignment_direction.json: $(VS)/plates.parquet analysis/alignment_direction.py $(ALIGN)
	POSS_SAMPLE=$(VS) $(PY) -u analysis/alignment_direction.py > $(VS)/alignment_direction.log
$(VS)/alignment_null_tune.json: $(VS)/plates.parquet analysis/alignment_null_tune.py $(ALIGN)
	POSS_SAMPLE=$(VS) $(PY) -u analysis/alignment_null_tune.py > $(VS)/alignment_null_tune.log
$(VS)/alignment_checks.json: $(VS)/alignment_direction.json $(VS)/alignment_null_tune.json $(VS)/alignment_sideband.json analysis/alignment_checks_summary.py
	POSS_SAMPLE=$(VS) $(PY) analysis/alignment_checks_summary.py > /dev/null
$(VS)/alignment_fit.json: $(VS)/alignments.npz analysis/alignment_fit.py analysis/alignment_fit_merge.py
	for f in 3 4 5 6; do POSS_FIT_FLASHES=$$f POSS_SAMPLE=$(VS) $(PY) -u analysis/alignment_fit.py > $(VS)/alignment_fit_$$f.log & done; wait
	POSS_SAMPLE=$(VS) $(PY) analysis/alignment_fit_merge.py > $(VS)/alignment_fit.log
$(VS)/alignment_example.json: $(VS)/alignments.npz analysis/alignment_example.py
	POSS_SAMPLE=$(VS) $(PY) analysis/alignment_example.py > /dev/null

tests: $(addprefix $(VS)/,shadow_fit.json replicas.json permutation_null.json shadow_total.json shadow_geometry.json shadow_neighbours.json site_data.json \
         alignment_fit.json alignment_checks.json alignment_example.json nuclear.json nuclear_window-1_0_1.json nuclear_vascolist.json \
         vasco_repeat.json vasco_ratings.json) \
       $(addprefix $(VR)/,permutation_null.json shadow_fit.json replicas.json) $(addprefix $(VF)/,shadow_fit.json replicas.json) \
       $(addprefix $(VL)/,shadow_fit.json replicas.json)

# the numbers the page quotes on VASCO's catalogue and on the plates' dates and times
report/vasco_stats.json: report/vasco_stats.py $(VS)/plates.parquet $(D)/plate_dust.csv analysis/expected.py vasco/vasco_sample.py $(D)/external/XE_final_random_sources.csv analysis/geometry.py
	$(PY) report/vasco_stats.py $(VS) $@ > /dev/null
report/dates_stats.json: report/dates_stats.py $(D)/plates_corrected.parquet $(D)/plate_date_corrections.csv $(D)/date_check_summary.csv $(D)/plate_times_all.csv $(D)/plate_time_posterior.csv \
                         $(D)/plate_times_nominal.csv $(D)/asteroid_clock_matches.csv dates/plate_times.py
	$(PY) report/dates_stats.py data $(D)/plates_corrected.parquet $@

# the page: the registry's numbers and the figures' data, the text by section (report/pages.py), the template from the
# briefs, the page assembled with its figure references, its provenance, and the check that holds the copy to the numbers
STATS := report/vasco_stats.json report/dates_stats.json
TEXT := $(wildcard report/briefs/*.md report/copy/*.md) report/pages.py report/briefs.head.md report/copy.head.md
report/numbers.json: report/registry.py tests $(STATS)
	$(PY) report/registry.py
report/figures.json: report/figures.py tests $(STATS) $(D)/plate_date_corrections.csv $(D)/asteroid_clock_matches.csv $(D)/plate_time_posterior.csv
	$(PY) report/figures.py
report/copy.md: $(TEXT)
	$(PY) report/pages.py
report/briefs.md: report/copy.md ;
report/template.html: report/template.py report/briefs.md report/copy.md report/base_template.html
	$(PY) report/template.py
report: report/numbers.json report/figures.json report/template.html
	$(PY) report/check.py
	$(PY) report/assemble.py report/copy.md report/template.html report/numbers.json docs/index.html
	$(PY) report/figrefs.py docs/index.html
	cp report/figures.json docs/figures.json
	$(PY) report/provenance.py
