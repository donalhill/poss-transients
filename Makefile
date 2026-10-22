# The whole analysis, from VASCO's public catalogue to the page. `make` builds what is out of date; `make -n` shows
# what would run; `make help` lists the stages. Each stage is one script; see README.md.

PY ?= python
A := analysis
D := data

EXT      := $(D)/external/SUPERVIKTIG_HELAVASCO_validated_v4.csv
LOG      := $(D)/palomar_log_vi25.csv
PLATES   := $(D)/transients_with_plates.parquet
DATES    := $(D)/plates_corrected.parquet
DUST     := $(D)/plate_dust.csv
DEPTH    := $(D)/plate_depth.csv
POOLED   := $(D)/shadow_pooled.json
SOUTH    := $(D)/shadow_plates_log.csv
FIT      := $(D)/shadow_fit.json
NUCLEAR  := $(D)/nuclear.json
PROPS    := $(A)/plate_properties.py $(DUST) $(DEPTH) $(D)/xe_headers.parquet
ALIGN    := $(D)/alignments.json
VCAND    := $(D)/vasco_candidates.json
AFIT     := $(D)/alignment_fit.json
NULLCHK  := $(D)/alignment_null_checks.json
SHCHK    := $(D)/shadow_null_checks.json
TESTSV   := $(D)/external/nuclear_tests_1949_1957.csv $(D)/external/poss_i_observing_dates_strict.csv
ML       := $(D)/external/scored_catalog.zip
LABELS   := $(D)/external/combined_training_labels.csv
MLCHK    := $(D)/ml_checks.json
MLRT     := $(D)/ml_retrain.json
PYML     ?= $(PY)
SITE     := docs/data/site.json
PAGE     := docs/index.html

# inputs read off the plates by eye or from slow public services, kept in the repository (see README.md)
EVIDENCE := $(D)/label_readings.csv $(D)/extra_labels.csv $(D)/date_check_summary.csv $(D)/xe_headers.parquet $(D)/xe_plate_centres.csv
SAMPLE   := $(A)/sample.py $(DATES) $(LOG) $(EXT) $(ML) $(D)/external/poss_i_observing_dates_strict.csv

.PHONY: all help catalogue dates shadow nuclear alignments null-checks shadow-checks ml-checks ml-retrain site page
all: page

help:
	@echo "make            everything (the page, docs/index.html)"
	@echo "make catalogue  each transient on its plate              -> $(PLATES)"
	@echo "make dates      corrected plate dates                    -> $(DATES)"
	@echo "make shadow     Earth's shadow: pooled count and fit     -> $(POOLED) $(FIT)"
	@echo "make nuclear    nights near US nuclear tests             -> $(NUCLEAR)"
	@echo "make alignments track search over 1,000 shifted copies   -> $(ALIGN)  (slow)"
	@echo "make null-checks the alignment null made other ways      -> $(NULLCHK)  (slow)"
	@echo "make shadow-checks the pooled shadow null: calibration, leak, injection -> $(SHCHK)"
	@echo "make ml-checks  VASCO-ML's ratings: labels, plate-level weighting, shuffles -> $(MLCHK)"
	@echo "make ml-retrain VASCO-ML's model retrained on other labels -> $(MLRT)  (slow; PYML=python with analysis/requirements-vasco-ml.txt)"
	@echo "make site       charts and numbers for the page          -> $(SITE)"
	@echo "make page       the page                                 -> $(PAGE)"

catalogue: $(PLATES)
dates: $(DATES)
shadow: $(POOLED) $(FIT) $(SOUTH)
nuclear: $(NUCLEAR)
alignments: $(ALIGN) $(AFIT) $(VCAND)
shadow-checks: $(SHCHK)
ml-checks: $(MLCHK)
ml-retrain:
	$(PYML) $(A)/ml_retrain.py
	$(PY) $(A)/ml_retrain_eval.py
null-checks:
	$(PY) $(A)/alignment_null_checks.py
site: $(SITE)
page: $(PAGE)

# downloaded inputs: fetched once, when missing, from fixed commits (checksums in fetch_external.py)
$(EXT) $(TESTSV) $(ML) $(LABELS):
	$(PY) $(A)/fetch_external.py

$(LOG):
	$(PY) $(A)/fetch_palomar_log.py

$(PLATES): $(A)/plates.py $(A)/geometry.py $(EXT)
	$(PY) $<

$(DATES): $(A)/date_corrections.py $(A)/date_check.py $(PLATES) $(LOG) $(EVIDENCE)
	$(PY) $<

$(DUST): $(A)/plate_dust.py $(SAMPLE)
	$(PY) $<

$(POOLED): $(A)/shadow_pooled.py $(A)/shadow.py $(A)/geometry.py $(SAMPLE) $(PROPS)
	$(PY) $<

$(SOUTH): $(A)/shadow_plates_log.py $(A)/shadow.py $(DATES) $(LOG)
	$(PY) $<

$(FIT): $(A)/shadow_fit.py $(A)/shadow.py $(POOLED) $(SAMPLE) $(LABELS)
	$(PY) $<

$(DEPTH): $(A)/plate_depth.py $(A)/stars.py $(A)/geometry.py $(SAMPLE)
	$(PY) $<

$(NUCLEAR): $(A)/nuclear.py $(SAMPLE) $(PROPS) $(D)/us_shot_times.csv $(TESTSV)
	$(PY) $<

$(ALIGN): $(A)/alignments.py $(PLATES)
	$(PY) $<

$(VCAND): $(A)/vasco_candidates.py $(D)/vasco_2025_candidates.csv $(EXT) $(D)/xe_headers.parquet $(D)/xe_plate_centres.csv
	$(PY) $<

$(AFIT): $(A)/alignment_fit.py $(A)/alignments.py $(ALIGN) | $(NULLCHK)
	$(PY) $<

# slow (about 40 minutes on 8 cores); its output is kept in the repository, so `make` uses it as it is (no
# prerequisites: built only when missing) and only `make null-checks` reruns it
$(NULLCHK):
	$(PY) $(A)/alignment_null_checks.py

$(MLCHK): $(A)/ml_checks.py $(A)/shadow_fit.py $(A)/shadow.py $(SAMPLE) $(PROPS) $(LABELS) $(POOLED) $(D)/us_shot_times.csv $(TESTSV)
	$(PY) $<

# slow, and needs VASCO-ML's own versions; its output is kept in the repository, so `make` uses it as it is (built only
# when missing) and only `make ml-retrain` reruns it
$(MLRT):
	$(PYML) $(A)/ml_retrain.py
	$(PY) $(A)/ml_retrain_eval.py

$(SHCHK): $(A)/shadow_null_checks.py $(A)/shadow_pooled.py $(A)/shadow.py $(SAMPLE) $(PROPS)
	$(PY) $<

$(SITE): $(A)/build_site.py $(A)/shadow.py $(A)/asteroid_clock.py $(A)/fetch_labels.py $(SAMPLE) $(EVIDENCE) $(DUST) $(PROPS) $(POOLED) $(SOUTH) $(FIT) $(NUCLEAR) $(ALIGN) $(AFIT) $(VCAND) $(MLCHK) | $(NULLCHK) $(MLRT)
	$(PY) $<

$(PAGE): site/assemble.py site/copy.md site/template.html $(SITE)
	$(PY) $<
