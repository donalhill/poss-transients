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
NUCFIT   := $(D)/nuclear_fit.json
NUCCHECK := $(D)/nuclear_checks.json
BLUE     := $(D)/plate_blue.csv
PROPS    := $(A)/plate_properties.py $(DUST) $(DEPTH) $(D)/xe_headers.parquet
ALIGN    := $(D)/alignments.json
VCAND    := $(D)/vasco_candidates.json
SHAPES   := $(D)/tail_shapes.json
SECOND   := $(D)/second_copy.json
TRAIL    := $(D)/xe492_trail.json
SITE     := docs/data/site.json
PAGE     := docs/index.html

# inputs read off the plates by eye or from slow public services, kept in the repository (see README.md)
EVIDENCE := $(D)/label_readings.csv $(D)/extra_labels.csv $(D)/date_check_summary.csv $(D)/xe_headers.parquet $(D)/xe_plate_centres.csv
SAMPLE   := $(A)/sample.py $(DATES) $(LOG) $(EXT)

.PHONY: all help catalogue dates shadow nuclear alignments tracks site page
all: page

help:
	@echo "make            everything (the page, docs/index.html)"
	@echo "make catalogue  each transient on its plate              -> $(PLATES)"
	@echo "make dates      corrected plate dates                    -> $(DATES)"
	@echo "make shadow     Earth's shadow: pooled count and fit     -> $(POOLED) $(FIT)"
	@echo "make nuclear    nights near US nuclear tests             -> $(NUCLEAR) $(NUCFIT) $(NUCCHECK)"
	@echo "make alignments track search over 1,000 shifted copies   -> $(ALIGN)  (slow)"
	@echo "make tracks     the strongest tracks on the plates       -> $(SHAPES) $(SECOND) $(TRAIL)"
	@echo "make site       charts and numbers for the page          -> $(SITE)"
	@echo "make page       the page                                 -> $(PAGE)"

catalogue: $(PLATES)
dates: $(DATES)
shadow: $(POOLED) $(FIT) $(SOUTH)
nuclear: $(NUCLEAR) $(NUCFIT) $(NUCCHECK)
alignments: $(ALIGN) $(VCAND)
tracks: $(SHAPES) $(SECOND) $(TRAIL)
site: $(SITE)
page: $(PAGE)

# downloaded inputs: fetched once, when missing
$(EXT):
	$(PY) $(A)/fetch_external.py

$(LOG):
	$(PY) $(A)/fetch_palomar_log.py

$(PLATES): $(A)/plates.py $(A)/geometry.py $(EXT)
	$(PY) $<

$(DATES): $(A)/date_corrections.py $(PLATES) $(LOG) $(EVIDENCE)
	$(PY) $<

$(DUST): $(A)/plate_dust.py $(SAMPLE)
	$(PY) $<

$(POOLED): $(A)/shadow_pooled.py $(A)/shadow.py $(A)/geometry.py $(SAMPLE)
	$(PY) $<

$(SOUTH): $(A)/shadow_plates_log.py $(A)/shadow.py $(DATES) $(LOG)
	$(PY) $<

$(FIT): $(A)/shadow_fit.py $(A)/shadow.py $(POOLED) $(SAMPLE)
	$(PY) $<

$(DEPTH): $(A)/plate_depth.py $(A)/stars.py $(A)/geometry.py $(SAMPLE)
	$(PY) $<

$(BLUE): $(A)/plate_blue.py $(A)/geometry.py $(SAMPLE)
	$(PY) $<

$(NUCLEAR): $(A)/nuclear.py $(SAMPLE)
	$(PY) $<

$(NUCFIT): $(A)/nuclear_fit.py $(A)/nuclear_model.py $(SAMPLE) $(PROPS) $(D)/us_shot_times.csv
	$(PY) $<

$(NUCCHECK): $(A)/nuclear_checks.py $(A)/nuclear_model.py $(A)/alignments.py $(SAMPLE) $(PROPS) $(BLUE) $(D)/us_shot_times.csv
	$(PY) $<

$(ALIGN): $(A)/alignments.py $(PLATES)
	$(PY) $<

$(VCAND): $(A)/vasco_candidates.py $(D)/vasco_2025_candidates.csv $(EXT) $(D)/xe_headers.parquet $(D)/xe_plate_centres.csv
	$(PY) $<

$(SHAPES): $(A)/tail_shapes.py $(A)/stars.py $(A)/geometry.py $(ALIGN) $(DATES)
	$(PY) $<

$(SECOND): $(A)/second_copy.py $(A)/geometry.py $(ALIGN) $(DATES)
	$(PY) $<

$(TRAIL): $(A)/xe492_trail.py $(A)/geometry.py $(ALIGN) $(DATES)
	$(PY) $<

$(SITE): $(A)/build_site.py $(A)/shadow.py $(A)/asteroid_clock.py $(A)/fetch_labels.py $(SAMPLE) $(EVIDENCE) $(DUST) $(POOLED) $(SOUTH) $(FIT) $(NUCLEAR) $(NUCFIT) $(NUCCHECK) $(ALIGN) $(VCAND) $(SHAPES) $(SECOND) $(TRAIL)
	$(PY) $<

$(PAGE): site/assemble.py site/copy.md site/template.html $(SITE)
	$(PY) $<
