# POSS-I transients re-examined

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23549537.svg)](https://doi.org/10.5281/zenodo.23549537)

A re-examination of three statistical tests reported by the VASCO project on transients in the red plates of the first
Palomar Observatory Sky Survey (1949 to 1957): fewer transients inside Earth's shadow, more transients near the dates of
US nuclear tests, and transients lying along straight lines. The tests are repeated on VASCO's public catalogue, with
every plate's date checked and corrected where required and its time calibrated against the asteroid trails it
recorded.

**The page: https://donalhill.github.io/poss-transients/** (built into `docs/index.html`). Every chart and every number
in its text comes from the chain below.

## Running it

```
python -m venv .venv && .venv/bin/pip install .     # the pinned versions the page was built with
make tests report PY=.venv/bin/python
```

`make -n` shows what would run and why. A run from scratch downloads VASCO's public data (`vasco/fetch_external.py`,
from fixed commits, checked against the files analysed), Palomar's log, and the SFD dust map. On an 8-core laptop it
takes a few hours; the slowest steps are the replicas behind the shadow test and the alignment search's 1,000 random
catalogues.

## The chain

1. **Plates and dates** (`vasco/`, `dates/`). Each transient of VASCO's catalogue is assigned to its plate
   (`vasco/plates.py`). Each plate's date is checked against its handwritten label, the asteroid trails it recorded,
   and Palomar's log, and corrected where required (`dates/date_corrections.py`). Each plate's mid-exposure time,
   calibrated against its own asteroid trails or the survey's measured clock offset, is in `data/plate_times_all.csv`,
   with its posterior in `data/plate_time_posterior.csv` and draws from it in `data/plate_time_draws.parquet`.
2. **Samples** (`build_sample.py`). The catalogue on each plate's corrected date and time, as the tests read it
   (`analysis/sample.py`): each transient counted by VASCO-ML's rating (`release/vasco_rated`, the tests' sample) or
   once (`release/vasco_raw`), with the average density of transients over the face of the plates the shadow never
   reaches as the acceptance (or a flat one, `release/vasco_flat`, as a check).
3. **Tests** (`make tests`).
   - Earth's shadow: the transients observed in shadow against those expected from each plate's own total
     (`analysis/shadow_fit.py`, `analysis/shadow_stats.py`), its null and the limit on glints from replicas of the
     real shadows on other plates (`analysis/replicas.py`), the cross-check with shuffled plate times
     (`analysis/permutation_null.py`), the plate times' systematic (`diagnostics/shadow_total.py`), and VASCO's count
     repeated step by step (`diagnostics/vasco_repeat.py`).
   - Nuclear-test dates: `analysis/nuclear.py`, with each plate's expected count from where it pointed
     (`analysis/expected.py`, the pointing model).
   - Alignments: the track search against random catalogues built from the data (`analysis/alignments.py`), its checks,
     and the fit (`analysis/alignment_*.py`).
   - VASCO-ML retrained eight ways (`diagnostics/vasco_ratings.py`), from ratings in `data/external/vasco_ml/ratings/`
     made with VASCO's own training code.
4. **Page** (`make report`). `report/registry.py` reads every number the page quotes from the tests' outputs and
   `report/figures.py` every chart's data; the text, by section in `report/copy/` with a brief per slot in
   `report/briefs/`, quotes numbers as `{n:key}`, and `report/check.py` holds it to them, so the text cannot quote a
   number the chain did not produce. The page is assembled into `docs/`, with the hashes of its inputs in
   `docs/PROVENANCE.json`. `python report/edit.py` serves the page locally with its text editable in place.

## Data

The small input tables are in `data/`: the plates' date records and corrections, their calibrated times, the asteroid
matches behind the calibration, the plates' depth and dust, Palomar's log (VizieR VI/25), and the US test times
(DOE/NV-209). The plate times (`data/plate_times_all.csv`, with their posterior, draws, and the asteroid matches behind
them) were measured from the asteroid trails on each plate's scan, as the page describes, and `make plate-times`
remakes them: `dates/asteroid_list.py` lists the numbered asteroids on each plate from the MPC's orbits (it needs
ASSIST and REBOUND, most easily from conda-forge), `dates/asteroid_clock.py` takes each asteroid's track from JPL
Horizons and measures its trail on the plate's own scan, and `dates/plate_times.py calibrate` turns the trails into
each plate's time. The label readings (`data/label_readings.csv`) were read from crops of each plate's label
(`dates/fetch_labels.py`) and checked by eye; the photographs on the page are cut from the DSS scans, and
`report/plate_image.py` makes the whole-plate image from plate XE344's scan. VASCO's catalogue and its companion files
are fetched, not kept here; the retrained VASCO-ML ratings are remade by `make` with `vasco/ml_retrain.py`.

## Licence

Code under the MIT License (`LICENSE`); the write-up, figures, and data tables under CC BY 4.0
(`LICENSE-CONTENT.md`).
