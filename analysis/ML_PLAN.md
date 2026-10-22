# Using VASCO-ML's ratings: analysis plan

Agreed and carried out: the checks are in `ml_checks.py`, `ml_retrain.py`, `ml_retrain_eval.py`, `shadow_fit.py`
(bands, joint fit, limits) and `nuclear.py` (bands, limit, the rated block); section 9 records what they found.

The aim: VASCO's 2026 analyses count each transient by a machine-learning rating. The page must show that we use the
ratings, understand what they encode, and demonstrate whether they change the shadow and nuclear results, with every
claim backed by a test that could have come out the other way.

## 0. What has already been looked at

Honesty about this goes on the page. Already run on the real data, in this order:

| Test | Defined before seeing its result? |
|---|---|
| Shadow and nuclear, every transient weighted by its rating (pooled and plate by plate) | Yes: VASCO's own weighting |
| Shadow, plate-by-plate fit in each tenth of the ratings (VASCO's decile edges) | Yes: VASCO's own split |
| Shadow, top tenth singled out (pooled and plate by plate) | No: chosen because VASCO highlight it; trials factor computed (19%) |
| Shadow, ten bands fitted at once with a free purity each | Defined after seeing the per-band ratios |
| Shadow, trend of purity with rating | Defined after seeing the per-band ratios |

Everything new in this plan is fixed here before it is run.

## 1. What the rating is (facts, with sources)

From VASCO-ML at commit 418e8c8 (the commit `fetch_external.py` pins):

- **Model:** an ensemble of four tree models (RF, GBM, XGBoost, LightGBM), mean probability (`MLscripts/pipeline_train.py`,
  `config_red_allcandidates.yaml`). Cross-validated AUC 0.81, sensitivity and specificity about 0.71 (README).
- **Not calibrated:** isotonic calibration "collapsed the probability distribution"; the raw ensemble mean is used
  (README). So a rating orders transients; it is not a probability of being real.
- **Labels:** visual judgements of real point source (1) or plate flaw (0). `data/combined_training_labels.csv` holds
  393 labels (141 real, 252 flaw) from four files: `training_data_r1` (79/32), `beatriz_r4_labels_conservative`
  (51/47), `bruehl_r2_labels` (4/37), `steve_batch_xlsx` (7/136). The README and STATISTICAL_METHODS.md say 250 labels
  (134 real, 116 flaw). To ask.
- **Label selection used nuclear dates:** `pipeline_active_learning.py` selects candidates for labelling from three
  buckets, the third "Nuclear-date hard negatives (low-prob on test dates)" (25% of each batch, `frac_nuclear_negative`).
  The README states the shadow and nuclear signals were "never used during training".
- **Plate-level inputs are per night:** `common.compute_plate_features(df, group_col='obs_date')` computes
  `plate_low_snr_frac`, `plate_snr_std`, `plate_n_high_snr`, and `plate_elongation_mean` over all transients on a
  night, across all of that night's plates. `plate_n_candidates` in the scored catalogue equals the night's transient
  count exactly.
- **What drives each rating** (each transient's top SHAP feature, from the scored catalogue): red gradient magnitude
  29%, plate quality flag 27% (GOOD 20%, MODERATE 7%), aperture flux 16%, night low-SNR fraction 10%, distance to scan
  edge 7%, night mean elongation 6%. A plate- or night-level input is the top driver for about 43% of transients.
- **Rating against plate and night:** about a third of rating variance lies between plates. A night's transient count
  and its mean rating have Spearman 0.54. Rating against |Galactic latitude| 0.26, distance to scan edge 0.24.
- **Run to run:** the released `scored_catalog.zip` is not the run behind VASCO-ML's decile table (top-decile edge 0.751
  in the released file, 0.781 in theirs). `reproduce.py` retrains and checks results within tolerances.

## 2. Circularity and leakage: tests, with what would prove each

Each test states its pass/fail criterion now.

**2a. Nuclear-selected labels.** Find the labelled transients on nights within a day of a US test (VASCO's list and the
full list), and their labels. Retrain VASCO's model (if feasible, see 5) with and without them, and rerun the weighted
nuclear test on each set of ratings.
*Proven to matter* if removing them moves the weighted nuclear excess (or its p-value) by more than the spread from
retraining with different seeds.

**2b. The nuclear test sees only plate-level ratings.** Every nuclear statistic (VASCO's per night, ours per photograph
and against pointing) uses each plate's or night's summed rating, which a shuffle of ratings within a plate leaves
unchanged. So, as an identity, weighting the nuclear test by rating is the same as weighting each plate by its mean
rating: the rating's power to tell a real transient from a flaw on the same plate plays no part. What remains to test:
(i) the weighted excess as the unweighted excess times the near-test plates' higher mean rating, step by step;
(ii) whether that higher mean rating is predicted by plate and night properties alone (pointing, the night's transient
count, plate quality, depth), from a model of plate mean rating fitted without the plates near tests, against shifted
calendars;
(iii) how the night-level inputs tie ratings to a night's richness (night count against night mean rating, with and
without pointing).
*Proven* if (ii) predicts the near-test plates' mean rating within the shifted calendars' spread: the ratings then add
nothing to the nuclear test beyond plate and night properties that the pointing allowance and counting per photograph
already address.

**2c. The information test, for the shadow.** For each weighted shadow result (pooled, plate by plate), the same with
ratings shuffled within each plate (all plates, hosts included, so the density maps are rebuilt too).
*Ratings carry information for that test* if the real-minus-shuffled difference lies outside the shuffles' spread. If it
does not, the page says the ratings change nothing in that test beyond plate-level weighting.

**2d. Training examples in the tests.** Which labelled transients lie on the shadow plates or the near-test plates; their
ratings are in-sample (inflated towards their labels). Rerun each test without them.
*Matters* if the result moves by more than its trial spread.

**2e. Sculpting.** On the plates the shadow never crosses, lay the real shadow outlines and compare ratings inside and
outside (per band, and the mean). On nights away from tests, compare ratings on shifted-calendar nights with the rest.
*Sculpting present* if the background's ratings differ by position or by calendar beyond trial spread. The trials already
absorb it; the check shows it explicitly.

## 3. Fixed hierarchy of tests

The same structure for the shadow and the nuclear test.

**Context:** VASCO's construction reproduced (their dates, their counting); counts against pointing.

**Model-independent search (no assumption about how a signal sits among the ratings):**
- Shadow: ten bands fitted at once, a free purity each (signal vanishes in shadow, so r_k = 1 - purity_k); test of no
  signal in any band, null from host plates with all bands on the same hosts.
- Nuclear: ten bands at once, a free excess each near tests, against pointing; null from shifted test calendars with all
  bands shifted together.

**Benchmarks**, each with an observed 95% upper limit on the signal and the expected limit with its ±1σ and ±2σ band
from no-signal trials:
- A. Signal in proportion to rating (the weighted fit).
- B. Signal rising with rating (the trend; a different hypothesis from A and C, reported separately).
- C. Signal confined to the top tenth (with the trials factor, since the band was singled out after the fact).
- D. Signal shaped like VASCO's labelled-real examples, "your labels at face value". Included (see below).

**Primary result:**
- Option 1, the model-independent search: assumes nothing; weaker reach (a VASCO-sized shortfall in every band is found
  38% of the time).
- Option 2, benchmark A (weighted): the most reach if the signal follows the rating (about 52% for a VASCO-sized
  shortfall); assumes the rating means what VASCO say, which section 2 may undermine.
- Recommendation: option 1 as primary, A to D as benchmarks, all reported.
- **Decided:** option 1 primary; A to D as benchmarks; D included, labelled throughout as VASCO's labels taken at face
  value.

## 4. Effect or no effect, operationally

For each test on the page, three numbers side by side: counts, real ratings, and ratings shuffled within plates. The
ratings "have an effect" on a test only if real differs from shuffled beyond the shuffles' spread; "the weighting moves
the result" only if real differs from counts beyond that spread.

## 5. Systematics from the classifier

- **Retraining:** `pipeline_train.py` reads `feature_cache.parquet` and the label files, and needs no plate scans (to
  confirm by running it in isolation, with paths adjusted). If it runs: 20 retrains with different seeds, every benchmark
  rerun on each, and the spread quoted. If not: the released-against-decile-run difference bounds it, stated as such.
- **Binning:** VASCO's decile edges against equal tenths of our own sample; five bands against ten.
- **Radius:** 8.69° against VASCO-ML's 8.50°.
- **Pointing curve:** fitted to counts against summed ratings; terms as in `shadow_null_checks.py`.

## 6. On the page

- A new section after "Number of transients varies greatly across plates": what the ratings are, what drives them
  (facts from section 1), the circularity and leakage results (section 2), and how the tests use them.
- Shadow and nuclear: the hierarchy of section 3, with one expected-against-observed limit chart per test.
- Neutral wording throughout; facts from VASCO-ML's own files, cited by path.
- The questions for VASCO (label counts, label selection, the run behind the deciles, the shadow geometry behind their
  totals) in the email, not on the page.

## 7. Order of work

1. Section 2 checks that need no retraining: 2b, 2c, 2d, 2e, and the label lookup for 2a.
2. Retraining feasibility; then 2a and section 5.
3. Section 3 for the nuclear test (the shadow side mostly exists).
4. Expected-against-observed limits for every benchmark.
5. The page.

## 8. Rule for circular ratings (decided)

Judged test by test. Where the ratings are shown to be circular for a test, that test is based on counts, with the
evidence stated, and the rated result kept as one line ("taken at face value, the ratings give X"). They are not dropped.

## 9. Findings so far (scratch runs; to be moved into the pipeline)

- **Labels** (2a lookup, `ml_checks.py`): 393 labels; 336 match a catalogue transient within 2". 49 lie on nights within
  a day of a test on VASCO's list (14 real, 35 flaw); labelled transients make up 0.4% of transients on those nights and
  0.3% overall. 10 lie on the shadow plates (2 real, 8 flaw).
- **The steve batch's "real" labels are random**: `pipeline_train.load_labels` sets every row of that batch to flaw and
  marks 7 rows chosen by a seeded random draw (`steve_labels: seed 50, n_positive 7`) as real.
- **Which labels trained the released model** is unclear: the documented 250 (134 real, 116 flaw) is exactly the 393
  without the steve batch, but retraining with VASCO's own code on the 250 matches the released ratings less well (rank
  correlation 0.88) than on all 393 (0.95). Neither reproduces them exactly.
- **2b, nuclear**: the ratings raise the near-test plates' mean rating by 8.5% (per photograph 1.27 → 1.37; against
  pointing 1.01 → 1.04). Against a model of plate mean rating from plate and night properties (70% of deviance
  explained, refitted without each calendar's plates), the near-test plates stand 3.5% above prediction; shifted
  calendars centre at 0.95 (a seasonal effect the model misses), so the near-test plates are high at p = 0.077.
  A night's transient count and its mean rating: Spearman 0.54 raw, 0.18 after pointing.
- **2c, shadow information test**: shuffling ratings within plates, the pooled weighted count in shadow is 124 against
  a shuffled median of 126 (p = 0.16), and the plate-by-plate ratio 0.96 against 1.00 (p = 0.17): the ratings' placement
  within plates adds nothing significant to the shadow test.
- **2e, sculpting**: on plates the shadow never crosses, transients inside the real outlines' positions are rated 2%
  lower than outside (0.98; shuffled 0.995 to 1.002): the rating depends on position on the plate face (it rises with
  distance from the scan edge). The page's weighted tests absorb this (density maps and host trials built from rated
  transients); VASCO's binomial tests do not.
- **Retraining** with VASCO's own code (pinned versions) works from the shipped feature cache and labels: 40 retrains
  running (all labels, without the 49 window-night labels, the 250, the 250 without window nights; seeds 1 to 10).
- **2a, retraining results** (medians of 10 seeds each; VASCO's code and pinned versions; scratch):

  | Labels | VASCO nuclear p | pointing ratio | its calendar p | shadow r |
  |---|---|---|---|---|
  | all 393 | 0.016 | 1.050 | 0.105 | 0.942 |
  | 393 without the 49 window-night labels | 0.020 | 1.029 | 0.152 | 0.935 |
  | documented 250 | 0.029 | 1.018 | 0.182 | 0.928 |
  | 250 without window-night labels | 0.021 | 1.034 | 0.168 | 0.909 |
  | released ratings | 0.020 | 1.044 | 0.146 | 0.961 |

  Seed spread within a label set: about ±0.001 in p, ±0.003 in ratios. Removing the window-night labels moves the nuclear
  result well beyond seed spread, but in opposite directions on the two base sets: circularity from label selection is
  **not proven**; dependence on the choice of labels **is** (a systematic: VASCO's p from 0.016 to 0.029 on their own
  construction). Under section 8 the ratings stay as the weighting, with the label-set systematic added.
- **Nuclear by rating band** (`nuclear.py`, released ratings): excess near tests against each band's pointing curve, bands
  1 to 9 between 0.93 and 1.06; top band 1.24, local p 0.02, 0.057 allowing for ten bands; every band at once p 0.38;
  trend p 0.15. 95% upper limit on an excess near tests: 41% (expected with none: 24%, 68% range 10% to 40%).
- **The top band's nuclear excess across label sets** (10 seeds each): all 393 labels 1.26 (released 1.24); without the
  49 window-night labels 1.05; documented 250 1.03; 250 without window nights 1.09.
- **Ablation, per-night inputs held at zero** (3 seeds each): all 393 labels 1.15; without window-night labels 1.13. With
  per-night inputs the 49 window-night labels swing the top band's excess by 0.21; without them, by 0.02. **Demonstrated:**
  labels on test-window nights, through inputs that identify nights, create most of the top band's excess near tests.
  VASCO's per-night headline barely depends on them (p 0.018 without night inputs, 0.016 with).
- **Decision (section 8):** ratings shown circular for the nuclear test, so it returns to counts as its basis, the rated
  result kept as one line with this demonstration as the reason. No circularity shown for the shadow test: it keeps the
  ratings, as VASCO's 2026 work uses them.
- **Decision revised (by the author):** both tests use the ratings by default, as VASCO do; the dependence on the
  label set is shown as a systematic, with the top-band mechanism as its explanation, and counts as a comparison
  (`nuclear.py` section 6).
