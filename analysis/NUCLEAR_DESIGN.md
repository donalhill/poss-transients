# Nuclear tests: the analysis design

Written and committed before the final run of `nuclear_fit.py`. The commit holding this file and the code is quoted on
the page. Everything below was fixed before that run; the list at the end says which choices were made after seeing
earlier results, and which tests had never been run on the real window.

## The claim under test
VASCO (Bruehl & Villarroel 2025, Sci Rep 15:34125; Bruehl et al. 2026, Sci Rep 16:31890): observing nights within
one day of a US nuclear test hold more transients than other nights. The window (the test date, the day before and the
day after) was fixed by the authors before they looked at the data. They read the association's sharpness in time
("temporal specificity") as evidence that it is real rather than a quirk of the observing schedule.

## Unit, outcome, exposure
- Unit: the photograph. Transients are found plate by plate; a night's count is the sum over however many plates were
  taken, which VASCO's test does not allow for.
- Outcome: the plate's transients in VASCO's catalogue. Clump-robust version, reported alongside: the number of
  2' x 2' cells on the plate holding at least one transient, so a dense patch of plate damage counts once per cell.
- Exposure (primary, VASCO's): the plate's corrected UT date is a US test's date, the day before or the day after.
  VASCO's test dates are the shots' UT dates (97 of 99 agree with DOE/NV-209).
- Exposure (secondary, physical): the plate's mid-exposure lies 0-24 h after a shot; 24-48 h after; 0-24 h before.

## Design
- Plates are compared only with plates from the same observing run (consecutive dark nights; a new run after more
  than 4 nights with no plates). Each run has its own free level, so anything that drifts across the survey
  (emulsion batches, processing, season, prevailing weather) drops out with no assumed shape. This is the
  time-stratified design of environmental epidemiology (case-crossover).
- Within a run, what differs from plate to plate enters the model: distance from the Milky Way (a smooth curve),
  dust, exposure time, airmass, moonlight, distance from the Moon, VASCO's plate-quality flag, ecliptic latitude,
  completeness (the share of Gaia stars with 18 < G < 20 found on the plate) and the width of Gaia-matched stars.
  Completeness and width are measured on Gaia stars only, so a plate's own specks and transients cannot change them.
- Negative binomial counts. The estimate: r, the ratio of transients on exposed plates to those on unexposed plates
  of the same run, after the plate properties. r = 1 means no effect. Runs holding no unexposed plates, or no
  exposed ones, carry no information about r; the identification table lists every run with an exposed plate.

## Inference
- Two ways of sampling the null (the test schedule unrelated to transient counts), fitted exactly like the data,
  with the scatter between plates re-estimated in each trial:
  (a) within each run, the run's exposed nights reassigned at random among its nights;
  (b) the real test calendar slid by 1 to 5 whole years and up to 20 days either way (keeping seasons and the
      clustering of tests), leaving out slides whose window falls largely on the real one.
  p = (k + 1) / (n + 1). If the two disagree by more than a factor of two, that is reported as a finding.
- Interval on r: Neyman, from effects injected into the trials of (a). Power: how often an excess as large as
  VASCO's (1.69) would give a value as low as the data's.

## Hypotheses and look-elsewhere
- Primary: VASCO's window. Local p, no correction: the window is theirs and was fixed in advance.
- Secondary: the physical windows above.
- Look-elsewhere, by a fixed rule: every window of consecutive days [a, b] within 3, 7 and 14 days of a test, with
  at least 5 plates; the global p at each size is the share of trials whose best window does at least as well as
  the data's best. And the rank of VASCO's window among them, against how often chance ranks it that high.

## Robustness (reported as they come out)
- Each run, and each test series, left out in turn.
- Run definition: a new run after more than 3, 4, 6 or 10 empty nights.
- Raw counts against the clump-robust count; the central 5 x 5 deg of each plate only.

## Mechanism (run last, only on the frozen analysis)
- The blue plates of the same fields, taken the same nights: detections with no Gaia star within 5", exposed against
  unexposed nights within runs. Anything in the emulsion, or settling on the plates, would add such detections to
  blue plates too; VASCO's red-plate transients would not appear there.
- Nevada against Pacific tests.

## Decided after seeing earlier results, or not yet run on the real window
- Decided after seeing results: comparing within observing runs (adopted after a smooth time trend proved to move
  the result: r from 1.14 to 1.48 with its flexibility); adding plate depth and image width (adopted after the
  excess survived the other plate properties; the first depth measure, which counted every detection, was replaced
  by completeness because specks on a plate inflated it); the physical windows in hours (after reading VASCO's
  remark that plates are exposed the night before morning tests).
- Never run on the real window before the frozen run: the clump-robust outcome, the within-run night reassignment
  null, leave-one-run-out, leave-one-series-out, the 5 x 5 deg cut within runs, and the blue-plate test.
- Seen before freezing: within runs, with the earlier depth measure, r was about 1.7 at a local p of about 0.006,
  and a global p of 0.04, 0.08 and 0.26 for windows within 3, 7 and 14 days.

## Addendum, specified after the freeze and before its result existed (19 October 2026)
- Prediction out of sample: each observing night held out in turn; its plates predicted from the within-run model
  fitted to every other night (the run's level from its other nights). Held-out negative binomial log-likelihood
  with and without the test-date feature (VASCO's window as one feature; and the day before, the day and the day
  after as three). The gain is calibrated against the same feature assigned to random nights within runs.
- The blue-plate measurement (plate_blue.py) was taken while the frozen run was going; its comparison of exposed and
  unexposed nights is run only after the frozen results.
- Fix, after the first frozen run stopped on a numerical failure and before any of its results were read: a trial
  whose window holds no plate sharing a run with an unexposed plate leaves r undetermined; such trials are dropped
  and the number kept is reported. Nothing else changed.
- Speed, before any results were read: the data are fitted once in full; every trial keeps the plate-property slopes
  of that fit and refits the run levels, r and the scatter (a negative binomial fit vectorised over trials). On
  random windows it agrees with full refits to within 1-2.5%. Trials raised to 2,000 per null and 500 per injected
  effect, which this makes cheap.
- Tracks, specified after the blue-plate and prediction results and before this one: the outcome is each plate's
  number of straight track candidates (alignments.py's track score at 1 in 10 or better, and separately 1 in 100 or
  better; each transient in at most one track, best first), with the plate's own number of transients as a further
  covariate, so the question is whether transients near tests line up more often than their numbers imply. Same
  within-run model and null (a).
- Diagnostic, found after the frozen results (from the map of every window): the hot cells at days -14 and -6 and
  VASCO's window are the same two nights, 5 November 1951 and 23 March 1955, each at day 0 from one shot and 14 days
  before another of the same series. Reported: each night's transients over expected; how many nights in the survey
  are as rich, and how many of them fall near tests (with the hypergeometric chance); r without the two nights.
- Diagnostic, found after the frozen results (from the chart of nights near tests): the median, over the nights near
  tests, of each night's transients over expected, with the run levels fitted with the window as a factor in the data
  and in every trial alike; also with each set's two richest nights left out. The statistic was chosen after seeing
  the chart.
- Correction: the "typical night" statistic above borrowed the fitted overall ratio through each run's level (a run
  with no other plates showed exactly the fitted ratio), so it was not independent evidence and is withdrawn. In its
  place, each informative run on its own: plates near tests against the run's other plates, each against its expected
  count from the plate properties alone (no run levels, no window), with nothing shared between runs.
- Added after the ground-up results and before these were run (nuclear_basic.py): (i) a within-plate ABCD test on
  object classes, edge (more than 2.5 deg from the centre in either axis) against central, and crowded (another
  transient within 5') against isolated; the ratio of class shares, plates in the window against the other plates of
  their runs (Mantel-Haenszel over runs), calibrated by the same within-run reassignment; (ii) a 95% upper limit on the
  Mantel-Haenszel ratio, from effects injected into the within-run reassignments.
- Shapes, specified before any were measured (plate_shapes.py, nuclear_basic.py): on every plate of the eight runs that
  hold plates near tests and other plates, up to 60 transients picked at random (seeded) get a 4' DSS cutout; each is
  fitted with an elliptical Gaussian, as are the stars on the same cutout, and its width is compared with the median
  width of stars of similar peak brightness (0.5 to 2 times). Star-like: within the 10-90% range of the same ratio for
  the stars themselves. The share of star-like transients, plates near tests against the other plates of their runs
  (Mantel-Haenszel odds over runs), calibrated by within-run reassignment. Glints imaged by the telescope predict more
  star-like transients near tests; flaws in the emulsion fewer.
