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
