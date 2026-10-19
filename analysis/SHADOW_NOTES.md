# Earth's shadow: method notes

VASCO's catalogue only; geocentric shadow (8.5 deg about the anti-sun point at mid-exposure); corrected dates; nights in
VASCO's list; transients in the central 6 x 6 deg of each plate (|xi|, |eta| < 3 deg). No transients removed. In the plate-level fit, each plate's shadowed share is weighted by the average density of
transients over the face of a plate (other plates counted equally), so the excess at edges and corners is allowed for.

- Pooled shortfall (shadow_test.py, shadow_paste_null.py): 330 in shadow against 481 from an average plate; the nine
  footprints pasted on random plates give a count this low 25% of the time (0.7 sigma), not 7.6 sigma.
- Plate-level fit (shadow_fit.py): O ~ Bin(n, r f / (r f + 1 - f)) per plate, one common r over the eight plates partly
  in shadow, f the density-weighted shadowed share. r = 1.02; counting-only 68% 0.94-1.11.
- Calibration on real plates: each real footprint laid on other VASCO plates with a similar number of transients
  (no-effect centre 1.06, 68% 0.84-1.32); effects injected by random removal; Neyman interval 0.78-1.20 (68%),
  0.68-1.48 (90%). Cross-checks in data/shadow_fit.json: other host choices, the shadow moved over its own plate,
  each plate left out, the fit on its own binomial model.

Why not the earlier routes: shadow-sized discs on the same plate leak a real effect into any correction built from them
(closure: injected x0.68 came back x0.70), and their per-plate factors follow each plate's own clump layout. On the 6 x 6 deg square, weighting
by the density map moved r from 1.07 (area shares) to 1.02 and the no-effect centre from 1.10 to 1.06. Removing clumped transients is not needed: the trials keep the
clumps in the data and in the null alike.
