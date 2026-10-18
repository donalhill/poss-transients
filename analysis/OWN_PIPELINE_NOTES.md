# Our own transient search: first plate (XE320), 18 October 2026

Goal: check whether an independent pipeline recovers VASCO's transients from the DSS scans, before using it on the
southern plates crossed by Earth's shadow that are missing from VASCO's catalogue (`data/shadow_plates_log.csv`).

## Steps (scripts in `analysis/`)
1. `fetch_scan.py` - full plate scan from STScI (392 MB), parallel byte ranges, resumable.
2. `gaia_field.py`, `catalog_field.py` - Gaia DR3 (36 boxes), AllWISE (IRSA TAP), Pan-STARRS DR2 (MAST API cones,
   objects with >= 3 detections; single detections are mostly spurious and would veto 77% of random positions).
3. `own_pipeline.py` - photutils detection in 2000-px tiles (3 sigma on the smoothed image, >= 3 pixels, deblended),
   GSC-1 plate solution re-calibrated against Gaia at the plate epoch (median offset 0.4-0.5").
4. `vasco_stamps.py` - significance on the scan at VASCO's positions, against faint Gaia stars and random positions.
5. `compare_vasco.py` - Solano et al. (2022) cuts as far as they can be reproduced, and recovery of VASCO's list.

## Results on XE320 (152 VASCO transients)
- At VASCO's positions the scan shows a peak above 5 sigma for all 152 (median 30 sigma in a fixed-aperture stamp
  measurement, `vasco_stamps.py`; median 52 sigma as the detected segment's peak over the smoothed-image noise,
  `compare_vasco.py`); random positions 3%. The transients are bright, compact spots, about 2x2 to 3x3 pixels.
- Our detection finds all 152 (median offset 0.5"). With Solano-like cuts (FWHM 2-7 px, elongation < 1.3, no Gaia or
  Pan-STARRS within 5", bright-star rule, plate edge band beyond 3.1 deg excluded, as VASCO has no transients there)
  plus a total-light threshold chosen at VASCO's 2nd percentile (tuned on VASCO's list, not a Solano cut), we keep
  134 of 137 and find 522 further candidates: at the brightness of VASCO's faintest, about 4x as many.
- The earlier 14,342 was before the Pan-STARRS, bright-star and brightness cuts.
- The extras are slightly sharper (0.79 of a star's width at the same peak, VASCO's 0.87) and fainter in total light.
  The remaining gap is likely SExtractor's SPREAD_MODEL and FLAGS, the bounding-box rule, and the undocumented step
  from Solano's 298,165 candidates to the 107,875 catalogue (about 318 per plate down to 170 per plate).

## Lessons
- Decide whether sources are bright or dark once per plate: per tile it fails at the plate edges (unexposed border).
- Threshold on the noise of the smoothed image, not the raw one.
- Downloads from US archives were capped at ~0.1-0.5 MB/s on the home connection; over a US VPN 3-4 MB/s.
  Cloud copies of Gaia, Pan-STARRS (HATS) and AllWISE (HEALPix Parquet) exist on AWS for larger runs.

## Closing the gap with all of Solano's cuts (`solano_cuts.py`, `solano_compare.py`, `eso_check.py`)
SExtractor-style re-measurement: 5-sigma detection on a 3x3-filtered image, SNR_WIN, blended flag, half-maximum FWHM,
FWHM and elongation within 2 sigma (MAD x 1.4826) of each 2000-px tile's median, bounding box, SPREAD_MODEL against a
PSF stacked from PSFEx-like stars. Inside 2.9 deg of the plate centre (VASCO's transients on XE320 stop there; ours
pile up at 2.9-3.1 deg, the scan's ragged edge). Output in `logs/solano_stages_XE320.txt`:

| stage | VASCO | others |
|---|---|---|
| inside 2.9 deg | 146 | 42,876 |
| no Gaia/Pan-STARRS within 5", not near a bright star | 144 | 1,851 |
| SNR_WIN > 30, FLAGS = 0 | 144 | 1,027 |
| shape cuts (FWHM, elongation, box) | 113 | 146 |
| S/N > 30 in the ESO compressed cutout VASCO searched | 113 | 126 |
| + SPREAD_MODEL > -0.002, PSF from faint stars | 77 | 19 |
| + SPREAD_MODEL, relative to stars of the same S/N | 55 | 12 |

- The 4x came from cuts not yet applied (FWHM/elongation relative to the tile median, SNR_WIN, blending) and from the
  plate edge. Without SPREAD_MODEL we still have about 2x (113 VASCO + 126 others). The others are much sharper than
  VASCO's transients (median 2.26 px vs 2.99 px; 60% at the 2x2-pixel minimum vs 13%; sharpness relative to stars
  of the same S/N, median -0.0067 vs -0.0021), which is what SPREAD_MODEL > -0.002 removes. Any SPREAD_MODEL
  threshold that keeps 95% of VASCO's transients leaves 45 others (152 in total against VASCO's 146 inside 2.9 deg).
- SPREAD_MODEL can't be pinned down on photographic plates: faint stars are sharper than bright ones (median -0.005 at
  S/N 30-50, +0.02 at S/N 200-500 against one PSF), so the result depends on which stars PSFEx used. Our
  reproduction also drops 30-90 of VASCO's own transients, so it is not exact: the remaining others are within that
  mismatch.
- ESO's DSS-1 cutouts (what VASCO searched) are the compressed DSS, differing from the STScI full-resolution scan by
  about a third of the pixel noise. 95% of the remaining others still exceed S/N 30 there, so compression is not the
  explanation.
- The 107,875 catalogue is Solano's intermediate 298,165 list, re-filtered in 2025 for Gaia/Pan-STARRS within 5",
  duplicates, and southern plates; no WISE/VOSA, SuperCOSMOS, or visual check (Villarroel et al. 2026, arXiv:2602.15171).
  Live at http://svocats.cab.inta-csic.es/vanish-bruehl26/ with a `prob` column. So there is no undocumented step.

## Which stars PSFEx used (`psfex_reverse.py`, `logs/psfex_reverse_XE320.txt`)
- PSFEx's default selection (S/N > 20, ellipticity < 0.3, FWHM within 20% of the modal FWHM) settles on stars of
  FWHM ~2.6 px, i.e. the faint S/N 20-100 stars; bright stars (FWHM 4-8 px on these plates) fall outside the 20% band.
  A PSFEx-weighted PSF and plain PSFs from S/N 20-50 stars give the same answer; bright-star PSFs are ruled out.
- Our SPREAD_MODEL had a -0.003 bias (stack blurring); corrected so each PSF's own stars average zero (no VASCO input).
- Our measurement scatter at S/N 30-70 is 0.0018 (12% of genuine faint stars fall below -0.002), comparable to the
  cut, so VASCO's hard edge at -0.002 shows up blurred: with the PSFEx-like PSF their median is -0.0014, 5th pct -0.0052,
  as expected for a distribution cut at -0.002 and smeared by ~0.002.
- After all cuts: VASCO 63 of 113 pass in our measure, others 13 of 126 (PSFEx-like PSF); 81 and 15 with the S/N 20-50
  PSF. Others are removed five times as often as VASCO's transients; scaled to VASCO's count, about 0.2x extra.

## Frozen selection and the second plate (18 October)
- `frozen_selection.py`: all of the above, with the FWHM/elongation cuts at 3 sigma (our measures are noisier), local
  background noise per 64-px mesh (as SExtractor), and two thresholds fixed on XE320 from where VASCO's transients stop:
  SNR_WIN > 36.6 on our scale (their 30; our S/N reads ~15% high) and SPREAD_MODEL > -0.0051 (bias-corrected).
- XE320: keeps 81% of VASCO's transients, total 1.06x their count. XE448 (b = 23 deg, ten times the stars; not used for
  tuning): keeps 61%, total 1.03x. The objects that differ sit at the faint S/N edge.

## Southern shadow plates
- The Palomar log's plate numbers are not unique (1149 appears twice on 1954-07-28, at +72 and -17 deg) and are the
  blue-plate numbers; the DSS red scan is usually labelled one higher (E1150). Scans are matched to log entries by
  position (< 1 deg) and date (+-2 days), from every scan header (`data/xe_scan_centres.csv`). All 22 have scans
  (`data/southern_shadow_scans.csv`); 5 are 15-micron scans, resampled to 25 microns (`rebin_scan.py`).
- Crowding (`crowding.py`, `data/crowding.csv`): share of random positions with Gaia or Pan-STARRS within 5". XE320 11%,
  XE448 32%. Plates above 50% are left out: 1105 (XE810), 1071 (XE807), 172 (XE868), 158 (XE809; Pan-STARRS cones
  capped at 500,000 sources).
- ESA's Gaia archive returned errors from ~03:00 on 18 October; Gaia DR3 now comes from IRSA (`gaia_field.py`).
