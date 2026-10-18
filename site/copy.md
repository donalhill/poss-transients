### hero-title
Re-examining the evidence for glints from orbit in 1950s sky photographs

### hero-lede
Between 1949 and 1957, Palomar Observatory in California photographed the northern sky on glass plates. The plates show faint points of light that have no counterpart in modern sky catalogues. The VASCO project, led by Beatriz Villarroel, has catalogued more than 100,000 of these transients.

Villarroel’s team has proposed that some transients may be produced by sunlight glinting off flat, reflective objects in orbit around Earth in the pre-satellite era. In papers published in 2025 and 2026, the team reported three statistical patterns in support of this theory: fewer transients appear inside Earth’s shadow; more transients appear around the dates of nuclear weapons tests; some transients appear in straight patterns.

The team has kindly made its data catalogue and code public, enabling us to independently repeat all three tests.

### survey-h2
The Palomar Sky Survey

### survey-intro
The plates were taken with Palomar’s 48-inch Schmidt telescope, which was built to photograph wide areas of sky. Each plate is a sheet of glass 14 inches square, and covers a patch of sky 6.5 degrees across (about thirteen full Moons). Every patch was photographed through both a blue and a red filter, where red exposures lasted 40 to 60 minutes. Copies of the plates were scanned in the 1990s, and astronomers now work from those scans.

VASCO compared the red plates with three modern catalogues of the sky: Gaia, Pan-STARRS, and WISE. A point of light became a transient if no modern catalogue lists anything within 5 arcseconds of it. An arcsecond is a 3,600th of a degree, and 60 arcseconds make an arcminute. The public catalogue holds 107,875 such transients spanning 635 plates.

The VASCO team expects most transients to be backgrounds, arising from phenomena such as dust on the glass, specks in the photographic emulsion, or marks made in copying and scanning. To isolate transients from background events, they trained a machine-learning model and assigned probabilities to each candidate representing their probability of being a genuine transient.

### survey-map-caption
Each tile is one plate, coloured by the number of transients on a logarithmic scale. The map is centred on the north celestial pole and keeps areas true, so equal patches of sky occupy equal space. Plates along the Milky Way hold fewer transients, as do plates close to a patch near Taurus and Perseus. Plates shown in grey are present in the catalogue, but were not included in VASCO’s published tests.

### dates-h2
When was each plate taken?

### dates-intro
In the image above, 74 plates on the map are shown in grey. These are plates which were left out of VASCO’s tests, due to a problem with their dates.

The date of each plate is usually taken from the digitised archive at the Space Telescope Science Institute. The archive stores dates in a short block of text called a header, which are attached to each scan.

In our analysis, we have found that around 75% of plates have archive dates which are one day early. Since various statistical tests in the following work rely on knowing when each plate was taken, we have checked the date of every plate to make corrections where required.

### dates-plate-title
Plate labels

### dates-plate-text
The observer wrote each plate’s number and date by hand, and the writing was photographed onto the edge of the glass. Since observing nights spanned two calendar days, observers wrote both dates; for example, “Nov 29/30 1951”.

Astronomers keep time by Universal Time (the clock at Greenwich), which is seven or eight hours ahead of California depending on daylight saving time.

Exposures began after dark, when it was already past midnight in Greenwich, so the second date is the correct one. For the two example plates below, we find that the archive date incorrectly used the first rather than the second date.

### dates-label-caption-1
Plate XE356, marked “E-441 Nov 29/30 1951”. The archive date is incorrectly assigned as 29 November.

### dates-label-caption-2
Plate XE617, marked “E-1595 April 8/9 1956”. The archive date is incorrectly stored as 8 April.

### dates-ast-title
Using asteroid trails to check dates

### dates-ast-text
Asteroids move relative to stars, typically by around half an arcsecond a minute. During long image exposures, asteroids thus appear as short trails, while stars appear fixed since the telescope tracks them. Using NASA’s JPL Horizons service, we can determine the position of every known asteroid on any night in the 1950s with high precision.

For two example plates below, we indicate the trail a known asteroid would follow on the archive’s date (dashed line) and where it would be a day later (solid line). In one day, an asteroid moves around a quarter of a degree, so the two predictions are well separated on the plate.

On the archive’s date, we find no evidence of asteroids in the image. One day later, however, we find trails left by asteroids at the expected position. Across the analysis sample, we found 62 separate occurrences of this. All of them can be explained by a single-day date error in the archive.

### dates-asteroid-1
Plate XE356, showing asteroid (210) Isabella. The archive date for this plate is 29 November 1951, but the asteroid trail is known to occur on 30 November.

### dates-asteroid-2
Plate XE617, showing asteroid (163) Erigone. The archive date for this plate is 8 April 1956, but the asteroid trail is known to occur on 9 April.

### dates-fix-title
Correcting the dates

### dates-fix-text
Palomar numbered the red plates in the order they were taken. The one-day error is not scattered at random: it affects long, unbroken runs of plate numbers, listed in the first table below. Palomar’s own log of the survey also survives. It records each night by its evening date, so the correct date is the log’s date plus one.

Where the log, the asteroids, and the handwritten labels cover the same plates, they agree on 108 of 110. The two exceptions are typing errors in the log, which we correct by hand. For every other plate in a run with the one-day error, we add one day to the archive date.

Eight plates break this pattern, and are listed in the second table below. Six of them are wrong by months or years, and we correct each of them by hand. Every date correction can be found in our code repository.

### dates-why
The mix of date conventions also explains most of the grey plates on the map. VASCO’s tests use only plates taken on one of 370 nights listed in Palomar’s log. The log names each night by its evening date, so a plate taken on the night of 8–9 April is listed under 8 April.

A plate whose archive date is one day early, 8 April, matches the list. A plate with the correct archive date, 9 April, matches only if Palomar also observed on the night of 9–10 April. If it did not, the plate is dropped, which happened to 71 of the 74 grey plates. If it did, the plate is kept but counted as part of the following night, which happened to 88 plates.

### cat-h2
Some plates hold far more transients than others

### cat-intro
The 561 plates in the published tests show a large degree of variation in how many transients they exhibit, ranging from 1 to 4,166.

### cat-spread-title
More variation than expected by random chance

### cat-spread-text
VASCO’s Earth shadow test (see below) scores its result by treating every transient as an independent event: each has the same small chance of falling inside the shadow, whatever the others do.

If transients were spread uniformly over the sky and independently of one another, plates would differ only by counting noise. An average plate holds 178 transients, and the counting noise on that is about its square root (13), so almost every plate would hold between 150 and 205 transients.

However, only 92 of the 561 plates fall in that range. The actual counts spread smoothly from a handful to several hundred: 68 plates hold fewer than 30 transients, while 18 hold more than 500. Some of this variation follows from where each plate points on the sky, as the sky map above shows.

How many transients a plate holds depends on the plate itself, so transients cannot be treated as independent events. A count drawn from only a few plates, such as the nine crossed by Earth’s shadow, can therefore swing much further than counting noise suggests.

### cat-spread-caption
The blue histogram shows the distribution of transient numbers across the 561 plates. The final bar contains the three plates with more than 1,000 transients.

The orange shaded band shows the range expected if transients were spread uniformly and independently. In this case, each plate would exhibit a number of transients close to the average (178), varying only by counting noise.

### cat-edge-title
More transients at the edges of the glass

### cat-edge-text
Across the middle of a plate, transients are spread evenly. Towards the edges they grow denser, reaching two to two and a half times the central density. Random positions that the VASCO team generated as a check show no such rise.

The edges of a plate are set by where the telescope happened to point, and the sky has no reason to be busier there. The rise most likely comes from the plates themselves, or from how they were copied, scanned, and matched. Near the edge, star images blur and can fail to match a modern catalogue.

### cat-edge-caption
Density of transients at each distance from the plate centre, shown relative to the centre density. Blue: transients, with shading showing how far the curve moves when the set of plates is resampled at random (90% of resamples). Grey: random positions. Of the 561 plates, 88% show a higher transient density at the edge of the plate compared to the middle.

### cat-crowd-title
Fewer transients near the Milky Way

### cat-crowd-text
A point of light is kept as a transient only if no modern catalogue lists anything within 5 arcseconds. Near the Milky Way, stars are so crowded that an unrelated star often lies close, and the point is discarded.

A typical plate within 5 degrees of the Milky Way holds 27 transients. Counts rise steadily away from the plane, and beyond about 40 degrees a typical plate holds around 190, seven times as many.

Dust in space also appears to influence the number of transients. In the sky map shown above, the dark patch towards Taurus and Perseus lies behind nearby clouds of interstellar dust in our own Galaxy.

The 49 plates that look through the most interstellar dust hold a median of 34 transients, against about 200 for the 152 plates that look through the least. However, the dust accounts for only part of this patch of sky, and we have not found an explanation for the remaining deficit of transients in this region.

### cat-crowd-caption
Transients per plate against distance from the plane of the Milky Way. At each distance, the line marks the median of all plates within 6 degrees of it, and the shading spans the middle half of those plates.

### shadow-h2
Do transients avoid Earth’s shadow?

### shadow-intro
An object in orbit is seen only by the sunlight it reflects. When it passes into Earth’s shadow, it goes dark.

VASCO looked at objects in geostationary orbit, 36,000 kilometres up, where one orbit takes a day and a satellite stays above one spot on the equator. At that height, Earth’s shadow covers a disc of sky about 17 degrees across, centred on the point directly opposite the Sun. If some transients are glints from objects in such orbits, the disc should hold fewer transients than the same area of sky elsewhere.

Flaws in the photographic emulsion, dust on the glass, and marks from copying have no connection with where Earth’s shadow fell. A shortfall inside the disc would therefore point to something in orbit. VASCO found such a shortfall, and we find it too. Where we differ is in how likely it is to arise by chance.

### shadow-theirs-headline
About 70 trillion to one
VASCO’s figure of 7.6 standard deviations, treating every transient as independent.

### shadow-ours-headline
About one in five
The chance of a shortfall this large, allowing for how much plates differ.

### shadow-step1-title
Nine plates intersect Earth's shadow

### shadow-step1-text
Across the 561 survey plates, Earth's shadow happened to intersect with 9 of them. This figure has been independently confirmed by Brian Doherty, a co-author of the 2026 paper.

Using the corrected dates we describe above, the parts of those plates falling in Earth's shadow hold 343 transients. Together, those parts cover 2.9 plates’ worth of sky: all of plate XE422, two-thirds of two others, and smaller slices of the remaining six. An average plate holds 178 transients, so if transients were spread evenly across the survey, the shadowed parts would hold about 2.9 × 178 = 521.

### shadow-step1-caption
Of the 561 plates used in the analysis, 9 intersect with Earth's shadow. This includes one plate (XE422) which falls fully within Earth's shadow.

### shadow-step2-title
Measuring chance with the real plates

### shadow-step2-text
VASCO’s expected count allows for how much of each plate the shadow covers. Its significance comes from a separate step, a binomial test that treats every transient as an independent event. With our corrected dates, the shadow should hold 521 transients. Under the binomial model that count would vary by about its square root, 23, and the observed 343 lies nearly 8 standard deviations below it, close to VASCO’s 7.6.

Transients are not independent, however: plates differ far more than counting noise allows. The count inside the shadow therefore depends mainly on whether the nine plates it crosses happen to be rich or poor in transients. To measure how far that count can move, we took the nine shadow patches, each with its real shape and position, and laid them on nine plates chosen at random. A randomly chosen plate has no connection with the shadow, so any shortfall in these trials comes from chance alone. Across 20,000 trials, a count of 343 or fewer turned up 21% of the time, which corresponds to 0.8 standard deviations.

Wesley Watters and colleagues used the same method on smaller public parts of the catalogue. Drawing the random plates only from those that pointed close to the shadow on their own night gives a similar result: about one in six, or 1.0 standard deviations.

### shadow-step2-caption
Orange: how the count would vary if every transient were independent (an average of 521, with a standard deviation of 23). Blue: how the count varies when the real Earth shadow patches are laid on random plates. This represents the data-driven null distribution. The observed count of 343 is well below the orange distribution. However, when using the data-driven null distribution, no significant deviation is found.

### shadow-step3-title
Shadow and sunlight on the same plate

### shadow-step3-text
Eight of the nine plates lie partly in shadow and partly in sunlight. Both parts of each plate were exposed at the same time on the same glass, so they share the plate’s own level of transients. This allows a sharper test. If some fraction of transients are glints, that fraction should be absent from the shadowed part of each plate, and only from there. We measured this missing fraction directly.

Across the eight plates, our best estimate of the fraction is zero: transients are as dense in shadow as in sunlight. The comparison allows for the rise in transients towards plate edges, because the expected count in each shadowed part is weighted by where transients usually fall on a plate, not by area alone. Assuming an even spread instead changes the estimate by only 0.01.

Most of the pooled shortfall of 343 against 521 comes from three sparse plates. Two of them (XE617 and XE369) are sparse in their sunlit parts too, and the third (XE422) lies wholly in shadow, so it has no sunlit part to compare.

With only eight plates, the uncertainty is wide. The fraction cannot be negative, and our estimate sits at zero, so the result is best stated as an upper limit. We set it using the method that Gary Feldman and Robert Cousins devised for measurements close to such a boundary. At 90% confidence, fewer than 36% of transients vanish in the shadow. A Bayesian analysis gives almost the same limit.

VASCO’s shortfall implies that about a third of transients vanish. That falls just inside our 90% range, so we cannot rule it out. But if a third of transients did vanish in the shadow, a result as low as ours would happen only about one time in 14.

### shadow-step3-caption
Our measurement of the fraction of transients that vanish inside Earth’s shadow, from the eight plates that lie partly in shadow. The curve shows how well each fraction fits the data, relative to the best fit. The most likely fraction is zero, and the shading covers the values allowed at 90% confidence, up to 0.36. The orange line marks the fraction implied by VASCO’s shortfall, about 0.32.

### shadow-step4-title
Twenty-two more plates in the south

### shadow-step4-text
VASCO’s catalogue stops at the celestial equator: its southernmost plates are centred within about a degree of it, although the survey went on to 30 degrees south. Palomar’s log shows that Earth’s shadow crossed 22 of the plates further south, holding about 9.5 plates’ worth of shadowed sky. This extends the shadowed sky coverage by around a factor 4, compared to using only the nine plates from VASCO. VASCO’s shadow analysis is restricted to the northern hemisphere, without a stated reason, a gap also noted by Watters and colleagues.

We downloaded the scans of these 22 plates, and performed a search for transients following the published VASCO selection: a point of light with the shape of a star, well above the noise, and with no Gaia or Pan-STARRS source within 5 arcseconds.

To check that our search matches VASCO’s, we ran it on two plates from VASCO’s catalogue. We set its thresholds using the first plate only, and left them unchanged for the second. On both plates, our search finds within 6% of the number of transients in VASCO’s catalogue, but not always the same ones: 81% of VASCO’s transients on the first plate are also in ours, and 61% on the second.

The transients that differ lie close to the thresholds on brightness and shape, where two programs that measure slightly differently can place a faint source on opposite sides of the line. [Draft: how much this matters, from rerunning VASCO’s nine shadow plates through our search, once processed.]

### shadow-step4-caption
The plates crossed by Earth’s shadow: VASCO’s nine in the north, and the 22 southern plates, coloured by the number of transients our search finds. Each red ring is the shadow on the night that plate was taken.

### shadow-step5-title
The southern plates on their own

### shadow-step5-text
[Draft, numbers to follow once all plates are processed.]

### shadow-step5-caption
The fraction of transients that vanish inside Earth’s shadow, measured on the southern plates that lie partly in shadow.

### shadow-step6-title
A tighter limit

### shadow-step6-text
[Draft, numbers to follow once all plates are processed.]

### shadow-step6-caption
The fraction of transients that vanish inside Earth’s shadow, measured on all the plates that lie partly in shadow, northern and southern.

### shadow-verdict-title
No evidence of a shadow effect

### shadow-verdict
Once the differences between plates are allowed for, the shortfall inside the shadow is not significant: 0.8 standard deviations, rather than 7.6. On the eight plates the shadow only partly covers, transients are as dense in shadow as in sunlight. Nine plates are too few to rule out a modest effect, with up to about a third of transients vanishing, but they give no evidence for one.

### nuclear-h2
Do transients cluster around nuclear tests?

### nuclear-intro
During the sky survey, the United States tested nuclear weapons in Nevada and the Pacific. Of the survey’s observing nights, 26 fell within a day of a US test. VASCO found that these nights held 1.84 times as many transients as the rest, and estimated the chance of this happening by accident at about one in forty. VASCO does not suggest that the explosions caused the transients; the claim is that the two coincide in time.

We find the same excess. However, the nuclear tests were not spread evenly throughout the year. Most tests took place in spring, when the telescope pointed at parts of the sky where the survey picks up more transients. Once the comparison allows for this, the excess is no longer significant.

### nuclear-theirs-headline
About one in forty
VASCO’s chance for the excess on nights within a day of a test (2.0 standard deviations).

### nuclear-ours-headline
About one in six
The same test, after removing the effect of where the telescope pointed (1.0 standard deviation).

### nuclear-step1-title
Transients track where the telescope pointed

### nuclear-step1-text
The sky above Palomar changes through the year as Earth orbits the Sun. On spring nights the telescope faced away from the Milky Way: plates taken from March to May lie 46 to 58 degrees from it on average, against about 20 degrees from July to December.

Away from the Milky Way, there are fewer stars. A point of light is discarded if a catalogued star lies within 5 arcseconds of it, so fewer points are discarded there, and more are kept as transients. Spring nights therefore yield the most transients of the year: 435 per night in March, against 157 in July. The link is not exact: November and December nights hold more transients than their distance from the Milky Way alone would suggest.

### nuclear-step1-caption
By month: transients per night (blue) and the average distance of that month’s plates from the Milky Way (green). Corrected dates (described earlier) are used.

### nuclear-step2-title
Most nuclear tests fell in spring

### nuclear-step2-text
Most US tests during the survey took place in spring. Palomar observed within a day of a test on 17 nights, and 12 of them fall between March and May, when the telescope faced away from the Milky Way. These nights would hold more transients than the average night, even if the nuclear tests had nothing to do with it.

The other five fall in November 1951, during a series of Nevada tests. Those nights hold about as many transients as other November nights: 333 per night, against 308.

### nuclear-step2-caption
By month: number of observing nights falling within a day of a US nuclear test. Tests on dates when Palomar was not observing do not appear. Corrected dates (described earlier) are used.

### nuclear-step3-title
VASCO’s test

### nuclear-step3-text
VASCO’s test works in two stages. Firstly, it determines how many transients the nights within a day of a test hold per night, compared with all other nights. Then it picks the same number of nights at random, thousands of times over, and counts how often random nights do as well. If random nights rarely do as well, the excess is unlikely to be due to chance.

We ran the same test, counting each transient once rather than weighting it by VASCO’s machine-learning score. With the archive dates VASCO used, the 26 nights near tests hold 1.69 times the usual number of transients. We find that random nights do as well about one time in 34, or 1.9 standard deviations, close to VASCO’s one in forty (2.0). With our corrected dates (see above), 17 nights fall within a day of a nuclear test. They hold 1.65 times the usual number of transients, with random nights doing as well around one time in 18 (1.6 standard deviations). Even with corrected dates, the excess is thus still present.

The 2026 paper also finds the excess strongest on the day before a test. With corrected dates, we find the strongest day to be the day of the test itself.

### nuclear-step3-caption
Testing whether randomly chosen nights (grey) show higher levels of transients than those on nights close to nuclear tests (red).

### nuclear-step4-title
VASCO’s test with the Milky Way effect removed

### nuclear-step4-text
Due to telescope direction, spring nights hold more transients than other nights (see {1_circle_symbol} above). Given that most nights near nuclear tests fell in spring (see {2_circle_symbol} above), random nights drawn from across the whole year are not an appropriate yardstick.

A fairer comparison is the number of transients each night should hold, based on where the telescope pointed. Using a curve fitted to all 561 plates, we assign every plate an expected transient count based on its distance to the Milky Way.

Palomar usually photographed several fields a night, up to eight, so a night’s expected count is the sum over all of its plates.

Measured against expected counts, the nights near nuclear tests hold 1.19 times their expected number rather than 1.65. Random nights match them about one time in six (1.0 standard deviation).

Other search windows do not fare any better. We searched every window from a few days to two weeks either side of the tests. Once allowing for that search (look-elsewhere effect), no test reaches more than around 1 standard deviation.

### nuclear-step4-caption
Once the Milky Way effect is removed, the nights near nuclear tests move into the body of the random sets: one in six random tests does as well, which is not statistically significant.

### nuclear-verdict-title
No evidence of a link to the nuclear tests

### nuclear-verdict
Nights close to US nuclear tests do indeed hold more transients. However, most of this excess comes from the season. Nuclear tests fell mainly in spring, when the telescope was pointed away from the Milky Way (fewer stars, so more transients kept).

Once the Milky Way effect is removed, VASCO’s own test finds no significant link between test dates and transient counts. The excess falls from 65% to 19%, which random nights match about one time in six.

### align-h2
Do transients form geometric patterns?

### align-intro
A flat object tumbling in orbit would flash each time it caught the Sun, leaving a line of points during a long plate exposure. VASCO searched for such alignments and shortlisted five candidates. The best was quoted at 3.9 standard deviations, corresponding to around a one in 20,000 chance of happening by accident.

However, a catalogue of 100,000 unevenly spread points could contain many chance alignments. Two questions naturally arise:

- Do the alignment quality checks used introduce biases in the directions of the alignments they can find?
- Does the real transient catalogue contain more alignments than expected by random chance?

### align-theirs-headline
About one in 20,000
VASCO’s chance that randomly scattered points would form its best alignment, within the small area searched around it.

### align-ours-headline
About one in four
The chance that randomly placed transients would line up as often as the real ones do, across the whole survey.

### align-step1-title
A search that sees every direction

### align-step1-text
VASCO kept a set of transients only if its two sky coordinates, right ascension (east–west) and declination (north–south), rise and fall together almost perfectly: a correlation coefficient (Pearson’s r) above 0.99 in size. That works for diagonal lines. However, along an east–west line, declination does not change, so there is nothing to correlate, and the line is rejected; the same happens to a north–south line.

In simulations, the VASCO method keeps 99% of diagonal lines, 24% of lines 2 degrees from east–west, and 1% of lines running exactly east–west. A geostationary satellite stays fixed in the sky above the observatory while the telescope follows the stars, so on the plate it moves exactly east–west. VASCO’s selection rejects almost all of the tracks such an object would leave.

Particle physicists face the same task when they reconstruct the tracks of charged particles from hits in a detector, and the standard method has no blind direction. Fit a straight line through each set of points that could form a track (a track candidate), and measure how far the points lie from it in the perpendicular direction. That distance is the same whichever way the line runs. Applied to the same simulated lines, a straight-line fit keeps 97–99% of lines in every direction, east–west included.

### align-step1-caption
The share of simulated straight lines kept by each method, by the line’s direction on the sky: north at the top, east and west at the sides. VASCO’s selection (orange) collapses along east–west, where a geostationary satellite’s flashes would lie, and along north–south. A straight-line fit (blue) keeps nearly all lines in every direction.

### align-step2-title
Fitting every track candidate

### align-step2-text
VASCO’s five shortlisted alignments came from an earlier list of 298,165 candidate transients (Solano et al. 2022). None of them appears in the later, cleaner catalogue of 107,875 transients that VASCO uses for its other tests, so we searched that catalogue afresh.

We treated every set of three to six transients up to 15 arcminutes long as a track candidate, fitted a straight line through it, and recorded how closely its points sit to the line. Two sources closer than about 4 arcseconds, the width of a star’s image on these scans, blur into one, so we ignored candidates with points closer than 10 arcseconds.

We then repeated the same process for 1,000 copies of the catalogue, within which every transient is randomly shifted by up to 20 arcseconds. The random shifting breaks any genuine linear relationships, but keeps the transients bunched where they are. Such randomly varied copies of the data can be used to measure the number of track candidates that chance alone produces.

The position of each transient on a scan is uncertain by about an arcsecond. So even points from a genuine straight track will stray from the line by about that much. Each fit is scored with a χ² per degree of freedom (χ²/ndof), which compares the points’ actual distances from the line with that one-arcsecond uncertainty. A genuine track gives a value near 1 or below. Points scattered loosely about the line give a much larger value.

Genuine glints from objects in orbit would form a highly linear pattern, so they would show up as an excess of candidates with χ²/ndof around 1 or below. The real catalogue shows no such excess, for track candidates built from three to six transients (see chart below).

Bin by bin, three-point candidates scatter about chance somewhat more than the randomly shifted copies do, in both directions. The random shifting does not capture every detail of how transients sit together, but the scatter adds up to no excess of good fits.

### align-step2-caption
For track candidates of every size, from three to six transients, the real catalogue (points) follows the null distribution (grey). The differences between the true catalogue and null distribution (bottom row) stay mostly within two standard deviations. There is also no excess among the best track fits with χ²/ndof ≤ 1, where genuine tracks would appear.

### align-verdict-title
No more straight lines than chance

### align-verdict
The track-fit method we propose can see in every direction, and avoids rejecting the east–west tracks that a geostationary object would leave.

Applied to the full catalogue of transients, it finds no evidence of an excess of well-fitting straight lines beyond what chance produces.

### end-h2
Summary

### card-4-question
Are the plates correctly dated?

### card-4-ours
Around 75% of the plates have an archive date that is one day early. Asteroid trails, Palomar’s own logs, and the handwritten dates on plates independently corroborate this finding. We provide a correction to the dates in our code repository.

### card-1-question
Do transients avoid Earth’s shadow?

### card-1-theirs
Far fewer transients fall inside the shadow than expected, with odds of about 70 trillion to one against chance.

### card-1-ours
The shadow reaches only nine plates. Counted plate by plate, a shortfall this large happens about one time in five.

### card-2-question
Do transients cluster around nuclear tests?

### card-2-theirs
Nights within a day of a US test hold 1.84 times as many transients, with a one-in-forty chance of arising by luck.

### card-2-ours
Tests clustered in spring, when the telescope pointed away from the Milky Way and the selection keeps more transients. With that effect removed, VASCO’s own test gives about one in six odds.

### card-3-question
Do transients form linear patterns?

### card-3-theirs
Several short alignments lie on single plates; the best candidate has about a one-in-20,000 chance of being accidental.

### card-3-ours
A straight-line fit that sees in every direction finds no excess of well-fitting straight lines beyond what chance produces.

### end-text
The tests above do not address what each transient is. A 2026 paper by Ivo Busko finds that some share the slight optical blur of nearby stars, which suggests their light came through the telescope; those cases remain open.

A 2024 paper by Hambly and Blair argues that many are defects made when the plates were copied, which would fit the excess observed at plate edges. The shortage of transients towards Taurus and Perseus is also not fully explained.

The 2026 VASCO paper also reports stronger results for the tenth of transients that its machine-learning model rates most likely to be genuine. Those scores are not yet public. The model uses properties of each whole plate, so transients on the same plate share part of their score. Random scores built the same way, with 10% to 30% of their variation shared within a plate, reproduce the paper’s top-tenth shadow result in 7% to 25% of trials. We would be glad to put the final scores through the plate-by-plate tests used here.

### footer-about
Built from the VASCO team’s public catalogue, VASCO’s analysis code in Brian Doherty’s VASCO-ML repository, and our own scripts. Asteroid positions come from JPL Horizons, star positions from Gaia, plate times from Palomar’s log (VizieR VI/25), and test times from the US Department of Energy report DOE/NV-209. Every number here can be reproduced from our repository, which also holds the table of corrected plate dates. We thank Beatriz Villarroel and her colleagues for publishing their data and code; their openness made this work possible.
