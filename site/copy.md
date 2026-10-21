### hero-title
Examining the evidence for glints from orbit in 1950s sky photographs

### hero-lede
Between 1949 and 1957, Palomar Observatory in California photographed the northern sky on glass plates. The plates show faint points of light that have no counterpart in modern sky catalogues. The VASCO project, led by Dr Beatriz Villarroel, has catalogued more than 100,000 of these transients.

Villarroel’s team has proposed that some transients may be produced by sunlight glinting off flat, reflective objects in orbit around Earth in the pre-satellite era. In papers published in 2025 and 2026, the team reported three statistical patterns in support of this theory: fewer transients appear inside Earth’s shadow; more transients appear around the dates of nuclear weapons tests; some transients appear in straight patterns.

We are grateful to the team for publishing its catalogue and code. Building the catalogue took years of careful work, and releasing it openly made this re-examination possible. We repeat the shadow and nuclear-test analyses on their data, and apply both their alignment method and our own to the released catalogue.


### hero-author
PhD in particle physics from the University of Oxford and CERN, with over ten years of research experience on the LHCb experiment. Now a data scientist in industry, following VASCO's work with interest.

### survey-h2
The Palomar Sky Survey

### survey-intro
Images were taken on plates using Palomar’s 48-inch Schmidt telescope, which was built to photograph wide areas of sky. Each plate is a sheet of glass 14 inches square, and covers a patch of sky 6.5 degrees across (about thirteen full Moons). Every patch was photographed through both a blue and a red filter, with red exposures lasting 40 to 60 minutes. Copies of the plates were scanned in the 1990s, and astronomers now work from those scans.

VASCO compared the red plates with three modern catalogues of the sky: Gaia, Pan-STARRS, and WISE. A point of light was kept as a transient if, among other checks, no modern catalogue lists anything within 5 arcseconds of it. An arcsecond is a 3,600th of a degree, and 60 arcseconds constitutes an arcminute. The public VASCO catalogue holds {n:cat.transients} such transients spanning {n:cat.plates} plates.

Many transients are likely to be artefacts, such as dust on the glass, specks in the photographic emulsion, or marks made in copying and scanning. The catalogue as a whole has not been checked by eye. The VASCO team has since trained a machine learning model to assess how likely each transient is to be genuine, but the catalogue was built before this work and includes every transient regardless of its model rating. The machine learning ratings are not yet public, so our analysis herein uses the whole catalogue.

### survey-map-caption
Each tile shown represents one plate, coloured by the number of transients on a logarithmic scale. The map is centred on the north celestial pole and keeps areas true, so equal patches of sky occupy equal space. Plates along the Milky Way hold fewer transients, as do plates close to a patch near Taurus and Perseus. Plates outlined with a dashed line are present in the catalogue but were left out of VASCO’s published tests; we include them in our analysis.

### dates-h2
When was each plate taken?

### dates-intro
In the sky map above, {n:cat.left_out} plates are outlined with a dashed line. These are plates that were left out of VASCO’s tests, due to a problem with their dates.

The date of each plate is usually taken from the digitised archive at the Space Telescope Science Institute. This archive stores dates in short blocks of text called headers, which are attached to each scan.

In our analysis, we have found that around {n:dates.early_pct}% of plates have archive dates which are one day early. Since various statistical tests in the following work rely on knowing when each plate was taken, we have checked the date of every plate and made corrections where required.

### dates-plate-title
Plate labels

### dates-plate-text
The observer wrote each plate’s number and date by hand, and the writing was photographed onto the edge of the glass. Since observing nights spanned two calendar days, observers wrote both dates; for example, “Nov 29/30 1951”.

Astronomers keep time using Universal Time (the clock at Greenwich), which is seven or eight hours ahead of California depending on daylight savings.

Exposures began after dark in California, when it was already past midnight in Greenwich, so the second date shown on the plates is the correct one. For the two example plates below, we find that the archive dates incorrectly use the first rather than the second written date.

### dates-label-caption-1
Plate XE356, hand-labelled “E-441 Nov 29/30 1951” by the observer. In the archive header, the date is incorrectly assigned as 29 November.

### dates-label-caption-2
Plate XE617, hand-labelled “E-1595 April 8/9 1956” by the observer. The archive date is incorrectly stored as 8 April.

### dates-ast-title
Using asteroid trails to check dates

### dates-ast-text
Asteroids move relative to stars, typically by around half an arcsecond per minute. During long exposures, asteroids thus appear as short trails, while stars appear fixed since the telescope tracks them. Using NASA’s JPL Horizons service, we can determine the position of every known asteroid on any night in the 1950s with high precision.

For two example plates below, we show the trail a known asteroid would follow on the archive’s date (red dashed line) and where it would be a day later (red solid line). Over the course of a day, an asteroid moves around a quarter of a degree, so the two positions are well separated on the plate.

Using the archive’s date, we find no trail where the asteroid is predicted to be. One day later, however, we find trails left by asteroids at their expected positions. Across the plates we checked, {n:dates.ast_late_trails} asteroid trails across {n:dates.ast_late_plates} plates appear one day after the archive date. On another {n:dates.ast_ok_plates} plates, {n:dates.ast_ok_trails} trails confirm that the archive date is correct, so the method is able to distinguish between cases.

### dates-asteroid-1
Plate XE356, showing asteroid (210) Isabella. The archive date for this plate is 29 November 1951, but the asteroid trail actually appears on 30 November.

### dates-asteroid-2
Plate XE617, showing asteroid (163) Erigone. The archive date for this plate is 8 April 1956, but the asteroid trail appears on 9 April.

### dates-fix-title
Correcting the dates

### dates-fix-text
Palomar numbered the red plates in the order they were taken. The one-day dating error is not scattered at random: it affects long, unbroken runs of plate numbers, listed in the first table below. Palomar’s own log of the survey also survives. It records each night by its evening date, so the correct date is the log’s date plus one.

Where the log, the asteroids, and the handwritten labels cover the same plates, they agree on {n:dates.log_agree} of {n:dates.log_checked} cases. The two exceptions are typing errors in the log, which we correct by hand. For every other plate in a run with a one-day error, we add one day to the archive date.

Across the {n:dates.dated} plates we dated (which include all {n:cat.plates} in the catalogue), Palomar’s logs date {n:dates.by_log}, asteroid trails {n:dates.by_asteroids}, and the handwritten labels another {n:dates.by_label}. Runs of plate numbers settle a further {n:dates.by_runs}. That leaves {n:dates.unresolved} plates, XE723 and XE884, with no evidence for their date; neither is in the catalogue, so they play no part in our subsequent analysis.

{n:dates.exceptions} plates don't follow this pattern, and are listed in the second table below. {n:dates.gross} of these plates are wrong by months or years, and we correct each of them by hand. All of the date corrections can be found in our code repository.

### dates-why
The dates also explain most of the outlined plates on the sky map above. VASCO’s tests keep only plates whose date matches one of the {n:dates.vasco_nights} nights in Palomar’s log, and the log names each night by the date of its evening: a plate taken during the night of 8–9 April is listed under 8 April.

A plate whose archive date is one day early reads 8 April and matches its night. A plate with the correct date reads 9 April, and matches only if Palomar also observed on the night of 9–10 April. If Palomar did not observe that night, the plate is dropped; this is why {n:dates.dropped} of the {n:cat.left_out} outlined plates are missing from VASCO’s tests. If Palomar did observe that night, the plate is kept but assigned to the wrong night, as happened to {n:dates.shifted} plates.

With corrected dates, all {n:cat.plates} plates match a night in the log, and our analysis uses all of them.

### cat-h2
Number of transients varies greatly across plates

### cat-intro
The {n:cat.plates} plates in the catalogue show a large degree of variation in how many transients they exhibit, ranging from {n:cat.min} to {n:cat.max}.

### cat-spread-title
More variation than expected by random chance

### cat-spread-text
VASCO’s Earth shadow test (see below) scores its result by treating every transient as an independent event: each has the same small chance of falling inside the shadow, regardless of other transients.

If transients were independent of one another and spread uniformly over the sky, plates would differ only by counting noise. An average plate holds {n:cat.avg} transients, and the counting noise on that is around its square root ({n:cat.noise}), so most plates would hold between {n:cat.band_lo} and {n:cat.band_hi} transients.

However, only {n:cat.in_band} of the {n:cat.plates} plates fall within this range. The actual counts spread smoothly, from a handful of transients to several hundred. Some of this variation follows from where each plate points on the sky, as the sky map above shows.

How many transients a plate holds depends on the plate itself, so transients cannot be treated as fully independent events. A count drawn from only a few plates, such as the {n:sh.plates} crossed by Earth’s shadow, can therefore vary considerably more than counting noise suggests.

### cat-spread-caption
The blue histogram shows the distribution of transient numbers across the {n:cat.plates} plates. The final bar contains the {n:cat.over1000} plates with more than 1,000 transients.

The orange shaded band shows the range expected if transients were distributed uniformly and independently. In this case, each plate would exhibit a number of transients close to the average ({n:cat.avg}), varying only by counting noise.

### cat-edge-title
More transients at plate edges

### cat-edge-text
Across the middle of a plate, transients are spread evenly. Towards the edges they grow denser, reaching over two times the central density. Random positions that the VASCO team generated as a check show no such rise.

A plate’s edges are set by its size and shape, not by anything in the sky. Wherever the telescope points, real sources have no reason to gather towards the edges of the plate. The rise may come from the plates themselves, and how the catalogue was built. A detection counts as a transient only if no star in the Gaia or Pan-STARRS catalogues lies within 5 arcseconds of it (Solano et al. 2022). Towards the edges of a plate, star images blur and their measured positions are less accurate, so real stars may miss that match and enter the catalogue as transients. Watters et al. (2026) find the same rise towards plate edges and corners in other lists of transients, and find that transients clump together on most plates.

### cat-edge-map-caption
Left: density of transients across the face of a plate, averaged over the {n:cat.plates} plates with each plate counted equally. Density is shown relative to the density near the plate centre (north up, east left).

Right: the random positions VASCO generated as a check, which cover the plate evenly up to its edge.

### cat-edge-caption
Density of transients at each distance from the plate centre, shown relative to the centre density.

Blue: transients, with shading showing how far the curve moves when the set of plates is resampled at random (90% of resamples).

Grey: random positions. Of the {n:cat.plates} plates, {n:cat.edge_share}% show a higher transient density at the edge of the plate compared to the centre.

### cat-crowd-title
Fewer transients near the Milky Way

### cat-crowd-text
A point of light is kept as a transient only if no modern catalogue lists anything within 5 arcseconds. Near the Milky Way, stars are more crowded, causing a greater proportion of transients to fail the 5 arcsecond threshold.

A typical plate within 5 degrees of the Milky Way holds {n:crowd.inner} transients. Counts rise steadily away from the galactic plane, and beyond 40 degrees a typical plate holds around {n:crowd.outer}, {n:crowd.ratio} times as many.

### cat-crowd-caption
A smooth curve fit to the relationship between transient counts and distance from the galactic plane. Use the buttons to move plates between their position on the graph, and their position in the sky map.

### cat-pointing-title
Where a plate points shapes how many transients it keeps

### cat-pointing-text
The number of transients retained on a plate depends on where it pointed. Crowding near the Milky Way and dust both remove detections, while plates near the ecliptic, the band of sky that asteroids follow, gain some from asteroids missing from modern lists. These effects come from how the catalogue was built, and act on transients regardless of their origin.

Dust between the stars is also linked to fewer transients. Clouds of dust in our own Galaxy dim everything behind them, and the sparse patch near Taurus and Perseus on the sky map (see above) lies behind some of the closest dust clouds.

To separate dust from crowding, we compared plates more than 20 degrees from the Milky Way, where stars are not crowded. The {n:dust.n_dusty} plates that look through the most dust hold a median of {n:dust.median_dusty} transients, compared to {n:dust.median_clear} held by the {n:dust.n_clear} plates which look through the least dust. Since dust does not cover the full Taurus–Perseus patch, we have not found an explanation for the shortfall across the full region.

One further effect is tied to particular plates rather than to the sky. The {n:cat.eq_plates} plates centred on the celestial equator, the southern edge of the catalogue, hold about {n:cat.eq_share}% as many transients as other rows of plates once the factors above are allowed for. Where these plates overlap their northern neighbours, the same patch of sky yields only about {n:cat.eq_same_sky}% as many transients on the equatorial plates. The cause therefore lies in how these plates were processed, not in the sky. The shortfall is already present in the list published by Solano et al. (2022), and is not explained by any of the recorded properties of the plates.

### cat-pointing-fit
To assign each plate an expected number of transients based on pointing alone, we combine the above factors in a single fit across all plates. In this model, dust is averaged across each plate, and we include a term for the equatorial row. In the following sections on Earth’s shadow and nuclear tests, we use these expected counts to compare plates on equal terms. Where appropriate, the model is refitted without the plates being tested, so that they cannot shape their own expected counts.

### cat-pointing-caption
Averaged over groups of similar plates, the pointing model predicts transient counts well. Pointing accounts for {n:cat.pointing_share}% of the differences between plates, and the shadow and nuclear tests remove this component, so that it cannot confound their results.

### cat-rows-caption
For plates in the equatorial row, we find that the level of transients is significantly lower than expected given their pointing. We thus include this effect as an additional factor in our pointing model.

### shadow-h2
Do transients avoid Earth’s shadow?

### shadow-intro
An object in orbit is visible only via the sunlight it reflects. When orbiting objects pass into Earth’s shadow, they should therefore go dark.

VASCO considered objects in geostationary orbit (36,000 kilometres up), where a single orbit takes one day and a satellite stays above a point on the equator. At this height, Earth’s shadow covers a disc of sky around 17 degrees across, centred on a point directly opposite the Sun. If some transients are due to glints from objects in orbit, the shadow disc should exhibit fewer transients than the same area of sky elsewhere.

Flaws in the photographic emulsion, dust on the glass, and marks from copying have no connection with where Earth’s shadow fell. If fewer transients fall inside the shadow than the differences between plates can explain, that would point to something in orbit. We find the same shortfall of transients in Earth's shadow as VASCO. Where we differ is on the statistical interpretation of the deficit.

### shadow-step1-title
Plates intersecting Earth's shadow

### shadow-step1-text
Across the {n:cat.plates} plates in the catalogue, Earth's shadow happened to intersect with {n:sh.plates_n} of them. {n:sh.vasco_plates_cap} of these are among the plates in VASCO's tests, as Brian Doherty, a co-author of the 2026 VASCO paper, also found. The remaining plate, XE472, was taken on 18–19 November 1949, two nights before VASCO's list of nights begins.

Throughout this section, we use the central 6 by 6 degrees of each plate, matching the plate size VASCO use in their own analysis. VASCO’s random positions cover this square evenly on every plate, as shown above.

Using our corrected dates, the areas of those plates falling in Earth's shadow are found to hold a total of {n:sh.observed} transients.

Together, the shadowed areas cover {n:sh.worth} plates’ worth of sky (all of plate XE422, two-thirds of two others, nearly half of XE472, and smaller slices of the remaining six). An average plate holds {n:sh.avg} transients in this square, so if transients were spread evenly across the survey, the shadowed parts would be expected to hold about {n:sh.expected} transients ({n:sh.worth} × {n:sh.avg}).

VASCO’s archive dates place the shadow slightly differently. A date one day out moves the shadow by only about a degree, small compared to to its 17-degree width, so on the {n:sh.vasco_plates} plates in their tests the count is similar: {n:sh.vasco_obs} against {n:sh.vasco_exp} expected, a shortfall of about a third.

### shadow-step1-caption
Of the {n:cat.plates} plates in the catalogue, {n:sh.plates_n} intersect with Earth's shadow. This includes one plate (XE422) which falls fully within Earth's shadow.

### shadow-step2-title
Measuring the probability of this deficit

### shadow-step2-text
VASCO scores the shortfall with a binomial test, which treats every transient as an independent event. Their primary figure of 22 standard deviations compared the shadow with its share of the whole sky; their 7.6 standard deviation result used each plate’s actual coverage, as we do here. With that model and our corrected dates, the expected count of {n:sh.expected} would vary by about its square root ({n:sh.binom_sd}), which puts the observed {n:sh.observed} about {n:sh.binom_z} standard deviations below expectation.

Transients are not independent, however: as shown above, transient counts on plates differ considerably more than counting noise permits. The count inside a given shadow therefore depends on whether the plate it falls on happens to be rich or poor in transients. To measure how far the count can move, we took the {n:sh.plates} Earth shadow patches, each with its real shape and position, and laid them onto {n:sh.plates} plates chosen at random.

A randomly chosen plate has no connection to Earth's shadow, so any shortfall in such trials comes from chance alone. Plates also differ in where they point, so in each trial we scale the count to what the chosen plate would hold had it pointed where the shadow plate did, using the pointing curve from the previous section (fitted without the plates the shadow crosses). Across 20,000 such trials, a shadow transient count of {n:sh.observed} or fewer turned up {n:sh.p_point}% of the time, a statistical significance of {n:sh.sig_point} standard deviations. Without the allowance for pointing, the figure is {n:sh.p_random}%.

Where the shadow fell on plates also matters. Earth’s shadow falls near the point opposite the Sun, which lies on the ecliptic, so it tends to cross plates in the thinnest parts of the catalogue. Of the {n:sh.plates} shadow plates, {n:sh.taurus_lc} lie towards Taurus and Perseus, where the catalogue holds unusually few transients, and {n:sh.eq_plates} are on the equatorial row (see the section on transients per plate).

The pointing curve allows for the equatorial row and for dust across each plate. Even so, {n:sh.empty_n} of the shadow plates, {n:sh.empty_plates}, hold under a tenth of what their pointing predicts, and none of our modelled factors explain why. Wesley Watters and colleagues made the same point, and found with their own random-shadow test, on a larger published list of transients that overlaps VASCO’s, that a deficit this large arises by chance about 30% of the time. The shadow plates were all taken 10 to 45 degrees from the point opposite the Sun. Drawing the random plates only from that range gives a probability of about one in {n:sh.odds_near} ({n:sh.sig_near} standard deviations).

VASCO also repeated their test using only transients within 2 degrees of each plate’s centre, away from the dense edges. In that region we count {n:sh.centre_obs} transients in shadow against {n:sh.centre_exp} expected, a shortfall random plates match about one time in {n:sh.odds_centre}.

Seen from Palomar rather than from the centre of the Earth, the shadow at geostationary height shifts by several degrees and crosses {n:sh.topo_plates} plates instead of {n:sh.plates}. This version of our analysis yields {n:sh.topo_obs} transients against {n:sh.topo_exp} expected, a shortfall random plates match about one time in {n:sh.odds_topo} ({n:sh.sig_topo} standard deviations).

### shadow-step2-caption
Comparing the observed count of transients in Earth's shadow to two measurements of chance. Use the arrows to step through how each result is built.

Purple: VASCO’s model, in which every transient is independent.

Blue: the real Earth shadows laid on random plates across {n:sh.trials} trials.

### shadow-plates-title
Taking a closer look at the plates in shadow

### shadow-plates-text
Earth’s shadow crosses {n:sh.plates} of the catalogue’s plates. One of them, XE422, lies fully inside Earth's shadow, so it has no sunlit area to compare with. The other {n:sh.partial} plates are shown below, with the position of every transient from VASCO’s catalogue indicated.

The plates differ significantly in terms of total transient counts. XE415 holds {n:pl.XE415_n} transients, while XE355 only holds {n:pl.XE355_n}. Within a plate, transients bunch into clumps and gather along some edges, while some plates are almost empty in places. For example, XE640’s transients crowd along its northern edge and are missing from two large rounded regions. The sky has no reason to be empty in two round patches of a single plate, so these gaps most likely come from the plate itself, or from how it was scanned and searched.

Earth’s shadow splits each plate in two. On XE617 the shadow covers {n:pl.XE617_area}% of the plate, and holds {n:pl.XE617_O} of its {n:pl.XE617_n} transients. On XE640 it clips {n:pl.XE640_area}% of the plate, in a corner that holds none of its {n:pl.XE640_n} transients.

### shadow-plates-caption
The {n:sh.partial} plates that Earth’s shadow partly covers, showing the central 6 by 6 degrees of each. Each dot is a transient in VASCO’s catalogue. Transients in Earth's shadow are shown in red, and those in sunlight blue. The red line marks the edge of Earth's shadow, 8.5 degrees from the point opposite the Sun when the plate was taken.

Under each plate, “expected” shows how many transients the shadowed region would hold if transients there were as dense as in the sunlit region, allowing for the higher density towards plate edges.

The ratio is the observed number of shadow transients divided by expected. A ratio of 1 means no difference in transient density between the shadow and sunlit regions.

### shadow-fit-title
Fitting the shadow transient density

### shadow-fit-text
VASCO report that transients are rarer in Earth’s shadow. Expressed as a measurement: in the shadowed part of a plate, transients are r times as dense as in the sunlit part. With no Earth's shadow effect, the ratio r should be 1. If some transients genuinely vanish in the shadow, r should be below 1.

In practice, r compares two shares on each plate: the share of its transients that lie in the shadowed part, and the share the shadowed part would hold with no effect. Because transients are denser towards the edges and corners of every plate (see the plate maps above), this second share is not simply the shadowed area. We weight each part of the plate by the average density of transients there, measured only on the {n:sh.never} plates that Earth’s shadow never crossed, so a real shortfall in the shadow cannot affect the weights.

{n:fit.ex_plate}, for example, has {n:fit.ex_area}% of its area in shadow, which with no effect would hold {n:fit.ex_weighted}% of its transients. It actually holds {n:fit.ex_share}% ({n:fit.ex_O} of {n:fit.ex_n}). If transients vanish in the shadow, their share should fall below the expected share, and the ratio r is designed to measure by how much.

Because each plate is compared only with itself, it does not matter how many transients a plate holds in total. The variation between plates, large as it is, plays no part here.

The single value of r that best fits all {n:sh.partial} plates together is {n:fit.r}: transients are as dense in shadow as in sunlight. If every transient were an independent event, counting statistics would put r between {n:fit.count_lo} and {n:fit.count_hi}.

The plates themselves show that this uncertainty is too optimistic: single-plate ratios run from 0 on XE640 to {n:fit.max_r} on {n:fit.max_plate}, a much wider spread than counting allows. In {5_circle_symbol} we measure the true uncertainty, using trials built from the catalogue’s other plates.

### shadow-fit-caption
Grey points: the density of transients in shadow divided by the density in sunlight, for each plate that intersects Earth's shadow. Error bars are taken from counting statistics. Plates with no transients in shadow have only upper limits.

Blue: the best-fit ratio across all {n:sh.partial} plates, with a one-standard-deviation range taken from counting statistics.

### shadow-toys-title
Finding the true uncertainty

### shadow-toys-text
Counting statistics treat every transient as independent. To measure the true uncertainty, we ran the same fit on plates where no shadow effect is possible. Each trial lays the {n:sh.partial} real shadow outlines onto {n:sh.partial} other plates that Earth’s shadow never crossed, and counts the number of transients inside and outside the shadow. In these trials, random plates are selected to contain a similar number of transients to the corresponding plates in the real analysis.

Over 20,000 such trials, the fit gives a typical r value of {n:fit.null_median}, with the middle 68% of results falling between {n:fit.null_lo} and {n:fit.null_hi}. Our best fit result on the {n:sh.partial} genuine plates (r = {n:fit.r}) falls well within this range, indicating that there is no evidence of reduced transient signals in Earth's shadow beyond random chance.

The same trials show that the fit can detect a genuine shortfall. When we remove transients from inside the shadow outlines on the random plates, so that they are only {n:fit.inj_r} times as dense there (a shortfall as large as VASCO’s), the fit returns a value as high as {n:fit.r} only {n:fit.inj_pct}% of the time. VASCO did not measure the effect plate by plate, but a shortfall of the size their count implies is excluded at {n:fit.inj_cl}% confidence.

These results hold when the random plates are selected differently, when each shadow is instead moved around on its own plate, or when any one plate is left out.

### shadow-toys-caption
Blue: fit result to the real plates.

Grey: 20,000 trials using the real shadows applied to random plates.

Purple: trials with a shortfall of transients in Earth's shadow added to the data.

Use the arrows to step through how this comparison is built.

### shadow-verdict-title
No significant shadow deficit

### shadow-verdict
Judged against a data-driven null distribution rather than counting noise, the pooled shortfall falls from {n:sh.binom_z} standard deviations to {n:sh.sig_point}. Within each plate, transients are also found to be no rarer in shadow than in sunlight: the best-fit ratio is {n:fit.r}, and at 90% confidence level the true ratio lies between {n:fit.ci90_lo} and {n:fit.ci90_hi}.

Earth’s shadow also crosses {n:sh.south} plates south of the celestial equator, beyond the edge of VASCO’s catalogue. Adding these plates would allow the same test to be performed on around {n:sh.south_ratio} times as much shadowed sky.

### nuclear-h2
Do transients cluster around nuclear tests?

### nuclear-intro
During the sky survey, the United States tested nuclear weapons in Nevada and the Pacific. In two papers (Scientific Reports, 2025 and 2026), VASCO report that observing nights within a day of a US test held more transients than other nights. The 2026 paper puts the chance of this at about 1 in 40 (p = 0.024). The authors fixed this time window before looking at the data, and read the sharpness of the association in time as a sign that it is genuine.

Repeating their measurement, we find a similar excess. It reduces to a level consistent with chance once two further factors are taken into account:

- Nights near tests tended to have more photographs, and each photograph adds its own transients.
- On nights close to tests, the telescope tended to point at parts of the sky where more transients survive the catalogue’s selection, mostly away from the Milky Way. Our pointing model (see the section on transients per plate) expects these plates to hold about {n:nc2.exp_pct}% more transients than others.

The choice of time window also matters. One day either side of a nuclear test is one of several reasonable windows that could have been chosen. Changing the window moves the result up or down, as individual plates enter or leave the window. We find that the one-day window does not stand out among alternative choices.

### nuc-vasco-title
Repeating VASCO’s test

### nuc-vasco-text
VASCO’s test counts the transients on each observing night, compares the nights within a day of a US test with all other nights, and asks how often the same number of nights selected at random show at least as many transients. We ran the analysis counting each transient once, not weighted by VASCO machine-learning scores, as the scores were not available at the time of analysis.

On VASCO’s own construction, using the default archive dates and their list of test nights, the {n:nuc.vasco_nights} nights near tests hold {n:nuc.vasco_ratio} times the usual number of transients. Random nights exhibit at least as many transients around one time in {n:nuc.vasco_odds} ({n:nuc.vasco_sig} standard deviations), close to VASCO’s 1 in 40 result.

When we use our corrected dates (as described above), {n:nuc.nights} nights fall within a day of a test. They hold {n:nuc.ratio} times the usual number of transients, with random nights doing as well about one time in {n:nuc.odds} ({n:nuc.sig} standard deviations).

### nuc-vasco-caption
The excess of transients on the {n:nuc.nights} nights within a day of a US test (red line), compared to the same number of nights drawn at random from the survey, 20,000 times (grey). The analysis is performed using corrected dates.

### nuc-corr-title
Counting photographs, and considering telescope pointing

### nuc-corr-text
VASCO’s test counts transients per night, but a night’s transients come from its photographs. The {n:nuc.nights} nights near nuclear tests averaged {n:nuc.plates_window} plates per night, compared to {n:nuc.plates_other} on other nights. Part of this gap comes from nights with no plates in the catalogue: {n:nc2.empty_near}% of the nights near tests, against {n:nc2.empty_other}% of other nights. A night with no plates contributes zero transients, even though nothing was observed. Watters et al. (2026) make the same point: counted over only the days with plates north of the celestial equator, the link between transients and test days is no longer significant. Counted per plate rather than per night, the excess reduces from {n:nc2.night_r} to {n:nc2.photo_r}.

The sky above Palomar changes through the year as Earth orbits the Sun. On spring nights, the telescope faced away from the Milky Way, where there are fewer stars to rule transients out. Plates within 10 degrees of the Milky Way hold a median of {n:nc2.med_low} transients, compared to {n:nc2.med_high} in plates more than 40 degrees away. Of the plates near nuclear tests, {n:nc2.spring}% were taken between February and June, which on their own would contribute around {n:nc2.exp_pct}% more transients.

A like-for-like comparison uses the pointing curve detailed earlier: the number of transients each plate’s pointing alone would leave in the catalogue. Measured against expected counts, plates near nuclear tests hold {n:nc2.point_r} times their expected number, rather than {n:nc2.photo_r}. Fitting the curve without the plates near tests, so they cannot shape their own expected counts, also gives {n:nc2.held_out}.

To judge whether {n:nc2.point_r} is an unusual excess, it must be compared with what chance can produce. VASCO’s test compares the nights near tests with the same number of nights drawn at random from across the survey. Using corrected dates, we find that random nights show at least as many transients {n:nc2.point_pv}% of the time.

This comparison treats the {n:nuc.nights} nights near tests as independent of one another, but in reality they are not. The tests came in series, so the nights near them bunch into a few observing runs. Plates from the same run share batches of emulsion and similar observing conditions. A clustered set of nights varies more by chance than a scattered one, so random nights understate how often chance can produce at least as many transients.

Our comparison retains this clustering property of the data. We build a null distribution for the test using the real test calendar:

- Move every test date by the same amount, for example two years and three weeks later. The result has the same gaps between tests and the same times of year, but falls on dates when no tests took place.
- Repeat for every shift of one to five years, earlier or later, plus up to two months. This gives a total of {n:nc2.slides} shifted calendars.
- Leave out every plate taken within two weeks of a real test, so that no shifted calendar picks up plates from around the real tests.
- Score each shifted calendar exactly as the real one: per photograph, against what each plate’s pointing predicts.

About one in {n:nc2.odds} shifted calendars shows an excess at least as large as the real one ({n:nc2.point_pc}%, or {n:nc2.sig} standard deviations). Whether compared with shifted calendars or with random nights ({n:nc2.point_pv}%, or {n:nc2.point_pv_sig} standard deviations), the excess near nuclear tests is not found to be statistically significant.

### nuc-corr-caption
Once the number of photographs and the telescope pointing are accounted for, the plates near nuclear tests hold {n:nc2.point_r} times their expected number of transients.

Compared to a null distribution built from {n:nc2.slides} shifted test calendars, about one in {n:nc2.odds} random trials exhibit at least the same level of excess.

### nuc-checks-title
Additional factors tested, but do not significantly alter the result

### nuc-checks-text
We also allowed for a number of other measurable properties of each plate: how deep the plate goes (the share of faint Gaia stars it records), how sharp its stars are, exposure time, height above the horizon, position of the Moon, and VASCO’s own plate-quality flag. None of them moves the result by more than a few per cent; with all of them together, the result is {n:nc3.all}.

Adjusting what counts as a transient makes little difference either. Counting 2-arcminute cells rather than transients, keeping only the central part of each plate, or keeping only transients with no others within 5 arcminutes produces results within the range {n:nc3.defs_lo} to {n:nc3.defs_hi}.

Leaving out any single series of nuclear tests gives results within the range {n:nc3.series_lo} to {n:nc3.series_hi}, with similar levels of statistical significance to the default result with all of them included.

### nuc-checks-caption
How far each further check moves the ratio for plates near nuclear tests. Each version of the analysis is measured against its own set of shifted test calendars.

### nuc-global-title
Choice of time window

### nuc-global-text
The 2026 paper describes the association as highly specific in time, and interprets this as a sign of a physical link. However, one day either side of a test is one of several time windows that could equally well have been chosen. Allowing up to a week either side gives {n:nc4.windows} windows, each of which can be analysed in the same way.

Across these windows the results run from {n:nc4.r_lo} to {n:nc4.r_hi}, with values falling into two separate groups. The upper group is driven by three particularly transient-rich plates moving in and out of the windows:

- {n:nc4.rw_plate} ({n:nc4.rw_date}) holds {n:nc4.rw_n} transients where its pointing predicts about {n:nc4.rw_e}. It was taken {n:nc4.rw_days} days after a test, so every window that reaches that far includes it.
- {n:nc4.rd1} and {n:nc4.rd2}, both taken on test days, hold {n:nc4.rd_share}% of the transients on the {n:nc4.rd_plates} test-day plates. They dominate the narrowest windows; wider windows dilute them with ordinary plates.

Without these three plates, the two groups merge into a single group centred on {n:nc4.wo_median}, and VASCO’s one-day time window yields {n:nc4.wo_vasco}. These plates were identified by looking at the data, so this check is intended to describe where the excess comes from rather than testing it.

VASCO fixed their time window in advance, so their result needs no allowance for other windows. But results like it are common among the alternatives. Had the strongest of the {n:nc4.windows} windows been reported, chance alone would produce a result at least as strong as VASCO’s {n:nc4.chance}% of the time, a reflection of the look-elsewhere effect. When viewed alongside the other windows, VASCO’s result does not indicate that the association is highly specific in time.

### nuc-global-caption
Each block shows the result of the analysis run using a specific choice of time window. There are 64 blocks in total, one of which is the one-day window used by VASCO in their analysis. Hover over a block to see the time window it corresponds to.

Using the buttons, you can repeat the analysis with and without a subset of highly transient-rich plates. The bimodal structure is caused by these plates entering and exiting the analysis time window. With the plates removed, the distribution of results collapses to a single distribution centred around no effect.

### nuclear-verdict-title
No significant link to the nuclear tests

### nuclear-verdict
We reproduce VASCO’s excess: nights within a day of a US nuclear test hold more transients than other nights. Counted per photograph, and measured against where the telescope pointed, the excess is found to be {n:nc2.point_r}, which test calendars shifted to other years match about one time in {n:nc2.odds} ({n:nc2.sig} standard deviations). No other measured property of the plates changes this.

VASCO fixed the time window of their analysis in advance, but it is one of several reasonable windows. Results using different time windows are variable, and depend on a few unusually rich plates entering and leaving the window. Without three such plates, identified by looking at the data, every window sits close to no difference. Had the strongest of {n:nc4.windows} time windows been reported, chance alone would give a result at least as strong as VASCO’s {n:nc4.chance}% of the time.

### align-h2
Do transients form geometric patterns?

### align-intro
A flat object tumbling in orbit would flash each time it caught the Sun, leaving a line of points during a long plate exposure. VASCO searched for such alignments and shortlisted five candidates. Their best candidate was reported at 3.9 standard deviations, corresponding to an estimated 1 in 10,000 chance of occurring by chance.

A catalogue of 100,000 unevenly spread points may contain many chance alignments. Three questions naturally arise:

- Do the alignment quality checks used introduce biases in the directions of the alignments they can find?
- Does the real transient catalogue contain more alignments than expected by random chance?
- How many genuine glint tracks could the catalogue hold?

### align-step1-title
Searching for alignments in all directions

### align-step1-text
VASCO kept a set of transients only if its two sky coordinates, right ascension (east–west) and declination (north–south), rise and fall together almost perfectly: a correlation coefficient (Pearson’s r) above 0.99. This approach works well for diagonal lines. However, the closer a line runs to east–west, the less its declination changes along it. The small scatter in each measured position then swamps the trend, the correlation falls below the cut, and the line is rejected. The same happens to lines close to north–south.

In simulations (image below), the VASCO method keeps {n:align.pear_diag}% of diagonal lines, {n:align.pear_2}% of lines 2 degrees from east–west, and only {n:align.pear_ew}% of lines running exactly east–west. A geostationary satellite stays fixed in the sky above the observatory while the telescope follows the stars, so on the plate it moves exactly east–west. Objects at the same height on slightly tilted orbits also move close to east–west. VASCO’s selection would reject almost all of the tracks such an object would leave.

The same alignment task arises in particle physics, where the paths of charged particles must be reconstructed from hits they leave in a detector. The simplest tracking method fits a straight line through a set of points, resulting in a track candidate. The quality of each candidate can be assessed by checking how far the points lie from the line of best fit. Such a method works equally well in all directions, including east–west.

### align-step1-caption
The share of simulated straight lines kept by each method, shown as a function of direction in the sky. VASCO’s selection (orange) collapses along east–west, where genuine geostationary satellite glints would lie, as well as along north–south. A straight-line fit (blue) keeps nearly all lines in every direction.

### align-step2-title
Fitting every track candidate

### align-step2-text
VASCO’s five shortlisted alignments came from an earlier list of 298,165 candidate transients (Solano et al. 2022). Of their {n:align.vasco_transients} transients, {n:align.vasco_in_cat} appears in the later catalogue of {n:cat.transients} transients that VASCO uses for its other tests. We have applied our track-fit method to the full catalogue of {n:cat.transients} transients.

We treated every set of three to six transients up to 15 arcminutes long as a track candidate, fitted a straight line through it, and recorded how closely its points follow the line. Two sources closer than about 4 arcseconds (the width of a star’s image on the scans) blur into one, so we ignored candidates with points closer than 8 arcseconds. The released catalogue itself holds almost no transients closer together than this. {n:align.vasco_close_n} of VASCO’s five shortlisted alignments contain points only {n:align.vasco_close} arcseconds apart, so they could not appear in the catalogue, regardless of the method used to search for alignments. Watters et al. (2026) also find that 8 of the 24 transients in the five alignments lie within 5 arcseconds of sources seen by the NEOWISE infrared survey, so may be stars rather than transients.

The position of each transient on a scan is uncertain by about an arcsecond. As such, even points from a genuine straight track can stray from a true straight line. Each track fit is scored with a χ² per degree of freedom (χ²/ndof), which compares the points’ actual distances from the line with a one-arcsecond uncertainty. A genuine track gives a value near 1 or below. Points scattered loosely about the line give a much larger value.

A good fit alone does not make a track stand out: in crowded parts of a plate, pure chance can cause transients to align frequently. We therefore give each track candidate a track score: the number of tracks of the same length that would arise by chance in a box at least 10 arcminutes across around the candidate, given how many transients fall in the box.

This concept, and the formula we use, aligns with the method VASCO used to judge its own candidates (Edmunds and George, 1985). We apply our method in all directions, measure straightness from the fitted line, and score only track candidates that are straight to within the measurement error (a χ² probability of at least 1%). The lower the track score, the less often chance alone can produce such a track.

Transients often come in close pairs: the catalogue holds far more pairs of transients within 20 arcseconds of each other than chance would give. Two transients this close fix a direction only loosely, so a third lying roughly in line with them makes a short track that looks stronger than it is. Following Lezama and colleagues (2015), who developed a method for finding alignments among clustered points, we count a track’s transients in boxes 20 arcseconds long, the distance over which close pairs are unusually common: transients closer together than this along the track count once.

A track candidate built from six transients also contains all possible shorter track candidates inside it. Where track candidates share a transient, we keep only the best-scoring candidate, so each transient counts once.

### align-method-caption
How the track score rates a real track candidate on plate {n:align.ex_plate}. A candidate is scored only if it first passes the straightness test.

The track score, μ, is the number of tracks with r transients, at most L long and lying within w of a straight line, that we would expect if the N transients in the box were scattered at random, with n their density (Edmunds and George 1985):

{eq:track-score}

### align-step3-title
Measuring how often tracks arise by chance

### align-step3-text
Any catalogue of transients holds some tracks by chance, so finding track candidates proves nothing by itself. The proper test is whether the real data holds more strong track candidates than chance alone would lead to.

This test requires a background-only model: random catalogues which match the real data in every respect, except that they contain no genuine tracks. Since transients cluster on plates, and their surface density varies, a uniform random distribution cannot serve as this model.

We therefore created {n:align.copies} random catalogues based on the real data itself. Each random catalogue keeps all of the plates from the real data, including the number of transients on each and their tight groups of transients, while breaking any alignments. Every lone transient is moved at random by up to {n:align.shift} arcseconds: enough to break genuine lines, but small enough that each plate keeps its distribution of transients on larger scales.

Transients within {n:align.link} arcseconds of one another, directly or through a chain of neighbours, form a group. Each group is moved as one piece, turned at random and shifted by up to {n:align.shift} arcseconds, and no transient is moved to within {n:align.link} arcseconds of another group. A track consisting of three or more transients from one group is left out, both in the real and random catalogues alike, because a group keeps its shape when moved; as such, any line inside it would appear in every random catalogue. The real dataset and every random catalogue are processed using the same code, to find track candidates in an identical fashion.

To measure chance correctly, the random catalogues must meet two conditions:

- <strong>They must break up real tracks.</strong> To test this, we added 100 artificial tracks to the real catalogue and made random catalogues from it in the same way. For every 100 added tracks scoring 1 in 10,000 or better, the random catalogues gain {n:align.absorb4}; at 1 in 1,000 or better, they gain {n:align.absorb3}. When comparing the real data to random catalogues, any genuine excess of glint tracks would therefore still show at {n:align.keep4}% of its true size ({n:align.keep3}% for 1 in 1,000). The limits we set on glint tracks below allow for this partial absorption of signal.
- <strong>They must not understate chance.</strong> If the random catalogues held fewer chance tracks than the real data, chance alone could pass for a signal. In practice, the random catalogues err slightly the other way. Pairs of transients 20 to 30 arcseconds apart are rare in the real catalogue, and moving the transients fills that gap, giving the random catalogues {n:align.pair2030} times as many such pairs and thus slightly more chance tracks. This makes the background-only test conservative: it could slightly hide a real excess, but it cannot create a false one.

### align-vasco-title
Comparing our track score and VASCO’s score

### align-vasco-text
Our implementation of VASCO’s formula reproduces the scores in their published table of candidates. Where both methods find a track candidate, they agree closely (rank correlation of {n:align.rho}), and {n:align.top50} of the 50 top-scoring candidates found by the VASCO method are among the 50 ranked highest by the track score method. However, the methods diverge in two ways.

VASCO’s score rates {n:align.vasco_strong} candidates as 1 in 1,000 or better, counting overlapping ones. The track score sets aside {n:align.vasco_bent} of these, because their transients do not fall along a straight line within the measurement error. The track score also rates candidates in all spatial directions, including those close to north-south and east-west which are rejected by VASCO’s correlation-based screening.

To compare the methods, we measured their track-finding efficiency, ε: the proportion of glint tracks found. We added 100 artificial tracks to the catalogue in a set direction, each made of 4 to 6 transients spread over 2 to 10 arcminutes, and counted how many each method found as a 4-transient track. For a fair comparison, each score requires a pass mark that lets through the same number of false alarms: around {n:align.false_alarms} track candidates created purely by chance across the whole survey, measured using the random catalogues described in {3_circle_symbol}.

Within two degrees of east–west, where a genuine geostationary object’s track would lie, VASCO’s score has ε = {n:align.vasco_2}%, falling to about {n:align.vasco_ew}% along exactly east–west. By comparison, the track score method has ε = {n:align.track_ew}% in every direction.

For diagonal tracks, VASCO’s score is more efficient at the same pass mark: ε = {n:align.vasco_diag}%, against {n:align.track_diag}%. When two transients of an artificial track fall less than 20 arcseconds apart, the track score counts them once, and some tracks are then left with too few transients to be found as a 4-transient track. VASCO’s score counts both, which finds those tracks but also lets the catalogue’s many close pairs pose as strong tracks.

Counting close pairs twice has a cost: VASCO’s result then depends on how the random catalogues treat close pairs. A close pair could be a feature of the plate, which should move as one piece in the random catalogues, or two separate glints, which should move independently. Our method keeps close pairs together; the alternative splits them, moving every transient on its own by up to {n:align.alone_move} arcminutes.

For VASCO’s score, the two choices yield opposite answers. With pairs kept, chance gives {n:align.vasco_exp} tracks scoring 1 in 1,000 or better, and the real catalogue’s {n:align.vasco_real} is unremarkable ({n:align.vasco_p}% of random catalogues hold as many). With pairs split, chance gives {n:align.vasco_alone}, and the real {n:align.vasco_real} is {n:align.vasco_alone_phrase}. The track score expects {n:align.exp3} and {n:align.alone_exp3} under the two choices, so the answer is the same in either case.

This is why we set the track score’s close-pair distance to 20 arcseconds: it is the shortest distance at which the chance count no longer depends on how close pairs are treated. A shorter distance would raise ε (e.g. {n:align.eff10}% at 10 arcseconds), but the chance count would then depend on how close pairs are treated. Such a gain in efficiency is not worth an unreliable background-only model.

### align-vasco-caption
Left: comparing VASCO’s score and the track score for the {n:align.both} track candidates found by both methods in the real data. Darker cells hold more candidates, and the dashed line represents equal scores across both methods.

Right: the proportion of artificial tracks found by each method, as a function of spatial direction. The track score finds the same share in every direction. VASCO’s score finds more diagonal tracks (because it counts both transients of a close pair), but almost none along east–west and north–south.

### align-excess-title
Looking for an excess of tracks

### align-excess-text
With a reliable background-only model in hand, we compare it with the track candidates in real data that pass our selection. Genuine glints would appear as candidates scoring better than chance allows, above what the random catalogues predict. We fit the real track score distributions for candidates of 3, 4, 5, and 6 transients, first with the background-only model alone, then with glint tracks added. Each glint track is modelled on our artificial ones, with a set number of transients spread over 2 to 10 arcminutes. A track of six transients also contains tracks of five, four, and three, so the four distributions are fitted together.

The first chart below compares the real track score distributions with the background-only model. The model fits the real catalogue as well as it fits a typical random one: in fits to random catalogues, {n:align.fit_bg_p}% have a worse fit quality than the fit to real data. Adding glint tracks as a signal component does not improve the fit: for every track size, the best-fitting number of glint tracks is {n:align.fit_none}. In plain counts, {n:align.n3_lc} candidates score 1 in 1,000 or better in track score, against {n:align.exp3} expected from the random catalogues.

To set limits, we repeated the fit on random catalogues with known numbers of glint tracks added, allowing for the small share of signal that random catalogues absorb (see {3_circle_symbol}). With a flat prior on the number of tracks, the second chart below shows the posterior probability of each number.

The 95% upper limits on the number of glint tracks across all {n:cat.plates} plates are:

- <em>N</em>(3-transient tracks) ≤ {n:align.lim3_n}
- <em>N</em>(4-transient tracks) ≤ {n:align.lim4_n}
- <em>N</em>(5-transient tracks) ≤ {n:align.lim5_n}
- <em>N</em>(6-transient tracks) ≤ {n:align.lim6_n}

The limit is weakest for three transients, as chance alignments mimic these most easily.

### align-excess-caption
Top: track score distributions for all 3-, 4-, 5-, and 6-track candidates in the catalogue. The real data distribution is shown in red, and the expected random background distribution is shown in grey.

Bottom: residual histograms, showing bin-by-bin comparisons of the real data and expected random background. The closer a value is to zero, the better the agreement between the data and the random background.

### align-limit-caption
For glint tracks of 3, 4, 5, and 6 transients, we set 95% limits on their number in the full catalogue. In all cases, the most likely number of genuine glint tracks is zero.

### align-verdict-title
No excess of alignments

### align-verdict
Using our track score method, we searched the full catalogue for linear glints of 3, 4, 5, and 6 transients. We find no excess above the background level expected from random catalogues.

With 95% probability, the catalogue holds at most {n:align.lim6} 6-transient, {n:align.lim5} 5-transient, and {n:align.lim4} 4-transient glint tracks. To the best of our knowledge, these are the first such limits to be set on linear glints.

### end-h2
Summary

### card-4-question
Are the plates correctly dated?

### card-4-ours
Around {n:dates.early_pct}% of plates have an archive date that is one day early. Asteroid trails, Palomar’s own logs, and the handwritten dates on plates independently corroborate this finding. We provide a correction to the dates in our code repository.

### card-1-question
Do transients avoid Earth’s shadow?

### card-1-ours
Allowing for where each plate points, the shortfall of transients in shadow happens by chance about one time in {n:sh.odds_point}. Plate by plate, transients are as dense in shadow as in sunlight, and a shortfall as large as VASCO’s count implies is excluded at {n:fit.inj_cl}% confidence.

### card-2-question
Do transients cluster around nuclear tests?

### card-2-ours
Nights near tests had more photographs, pointed at parts of the sky where more transients survive the catalogue’s selection. Allowing for this, the excess per photograph is {n:nc2.point_r}, which shifted test calendars match about one time in {n:nc2.odds}. VASCO fixed their time window in advance; had the strongest of {n:nc4.windows} reasonable windows been reported, chance would give a result as strong as theirs {n:nc4.chance}% of the time.

### card-3-question
Do transients form linear patterns?

### card-3-ours
Our track score, which looks in every direction and counts each close group of transients once, finds as many track candidates as chance predicts. We set the first limits on the number of glint tracks in the catalogue.

### end-text
Our corrected plate dates, track score, and background-only model are in our code repository, free for anyone to use. We would welcome the VASCO team’s response to this analysis, and thank them again for making their data and code available.

### footer-about
Built from the VASCO team’s public catalogue, VASCO’s analysis code in Brian Doherty’s VASCO-ML repository, and our own scripts. Asteroid positions come from JPL Horizons, star positions from Gaia, plate times from Palomar’s log (VizieR VI/25), and test times from the US Department of Energy report DOE/NV-209. Every result of our own analysis can be reproduced from our repository, which also holds the table of corrected plate dates; figures quoted from VASCO and from other authors are taken from their papers. Watters et al. (2026, PASA) re-examined the shadow and nuclear-test results using other published lists of transients, and the VASCO team has responded (Villarroel et al. 2026, arXiv:2602.15171). We thank Beatriz Villarroel and her colleagues for publishing their data and code; their openness made this work possible.
