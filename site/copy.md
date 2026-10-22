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

Many transients are likely to be artefacts, such as dust on the glass, specks in the photographic emulsion, or marks made in copying and scanning. The catalogue has not, to our knowledge, been inspected transient by transient. The VASCO team has since trained a machine learning model to assess how likely each transient is to be genuine, but the catalogue was built before this work and includes every transient regardless of its model rating. The released catalogue has a column for these ratings, but {n:cat.ml_default}% of transients share the same default value, so the ratings cannot yet be used and our analysis uses the whole catalogue.

### survey-map-caption
Each tile shown represents one plate, coloured by the number of transients on a logarithmic scale. The map is centred on the north celestial pole and keeps areas true, so equal patches of sky occupy equal space. Plates along the Milky Way hold fewer transients, as do plates close to a patch near Taurus and Perseus. Plates outlined with a dashed line are present in the catalogue but were not matched to a night in VASCO’s list; we include them in our analysis.

### dates-h2
When was each plate taken?

### dates-intro
In the sky map above, {n:cat.left_out} plates are outlined with a dashed line. These are plates that VASCO’s tests left out, because their archive dates do not match a night in VASCO’s list.

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
Asteroids move relative to stars, typically by around half an arcsecond per minute. During long exposures, asteroids thus appear as short trails, while stars appear fixed as the telescope tracks them. Using NASA’s JPL Horizons service, we can determine the position of every known asteroid on any night in the 1950s with high precision.

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

Where the log covers the same plates as the asteroids or the handwritten labels, they agree on {n:dates.log_agree} of {n:dates.log_checked} cases. The two exceptions are typing errors in the log, where we take the date on the label. For every other plate in a run with a one-day error, we add one day to the archive date.

Across the {n:dates.dated} plates we dated (which include all {n:cat.plates} in the catalogue), Palomar’s logs date {n:dates.by_log}, asteroid trails {n:dates.by_asteroids}, and the handwritten labels another {n:dates.by_label}. Runs of plate numbers settle a further {n:dates.by_runs}, {n:dates.runs_cat}. That leaves {n:dates.unresolved} plates, XE723 and XE884, with no evidence for their date; neither is in the catalogue, so they play no part in our subsequent analysis.

{n:dates.exceptions} plates don't follow this pattern, and are listed in the second table below. {n:dates.gross} of these plates are wrong by months or years. Most are corrected from their labels or from Palomar’s log. The remaining two carry the same error as the plates numbered next to them, so they take those plates’ corrected dates. All of the date corrections can be found in our code repository.

### dates-why
Most of the outlined plates on the sky map above are a side effect of the dates. VASCO’s tests use a plate only if its date is on their list of {n:dates.vasco_nights} observing nights, taken from Palomar’s log. The log names each night by the evening it began: the night of 8–9 April is listed as 8 April.

A plate whose archive date is one day early therefore shows the evening date, and matches the list. A plate whose archive date is correct shows the morning date, 9 April, and does not match its own night:
- If Palomar did not observe on the night of 9 April, the plate matches no night and is left out. This accounts for {n:dates.dropped} of the {n:cat.left_out} outlined plates.
- If Palomar did observe that night, the plate is kept but placed on the wrong night. This happened to {n:dates.shifted} plates.

With corrected dates, every plate matches the night it was taken, and our analysis uses all {n:cat.plates}.

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

A plate’s edges are set by its size and shape, not by anything in the sky. Wherever the telescope points, real sources have no reason to gather towards the edges of the plate. The rise may come from the plates themselves, and how the catalogue was built. A detection counts as a transient only if nothing in the Gaia, Pan-STARRS, or WISE catalogues lies within 5 arcseconds of it (Solano et al. 2022). Towards the edges of a plate, star images blur and their measured positions are less accurate, so real stars may miss that match and enter the catalogue as transients. Watters et al. (2026) find the same rise towards plate edges and corners in other lists of transients, and find that transients clump together on most plates.

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

One further effect is tied to particular plates rather than to the sky. The {n:cat.eq_plates} plates centred on the celestial equator, the southern edge of the catalogue, hold about {n:cat.eq_share}% as many transients as other rows of plates once the factors above are accounted for. Where these plates overlap their northern neighbours, the same patch of sky yields only about {n:cat.eq_same_sky}% as many transients on the equatorial plates. The cause therefore lies in the plates themselves or how they were processed, not in the sky. The shortfall is already present in the list published by Solano et al. (2022), and is not explained by any of the recorded properties of the plates.

### cat-pointing-fit
To assign each plate an expected number of transients based on pointing alone, we combine the above factors in a single fit across all plates. In this model, dust is averaged across each plate, and we include a term for the equatorial row. In the following sections on Earth’s shadow and nuclear tests, we use these expected counts to compare plates on equal terms. Where appropriate, the model is refitted without the plates being tested, so that they cannot shape their own expected counts.

### cat-pointing-caption
Averaged over groups of similar plates, the pointing model predicts transient counts well. We find that pointing accounts for {n:cat.pointing_share}% of the differences between plates.

To avoid this pointing-based variation confounding measurements, we use our pointing model to remove this component in the Earth's shadow and nuclear tests to follow.

### cat-rows-caption
For plates in the equatorial row, we find that the level of transients is significantly lower than expected given their pointing. We thus include this effect as an additional factor in our pointing model.

### shadow-h2
Do transients avoid Earth’s shadow?

### shadow-intro
An object in orbit is visible only via the sunlight it reflects. When orbiting objects pass into Earth’s shadow, they should therefore go dark.

VASCO considered objects in geostationary orbit (36,000 kilometres up), where a single orbit takes one day and a satellite stays above a point on the equator. At this height, Earth’s shadow covers a disc of sky around 17 degrees across, centred on a point directly opposite the Sun. If some transients are due to glints from objects in orbit, the shadow disc should exhibit fewer transients than the same area of sky elsewhere.

Flaws in the photographic emulsion, dust on the glass, and marks from copying have no connection with where Earth’s shadow fell. If fewer transients fall inside the shadow than the differences between plates can explain, that would point to something in orbit. Overall, we find the same shortfall of transients in Earth's shadow as VASCO. Where we differ is on the statistical interpretation of the deficit.

### shadow-step1-title
Plates intersecting Earth's shadow

### shadow-step1-text
Across the {n:cat.plates} plates in the catalogue, Earth's shadow happened to intersect with {n:sh.plates_n} of them. {n:sh.vasco_plates_cap} of these are among the plates in VASCO's tests, as Brian Doherty, a co-author of the 2026 VASCO paper, also found. The remaining plate, XE472, was taken on 18–19 November 1949, two nights before VASCO's list of nights begins.

Throughout this section, we use the central 6 by 6 degrees of each plate, matching the plate size VASCO use in their own analysis. VASCO’s random positions cover this square evenly on every plate, as shown above.

Using our corrected dates, the areas of those plates falling in Earth's shadow are found to hold a total of {n:sh.observed} transients.

Together, the shadowed areas cover {n:sh.worth} plates’ worth of sky (all of plate XE422, two-thirds of two others, nearly half of XE472, and smaller slices of the remaining six). An average plate holds {n:sh.avg} transients in its central 6 by 6 degrees, so if transients were spread evenly across the survey, the shadowed parts would be expected to hold about {n:sh.expected} transients ({n:sh.worth} × {n:sh.avg}).

VASCO’s archive dates place the shadow slightly differently. A date one day out moves the shadow by only about a degree, small compared to its 17-degree width, so on the {n:sh.vasco_plates} plates in their tests the shortfall is similar: {n:sh.vasco_obs} transients against {n:sh.vasco_exp} expected, nearly a third fewer.

### shadow-step1-caption
Of the {n:cat.plates} plates in the catalogue, {n:sh.plates_n} intersect with Earth's shadow. This includes one plate (XE422) which falls fully within Earth's shadow.

### shadow-step2-title
Measuring the probability of this deficit

### shadow-step2-text
VASCO scores the shortfall with a binomial test, which treats every transient as an independent event. Their primary figure of 22 standard deviations compared the shadow with its share of the whole sky; their 7.6 standard deviation result used each plate’s actual coverage, as we do here. With that model and our corrected dates, the expected count of {n:sh.expected} would vary by about its square root ({n:sh.binom_sd}), which puts the observed {n:sh.observed} about {n:sh.binom_z} standard deviations below expectation.

Transients are not independent, however: as shown above, transient counts on plates differ considerably more than counting noise permits. The count inside a given shadow therefore depends on whether the plate it falls on happens to be rich or poor in transients. To measure how far the count can move, we took the {n:sh.plates} Earth shadow patches, each with its real shape and position, and laid each of them onto a plate chosen at random.

A randomly chosen plate has no connection to Earth's shadow, so any shortfall in such trials comes from chance alone. Plates also differ in where they point, so each trial count is scaled by the pointing model to the shadow plate’s pointing (with the model refitted without the shadow plates). Across 20,000 such trials, a shadow transient count of {n:sh.observed} or fewer turned up {n:sh.p_point}% of the time, a statistical significance of {n:sh.sig_point} standard deviations. Without the allowance for pointing, the figure is {n:sh.p_random}%.

Where the shadow fell on plates also matters. Earth’s shadow falls near the point opposite the Sun, which lies on the ecliptic, so it tends to cross plates in the thinnest parts of the catalogue. Of the {n:sh.plates} shadow plates, {n:sh.taurus_lc} lie towards Taurus and Perseus, where the catalogue holds unusually few transients, and {n:sh.eq_plates} are on the equatorial row (see the section on transients per plate).

The pointing curve allows for the equatorial row and for dust across each plate. Even so, {n:sh.empty_n} of the shadow plates, {n:sh.empty_plates}, hold under a tenth of what their pointing predicts, and none of our modelled factors explain why. Wesley Watters and colleagues made the same point, and found with their own random-shadow test, on a larger published list of transients that overlaps VASCO’s, that a deficit this large arises by chance about 30% of the time. The shadow plates were all taken with their centres {n:sh.shadow_lo} to {n:sh.shadow_hi} degrees from the point opposite the Sun. Drawing the random plates only from the {n:sh.near_hosts} plates within {n:sh.near_max} degrees of that point, with each count scaled for pointing as before, gives a probability of about one in {n:sh.odds_near} ({n:sh.sig_near} standard deviations).

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

Earth’s shadow splits each plate into two section. On XE617, the shadow covers {n:pl.XE617_area}% of the plate, and holds {n:pl.XE617_O} of its {n:pl.XE617_n} transients. On XE640, the shadow clips {n:pl.XE640_area}% of the plate, in a corner that holds none of its {n:pl.XE640_n} transients.

### shadow-plates-caption
The {n:sh.partial} plates that Earth’s shadow partly covers, showing the central 6 by 6 degrees of each. Each dot is a transient in VASCO’s catalogue. Transients in Earth's shadow are shown in red, and those in sunlight blue. The red line marks the edge of Earth's shadow, 8.5 degrees from the point opposite the Sun when the plate was taken.

Under each plate, “expected” shows how many transients the shadowed region would hold if transients there were as dense as in the sunlit region, allowing for the higher density towards plate edges (measured on plates the shadow never crosses).

The ratio is the observed number of shadow transients divided by expected. A ratio of 1 means no difference in transient density between the shadow and sunlit regions.

### shadow-fit-title
Fitting the shadow transient density

### shadow-fit-text
VASCO report that transients are rarer in Earth’s shadow. Expressed as a measurement: in the shadowed part of a plate, transients are r times as dense as in the sunlit part. With no Earth's shadow effect, the ratio r should be 1. If some transients genuinely vanish in the shadow, r should be below 1.

On each plate, r compares the share of its transients that fall in the shadowed part with the share that part would hold if the shadow made no difference. That expected share is not simply the shadowed fraction of the area, because transients are denser towards the edges and corners of every plate (see the plate maps above). We therefore weight each part of the plate by the typical density of transients there, measured on the {n:sh.never} plates the shadow never crosses, so that a real shortfall in shadow cannot affect the weights.

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

Over {n:fit.trials} such trials, the fit gives a typical r value of {n:fit.null_median}, with the middle 68% of results falling between {n:fit.null_lo} and {n:fit.null_hi}. Our best fit result on the {n:sh.partial} genuine plates (r = {n:fit.r}) falls well within this range, indicating that there is no evidence of reduced transient signals in Earth's shadow beyond random chance.

The same trials show that our plate-based fit can detect a genuine shortfall. When we remove transients from inside the shadow outlines on the random plates, so that they are only {n:fit.inj_r} times as dense there (a shortfall as large as VASCO’s), the fit returns a value as high as {n:fit.r} only {n:fit.inj_pct}% of the time. VASCO did not measure the effect plate by plate, but a shortfall of the size their count implies is excluded at {n:fit.inj_cl}% confidence (one-sided).

These results hold when the random plates are selected differently, when each shadow is instead moved around on its own plate, or when any one plate is left out.

### shadow-toys-caption
Blue: fit result to the real plates.

Grey: {n:fit.trials} trials using the real shadows applied to random plates.

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
During the sky survey, the United States tested nuclear weapons in Nevada and the Pacific. In two papers (Scientific Reports, 2025 and 2026), VASCO report that observing nights within a day of a nuclear test held more transients than other nights. The 2026 paper uses US tests and weights each transient by a machine-learning rating of how likely it is to be real; it puts the chance of the excess at about 1 in 40 (p = 0.024). It keeps the window of one day either side from the 2025 paper, as a hypothesis set before the analysis, and the authors read the sharpness of the association in time as a sign that it is genuine.

Repeating their measurement with the full list of US tests (VASCO’s list leaves out the {n:nuc.missing} tests of spring 1953), we find a similar excess. It reduces to a level consistent with chance once two further factors are taken into account:

- Nights near tests tended to have more photographs, and each photograph adds its own transients.
- On nights close to tests, the telescope tended to point at parts of the sky where more transients survive the catalogue’s selection, mostly away from the Milky Way. Our pointing model (see the section on transients per plate) expects these plates to hold about {n:nc2.exp_pct}% more transients than others.

### nuc-vasco-title
Repeating VASCO’s test

### nuc-vasco-text
VASCO’s test counts the transients on each observing night, compares the nights within a day of a US test with all other nights, and asks how often the same number of nights selected at random show at least as many transients. We count each transient once, since the released catalogue does not yet hold usable machine-learning ratings.

On VASCO’s own construction, using the default archive dates and their list of test nights, the {n:nuc.vasco_nights} nights near tests hold {n:nuc.vasco_ratio} times the usual number of transients. Random nights exhibit at least as many transients around one time in {n:nuc.vasco_odds} ({n:nuc.vasco_sig} standard deviations), close to VASCO’s 1 in 40 result.

Our version of this measurement applies two corrections. We use the corrected dates described above, and the full list of US tests from the Department of Energy (DOE/NV-209). VASCO’s 2026 paper cites the same report, but the list in their code repository, which gives the paper’s {n:nuc.vasco_nights} nights near tests, leaves out all {n:nuc.missing} tests of Operation Upshot–Knothole, in spring 1953.

As in VASCO’s test, a night counts based on the local date on which it began. With both corrections, {n:nuc.nights} nights fall within a day of a test. They hold {n:nuc.ratio} times the usual number of transients, with random nights doing as well about one time in {n:nuc.odds} ({n:nuc.sig} standard deviations).

### nuc-vasco-caption
The excess of transients on the {n:nuc.nights} nights within a day of a US test (red line), compared to the same number of nights drawn at random {n:nuc.trials} times from the survey (grey). The analysis is performed using corrected dates and the full list of US tests.

### nuc-corr-title
Counting photographs, and considering telescope pointing

### nuc-corr-text
VASCO’s test counts transients per night, but a night’s transients come from its photographs. The {n:nuc.nights} nights near nuclear tests averaged {n:nuc.plates_window} plates per night, compared to {n:nuc.plates_other} on other nights. Part of this gap comes from nights with no plates in the catalogue: {n:nc2.empty_near}% of the nights near tests, against {n:nc2.empty_other}% of other nights. A night with no plates contributes zero transients, even though nothing was observed. Watters et al. (2026) make the same point: counted over only the days with plates north of the celestial equator, the link between transients and test days is no longer significant. Counted per plate rather than per night, the excess reduces from {n:nc2.night_r} to {n:nc2.photo_r}.

The sky above Palomar changes through the year as Earth orbits the Sun. On spring nights, the telescope pointed away from the Milky Way. There, fewer stars crowd the field, so fewer detections are rejected for lying near a known star, and more survive as transients. Plates within 5 degrees of the Milky Way hold a median of {n:crowd.inner} transients, compared to {n:crowd.outer} in plates more than 40 degrees away. Of the plates near nuclear tests, {n:nc2.spring}% were taken between February and June. Allowing for where every plate pointed, our pointing model expects the plates near tests to hold around {n:nc2.exp_pct}% more transients than others.

A like-for-like comparison uses the pointing curve detailed earlier: the number of transients each plate’s pointing alone would leave in the catalogue. Measured against expected counts, plates near nuclear tests hold {n:nc2.point_r} times their expected number, rather than {n:nc2.photo_r}. Fitting the curve without the plates near tests, so they cannot shape their own expected counts, also gives {n:nc2.held_out}.

To judge how often chance alone gives a result like {n:nc2.point_r}, VASCO’s test draws the same number of nights at random from across the survey. Random nights do at least as well {n:nc2.point_pv}% of the time ({n:nc2.point_pv_sig} standard deviations).

Nights near tests are not a random set, however: the tests came in series, so those nights bunch into a few observing runs, which share batches of emulsion and similar conditions. As a stricter check that keeps this bunching, we slide the real test calendar to other years, leave out plates taken within two weeks of a real test, and score each shifted calendar in the same way. Shifted calendars match or exceed the result {n:nc2.point_pc}% of the time ({n:nc2.sig} standard deviations).

It is the allowance for pointing that removes the excess: counted per photograph without it, shifted calendars match or exceed the excess only {n:nc2.photo_pc}% of the time. Importantly, our pointing correction does not remove a genuine excess. To check this, we add a fake excess, multiplying the number of transients on every plate near a test by {n:nc2.inj_f}, and repeat the whole measurement. The pointing curve is refitted to all plates, fake excess included, so it has every chance to absorb the excess. It does not: the plates near tests still hold {n:nc2.inj_r} times their expected number, so the fake excess is recovered.

Using VASCO’s own US test list instead of the full list gives {n:nuc.vlist_r}, which shifted calendars match or exceed {n:nuc.vlist_pct}% of the time.

### nuc-corr-caption
Once the number of photographs and telescope pointing are accounted for, the plates near nuclear tests hold {n:nc2.point_r} times their expected number of transients.

Compared with {n:nc2.slides} shifted test calendars, each catching about as many plates as the real one, about {n:nc2.frac} show at least the same level of excess.

### nuclear-verdict-title
No significant link to the nuclear tests

### nuclear-verdict
We reproduce VASCO’s excess of transients near nuclear tests, but find that it is primarily due to where the telescope pointed. Once pointing is allowed for, no significant effect of the nuclear test schedule on transient counts is found.

### align-h2
Do transients form geometric patterns?

### align-intro
A flat object tumbling in orbit would flash each time it caught the Sun, leaving a line of points during a long plate exposure. VASCO searched for such alignments and shortlisted five candidates. Their best candidate was reported at 3.9 standard deviations, corresponding to an estimated 1 in 10,000 probability of arising by chance.

A catalogue of 100,000 unevenly spread points may contain many chance alignments. Three questions naturally arise:

- Do the alignment quality checks used introduce biases in the directions of the alignments they can find?
- Does the real transient catalogue contain more alignments than expected by random chance?
- How many genuine glint tracks could the catalogue hold?

### align-step1-title
Searching for alignments in all directions

### align-step1-text
VASCO kept a set of transients only if its two sky coordinates, right ascension (east–west) and declination (north–south), rise and fall together almost perfectly: a correlation coefficient (Pearson’s r) above 0.99. This approach works well for diagonal lines. However, a line running close to east–west barely changes in declination along its length. The small scatter in each measured position then hides that change, so the correlation falls below the cut and the line is rejected. The same happens to lines close to north–south.

In simulations (image below), the VASCO method keeps {n:align.pear_diag}% of diagonal lines, {n:align.pear_2}% of lines 2 degrees from east–west, and only {n:align.pear_ew}% of lines running exactly east–west. A satellite at geostationary height stays nearly fixed above the observatory while the telescope follows the stars, so on the plate it moves close to east–west. Seen from Palomar, the geostationary belt itself lies just south of the catalogue’s plates, so any such object on them is on a slightly tilted orbit, and its track runs within a few degrees of east–west. VASCO’s selection would reject almost all of the tracks such an object would leave.

The same alignment task arises in particle physics, where the paths of charged particles must be reconstructed from hits they leave in a detector. The simplest tracking method fits a straight line through a set of points, resulting in a track candidate. The quality of each candidate can be assessed by checking how far the points lie from the line of best fit. Such a method works equally well in all directions, including east–west when applied to this setting.

### align-step1-caption
The share of simulated straight lines kept by each method, shown as a function of direction in the sky. VASCO’s selection (orange) collapses along east–west, where genuine geostationary satellite glints would lie, as well as along north–south. A straight-line fit (blue) keeps nearly all lines in every direction.

### align-step2-title
Fitting every track candidate

### align-step2-text
VASCO’s five shortlisted alignments came from an earlier list of 298,165 candidate transients (Solano et al. 2022). Of their {n:align.vasco_transients} transients, {n:align.vasco_in_cat} appears in the later catalogue of {n:cat.transients} transients that VASCO uses for its other tests. We have applied our track-fit method to the full catalogue of {n:cat.transients} transients.

We build track candidates from pairs of transients up to 15 arcminutes apart. Each pair marks the two ends of a track, and the transients lying closest to the line between them are added, giving three to six transients in total. We then fit a straight line through each candidate, and measure how closely its transients follow it.

Two sources closer than about 4 arcseconds (the width of a star’s image on the scans) blur into one, so we ignore candidates with points closer than 8 arcseconds. The released catalogue itself holds almost no transients closer together than this. {n:align.vasco_close_n} of VASCO’s five shortlisted alignments contain points only {n:align.vasco_close} arcseconds apart, closer together than almost any pair in the released catalogue. Watters et al. (2026) also find that 8 of the 24 transients in the five alignments lie within 5 arcseconds of sources seen by the NEOWISE infrared survey, so may be stars rather than transients.

The position of each transient on a scan is uncertain by about an arcsecond. As such, even points from a genuine straight track can stray from a true straight line. Each track fit is scored with a χ² per degree of freedom (χ²/ndof), which compares the points’ actual distances from the line with a one-arcsecond uncertainty. A genuine track gives a value near 1 or below. Points scattered loosely about the line give a much larger value.

A good fit alone does not make a track stand out: in crowded parts of a plate, pure chance can cause transients to align frequently. We therefore assign each track candidate a track score, which measures how easily chance alone could produce it. The score uses four numbers, as explained in the diagram below: the candidate’s number of transients (r), its length (L), its width (w), and the number of transients in a box at least 10 arcminutes across around it (N).

This concept, and the formula we use, aligns with the method VASCO used to judge its own candidates (Edmunds and George, 1985). We apply our method in all directions, measure straightness from the fitted line, and score only track candidates that are straight to within the measurement error (a χ² probability of at least 1%). The lower the track score, the less often chance alone can produce such a track.

Transients often come in close pairs: the catalogue holds far more pairs of transients within 20 arcseconds of each other than chance would give. Two transients this close fix a direction only loosely, so a third lying roughly in line with them makes a short track that looks stronger than it is. Following Lezama and colleagues (2015), who developed a method for finding alignments among clustered points, we count a track’s transients in boxes 20 arcseconds long, the distance over which close pairs are unusually common: transients closer together than this along the track count once.

Candidates built from a set of six transients also yield candidates composed of any combination of three, four, or five of that set. Where track candidates share a transient, we keep only the best-scoring candidate, such that each transient counts once.

### align-method-caption
How the track score rates a real track candidate on plate {n:align.ex_plate}. A candidate is scored only if it first passes the straightness test.

The track score, μ, is the number of tracks with r transients, at most L long and lying within w of a straight line, that we would expect if the N transients in the box were scattered at random, with n their density (Edmunds and George 1985):

{eq:track-score}

### align-step3-title
Measuring how often tracks arise by chance

### align-step3-text
Any catalogue of transients will hold tracks by chance, so finding track candidates proves nothing by itself. The proper test is whether the real data holds more strong track candidates than chance alone can produce.

This test requires a robust background-only model: random catalogues which match the real data in every respect, except that they contain no genuine tracks. Since transients cluster on plates, and their surface density varies, a uniform random distribution cannot serve as this model.

We therefore created {n:align.copies} random catalogues based on the real data itself. Each random catalogue keeps all of the plates from the real data, including the number of transients on each and their tight groups of transients, while breaking any alignments. Every lone transient is moved at random by up to {n:align.shift} arcseconds; this is enough to break genuine lines, but small enough that each plate keeps its distribution of transients on larger scales.

Transients within {n:align.link} arcseconds of one another, either directly or through a chain of neighbours, form a group. Each group is moved as a single piece, turned at random and shifted by up to {n:align.shift} arcseconds, and no transient is moved to within {n:align.link} arcseconds of another group. Because a group keeps its shape when it moves, a line formed within a group would also appear in every random catalogue, and the comparison could not tell whether it arose by chance. We therefore leave out any track holding three or more transients from one group, in both the real and random catalogues.

To measure chance correctly, the random catalogues must satisfy two conditions:

- <strong>They must break up real tracks.</strong> To test this, we added 100 artificial tracks to the real catalogue and created random catalogues from it in the same way. For every 100 added tracks scoring 1 in 10,000 or better, the random catalogues gain {n:align.absorb4}; at 1 in 1,000 or better, they gain {n:align.absorb3}. When comparing the real data to random catalogues, any genuine excess of glint tracks would therefore still show at {n:align.keep4}% of its true size ({n:align.keep3}% for 1 in 1,000). The limits we set on glint tracks below allow for this, track by track.
- <strong>They must not understate chance.</strong> If the random catalogues held fewer chance tracks than the real data, chance alone could pass for a signal. In practice, the random catalogues err slightly the other way. Pairs of transients 20 to 30 arcseconds apart are rare in the real catalogue, and moving the transients fills that gap, giving the random catalogues {n:align.pair2030} times as many such pairs and thus slightly more chance tracks. This makes the background-only test conservative: it could slightly hide a real excess, but it cannot create a false one.

The real dataset and random catalogues are processed using the same analysis code, to ensure track candidates are found in an identical fashion.

### align-vasco-title
Comparing our track score and VASCO’s score

### align-vasco-text
VASCO’s method uses the same formula as our track score ({2_circle_symbol}), applied to candidates selected by a correlation screen rather than a straight-line fit, and counting both transients of a close pair. VASCO interpret each candidate’s score directly as its chance of arising, which assumes transients are scattered uniformly. We judge scores against random catalogues built from the real data instead, and we ran VASCO’s score through the same catalogues. Where both methods find a track candidate, they agree closely (rank correlation of {n:align.rho}, with {n:align.top50} of the top 50 candidates in common).

However, the track score rejects {n:align.vasco_bent} of the {n:align.vasco_strong} candidates that VASCO’s score rates as 1 in 1,000 or better. This is because their component transients do not fall along a straight line within the measurement error, and thus they are rejected by our χ² probability requirement. The track score method also finds tracks in all directions, avoiding the spatial biases shown in {1_circle_symbol}.

We measured each method’s efficiency as the share of 100 artificial glint tracks it finds, at a pass mark that lets through around {n:align.false_alarms} chance candidates across an entire random catalogue. The track score achieves an efficiency of {n:align.track_ew}% in every direction. VASCO’s method has a {n:align.vasco_diag}% efficiency for diagonal tracks, dropping to {n:align.vasco_2}% two degrees from east–west. VASCO’s method does better on diagonal tracks because it counts both transients of a close pair (under 20 arcseconds apart), whereas the track score counts only one.

This higher efficiency comes at a cost. A close pair could be a feature of the plate, such as a defect detected twice. That is part of the background, which the random catalogues must reproduce, so the pair should move as one piece. Or it could be two glints from one object. That is part of the signal, which the random catalogues should break up, so the two should move independently. For any individual close pair, it is not possible to know which of these is correct.

Using VASCO’s score, the two possibilities lead to opposite analysis outcomes. With close pairs kept together in the random catalogues, chance gives {n:align.vasco_exp} tracks scoring 1 in 1,000 or better. When measured against this, a count of {n:align.vasco_real} candidates in the real data is not significant ({n:align.vasco_p}% of random catalogues hold as many). When splitting close pairs in the random catalogues, only {n:align.vasco_alone} tracks scoring 1 in 1,000 or better are found; in this case, none of the random catalogues holds as many as the real data. Whether VASCO’s {n:align.vasco_real} candidates are significant thus depends on an unknown property of the data.

By contrast, the track score method expects {n:align.exp3} and {n:align.alone_exp3} candidates under the two choices, demonstrating that its background-only model is robust to this unknown aspect of the data. Twenty arcseconds is the shortest close-pair distance for which this equivalence holds: at shorter distances, a higher track-finding efficiency is possible (e.g. {n:align.eff10}% at 10 arcseconds), but at the cost of a background-only model which depends on that same unknown.

### align-vasco-caption
Left: comparing VASCO’s score and the track score for the {n:align.both} track candidates found by both methods in the real data. Darker cells hold more candidates, and the dashed line represents equal scores across both methods.

Right: the proportion of artificial tracks found by each method, as a function of spatial direction. The track score finds the same share in every direction. VASCO’s score finds more diagonal tracks (because it counts both transients of a close pair), but a much lower proportion along east–west and north–south.

### align-excess-title
Looking for an excess of tracks

### align-excess-text
With a reliable background-only model in hand, we compare it with the candidates in real data that pass our track score selection. Genuine glints would appear as candidates scoring better than chance allows, showing as an excess beyond what the random catalogues predict. We fit the real track score distributions for candidates of 3, 4, 5, and 6 transients, first with the background-only model, then with glint tracks added. Each glint track is modelled on our artificial ones, with a set number of transients spread over 2 to 10 arcminutes. A track of six transients also contains tracks of three, four, and five, so the four distributions are fitted together.

The chart below compares the real track score distributions with our background-only model. The model fits the real data as well as it fits a typical random catalogue: in fits to random catalogues, {n:align.fit_bg_p}% have a worse fit quality than the fit to real data. Adding glint tracks as a signal component does not improve the fit quality: for every track size, the best-fitting number of glint tracks is {n:align.fit_none}. In plain counts, {n:align.n3_lc} candidates score 1 in 1,000 or better in track score, against {n:align.exp3} expected from the random catalogues.

### align-excess-caption
Top: track score distributions for all 3-, 4-, 5-, and 6-track candidates in the catalogue. The real data distribution is shown in red, and the expected random background distribution is shown in grey.

Bottom: residual histograms, showing bin-by-bin comparisons of the real data and expected random background. The closer a value is to zero, the better the agreement between the data and the random background.

### align-limit-text
We set limits by asking how many glint tracks the catalogue could contain while still yielding a best fit of zero. For each possible number of tracks, we add that many artificial glint tracks to random catalogues, many times over, and fit each one exactly as we fit the real data.

As the number of added tracks grows, a best fit of zero becomes rarer. Its frequency at each number of tracks gives the likelihood of the observed result; with a flat prior, this yields the posterior probability for each number of tracks, shown in the chart below.

The 95% upper limits on the number of glint tracks across all {n:cat.plates} plates are:

- <em>N</em>(3-transient tracks) ≤ {n:align.lim3_n}
- <em>N</em>(4-transient tracks) ≤ {n:align.lim4_n}
- <em>N</em>(5-transient tracks) ≤ {n:align.lim5_n}
- <em>N</em>(6-transient tracks) ≤ {n:align.lim6_n}

The limit is weakest for three-transient tracks, which chance alignments mimic most easily.

The limits rest on two assumptions about how glints appear on the plates. Firstly, they apply to tracks like our artificial ones, with flashes spread over 2 to 10 arcminutes. Flashes within 20 arcseconds of one another count once, so a glint whose flashes all fall this close together would not be detected, and our test cannot determine whether the catalogue’s close pairs are themselves glints. Secondly, they assume each flash’s position is measured to within 1 arcsecond. At 2 arcseconds, tracks are more difficult to find, and the limits rise to {n:align.lim3_2}, {n:align.lim4_2}, {n:align.lim5_2}, and {n:align.lim6_2} for three to six transients.

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
Allowing for where each plate points, the shortfall of transients in shadow happens by chance about one time in {n:sh.odds_point}. Plate by plate, transients are as dense in shadow as in sunlight, and a shortfall as large as VASCO’s count implies is excluded at {n:fit.inj_cl}% confidence (one-sided).

### card-2-question
Do transients cluster around nuclear tests?

### card-2-ours
Nights near tests had more photographs, and pointed at parts of the sky where more transients survive the catalogue’s selection. Allowing for this, and using the full list of US tests, the excess per photograph is {n:nc2.point_r}, which shifted test calendars match or exceed {n:nc2.point_pc}% of the time.

### card-3-question
Do transients form linear patterns?

### card-3-ours
Our track score, which looks in every direction and counts each close group of transients once, finds as many track candidates in the released catalogue as chance predicts. We set the first limits on the number of glint tracks it holds; VASCO’s five shortlisted alignments come from an earlier list and are not tested here.

### end-text
Our corrected plate dates, track score, and background-only model are in our code repository, free for anyone to use. We would welcome the VASCO team’s response to this analysis, and thank them again for making their data and code available.

### footer-about
The analysis is our own code, run on the VASCO team’s public data from Brian Doherty’s VASCO-ML and VASCO-challenge repositories: their transient catalogue, random positions on each plate, and lists of observing nights and nuclear tests. Where we reproduce VASCO’s results, we follow their published method and code. Asteroid positions come from JPL Horizons, star positions from Gaia, plate times from Palomar’s log (VizieR VI/25), and test times from the US Department of Energy report DOE/NV-209. Every result of our own analysis can be reproduced from our repository, which also holds the table of corrected plate dates; figures quoted from VASCO and from other authors are taken from their papers. Watters et al. (2026, PASA) re-examined the shadow and nuclear-test results using other published lists of transients, and the VASCO team has responded (Villarroel et al. 2026, arXiv:2602.15171). Their response argues that Watters et al. used a reduced, filtered subset of transients, lacking complete timing information, too small for the shadow test, and without uncertainties. The analysis here uses VASCO’s full released catalogue, every plate’s corrected date and time, and uncertainties from random trials built on the real data. We thank Beatriz Villarroel and her colleagues for publishing their data and code; their openness made this work possible.
