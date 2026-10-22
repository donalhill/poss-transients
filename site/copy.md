### hero-title
Examining the evidence for glints from orbit in 1950s sky photographs

### hero-lede
Between 1949 and 1957, Palomar Observatory in California photographed the northern sky on glass plates. The plates show faint points of light that have no counterpart in modern sky catalogues. The VASCO project, led by Dr Beatriz Villarroel, has catalogued more than 100,000 of these transients.

Villarroel’s team has proposed that some transients may be produced by sunlight glinting off flat, reflective objects in orbit around Earth in the pre-satellite era. In papers published in 2025 and 2026, the team reported three statistical patterns in support of this theory: fewer transients appear inside Earth’s shadow; more transients appear around the dates of nuclear weapons tests; some transients appear along straight lines.

We are grateful to the team for publishing its catalogue and code. Building the catalogue took years of careful work, and releasing it openly made this re-examination possible. We repeat the shadow and nuclear-test analyses on their data, and apply both their alignment method and our own to the released catalogue.


### hero-author
PhD in particle physics from the University of Oxford and CERN, with over ten years of research experience on the LHCb experiment. Now a data scientist in industry, following VASCO's work with interest.

### survey-h2
The Palomar Sky Survey

### survey-intro
Images were taken on plates using Palomar’s 48-inch Schmidt telescope, which was built to photograph wide areas of sky. Each plate is a sheet of glass 14 inches square, and covers a patch of sky 6.5 degrees across (about thirteen full Moons). Every patch was photographed through both a blue and a red filter, with red exposures lasting 40 to 60 minutes. Copies of the plates were scanned in the 1990s, and astronomers now work from those scans.

VASCO compared the red plates with three modern catalogues of the sky: Gaia, Pan-STARRS, and WISE. A point of light was kept as a transient if, among other checks, no modern catalogue lists anything within 5 arcseconds of it. An arcsecond is a 3,600th of a degree, and 60 arcseconds constitutes an arcminute. The public VASCO catalogue holds {n:cat.transients} such transients spanning {n:cat.plates} plates.

Many transients are likely to be artefacts, such as dust on the glass, specks in the photographic emulsion, or marks made in copying and scanning. The catalogue has not, to our knowledge, been inspected transient by transient. The VASCO team has since trained a machine-learning model that rates each transient from 0 to 1 by how likely it is to be genuine. The ratings rank transients, rather than giving calibrated probabilities. VASCO’s 2026 analyses count each transient by its rating, so that a transient rated 0.8 counts as 0.8 of a transient. We call the total a rated count. The ratings are published in a separate file alongside the team's code, rather than in the catalogue itself. The model learned from a few hundred transients that the team labelled by eye as real or as flaws, and rates each transient from features of its image together with properties of its plate and of the night it was taken.

We count each transient by its rating in both the shadow and nuclear tests, as VASCO do. For the nuclear test, we also show how much the result depends on how the ratings were trained.

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

Astronomers keep time using Universal Time (the clock at Greenwich), which is seven or eight hours ahead of California depending on daylight saving time.

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

Using the archive’s date, we find no trail where the asteroid is predicted to be. One day later, however, we find trails left by asteroids at their expected positions. On the plates we checked, {n:dates.ast_late_trails} asteroid trails across {n:dates.ast_late_plates} plates appear one day after the archive date. On another {n:dates.ast_ok_plates} plates, {n:dates.ast_ok_trails} trails confirm that the archive date is correct, so the method is able to distinguish between cases.

### dates-asteroid-1
Plate XE356, showing asteroid (210) Isabella. The archive date for this plate is 29 November 1951, but the asteroid trail actually appears on 30 November.

### dates-asteroid-2
Plate XE617, showing asteroid (163) Erigone. The archive date for this plate is 8 April 1956, but the asteroid trail appears on 9 April.

### dates-fix-title
Correcting the dates

### dates-fix-text
Palomar numbered the red plates in the order they were taken. The one-day error affects long, unbroken runs of plate numbers, listed in the first table below. Palomar’s own log of the survey also survives. It records each night by its evening date, so the correct date in Universal Time is the day after the log’s date.

For the {n:dates.log_checked} plates also dated by an asteroid trail or a handwritten label, the log agrees on {n:dates.log_agree}. The two exceptions are typing errors in the log, where we take the date on the label.

Every one of the {n:cat.plates} plates in the catalogue is dated from its own record: Palomar’s log dates {n:dates.by_log}, asteroid trails {n:dates.by_asteroids}, and the handwritten labels the other {n:dates.by_label}.

{n:dates.exceptions} of these plates do not follow the pattern of their run, and are listed in the second table below. Of these, {n:dates.gross} are wrong by months or years, and {n:dates.one_day} by a single day. Each is corrected from its label, and all of the date corrections can be found in our code repository.

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

If transients were truly independent of one another and spread uniformly over the sky, plates would differ only by counting noise. An average plate holds {n:cat.avg} transients, and the counting noise on that is around its square root ({n:cat.noise}), so most plates would hold between {n:cat.band_lo} and {n:cat.band_hi} transients.

However, only {n:cat.in_band} of the {n:cat.plates} plates fall within this range. The actual counts spread smoothly, from a handful of transients to several hundred, and on {n:cat.over1000} plates more than 1,000. Some of this variation follows from where each plate points on the sky, as the sky map above shows.

How many transients a plate holds depends on the plate itself, and thus transients cannot be treated as fully independent events. A count drawn from only a few plates, such as the {n:sh.plates} crossed by Earth’s shadow, can therefore vary considerably more than counting noise suggests.

### cat-spread-caption
The blue histogram shows the distribution of transient numbers across the {n:cat.plates} plates. The final bar contains the {n:cat.over1000} plates with more than 1,000 transients.

The orange shaded band shows the range expected if transients were distributed uniformly and independently. In this case, each plate would exhibit a number of transients close to the average ({n:cat.avg}), varying only by counting noise.

### cat-edge-title
More transients at plate edges

### cat-edge-text
Across the middle of a plate, transients are spread evenly. Towards the edges they grow denser, reaching over two times the central density. Random positions that the VASCO team generated as a check show no such rise.

A plate’s edges are determined by its size and shape, not by anything in the sky. Wherever the telescope points, real sources have no reason to gather towards the edges of a plate. The observed rise may come from the plates themselves, and how the catalogue was built. A detection counts as a transient only if nothing in the Gaia, Pan-STARRS, or WISE catalogues lies within 5 arcseconds of it (Solano et al. 2022). Towards the edges of a plate, star images blur and their measured positions are less accurate, so real stars may miss that match and enter the catalogue as transients. Watters et al. (2026) find the same rise towards plate edges and corners in other lists of transients, and find that transients clump together on most plates.

### cat-edge-map-caption
Left: density of transients across the face of a plate, averaged over the {n:cat.plates} plates with each plate counted equally. Density is shown relative to the density near the plate centre (north up, east left).

Right: the random positions VASCO generated as a check, which cover the plate evenly up to its edge.

### cat-edge-caption
Density of transients at each distance from the plate centre, shown relative to the centre density.

Blue: transients on average across all plates, with shading showing how far the curve moves when the set of plates is resampled at random (90% of resamples).

Grey: random positions.

Of the {n:cat.plates} plates, {n:cat.edge_share}% show a higher transient density at the edge of the plate compared to the centre.

### cat-crowd-title
Fewer transients near the Milky Way

### cat-crowd-text
A point of light is kept as a transient only if no modern catalogue lists anything within 5 arcseconds. Near the Milky Way, stars are more crowded, causing a greater proportion of candidates to fail the 5 arcsecond threshold.

A typical plate within 5 degrees of the Milky Way holds {n:crowd.inner} transients. Counts rise steadily away from the galactic plane: beyond 40 degrees, a typical plate holds around {n:crowd.outer} ({n:crowd.ratio} times as many).

### cat-crowd-caption
A smooth curve fit to the relationship between transient counts and distance from the galactic plane. Use the buttons to move plates between their position on the graph, and their position in the sky map.

### cat-pointing-title
Where a plate points shapes how many transients it has

### cat-pointing-text
The number of transients retained on a plate depends on where it pointed. Crowding near the Milky Way and dust both remove detections, while plates near the ecliptic, the band of sky that asteroids follow, gain some from asteroids missing from modern lists. These effects come from how the catalogue was built, and act on transients regardless of their physical origin.

Clouds of dust in our own Galaxy dim everything behind them, and the sparse patch near Taurus and Perseus on the sky map (see above) lies behind some of the closest dust clouds.

To separate dust from crowding, we compared plates more than 20 degrees from the Milky Way, where stars are not crowded. The {n:dust.n_dusty} plates that look through the most dust hold a median of {n:dust.median_dusty} transients, compared to {n:dust.median_clear} held by the {n:dust.n_clear} plates which look through the least dust. Dust covers only part of the Taurus–Perseus patch, so it cannot explain the whole shortfall there; we have not identified what does.

One further effect is tied to particular plates rather than to the sky. The {n:cat.eq_plates} plates centred on the celestial equator, the southern edge of the catalogue, hold about {n:cat.eq_share}% as many transients as other rows of plates once the factors above are accounted for. Where these plates overlap their northern neighbours, the same patch of sky yields only about {n:cat.eq_same_sky}% as many transients on the equatorial plates. The cause therefore lies in the plates themselves or how they were processed, not in the sky. The shortfall is already present in the list published by Solano et al. (2022), and is not explained by any of the recorded properties of the plates.

### cat-pointing-fit
To assign each plate an expected number of transients based on pointing alone, we combine the above factors in a single fit across all plates. In this model, dust is averaged across each plate, and we include a term for the equatorial row. In the following sections on Earth’s shadow and nuclear tests, we use these expected counts to compare plates on equal terms. Where appropriate, the pointing model is refitted without the plates being tested, so that they cannot shape their own expected counts. In those tests, the pointing model is fitted to each plate's rated count.

### cat-pointing-caption
Averaged over groups of similar plates, the pointing model predicts transient counts well. Pointing accounts for {n:cat.pointing_share}% of the variation in plates’ transient counts, measured on a logarithmic scale.

The Earth's shadow and nuclear tests that follow use the pointing model to allow for this variation.

### cat-rows-caption
For plates in the equatorial row, we find that the level of transients is significantly lower than expected given their pointing. We thus include this effect as an additional factor in our pointing model.

### shadow-h2
Do transients avoid Earth’s shadow?

### shadow-intro
An object in orbit is visible only via the sunlight it reflects. When orbiting objects pass into Earth’s shadow, they should therefore go dark.

VASCO considered objects in geostationary orbit (36,000 kilometres up), where a single orbit takes one day and a satellite stays above a point on the equator. At this height, Earth’s shadow covers a disc of sky {n:sh.diam} degrees across, centred on a point directly opposite the Sun. If some transients are due to glints from objects in orbit, the shadow disc should show fewer transients than the same area of sky elsewhere.

Flaws in the photographic emulsion, dust on the glass, and marks from copying have no connection with where Earth’s shadow fell. If fewer transients fall inside the shadow than the differences between plates can explain, that would point to something in orbit. Overall, we find the same shortfall of transients in Earth's shadow as VASCO. We differ on our interpretation of how significant this deficit is.

### shadow-step1-title
Plates intersecting Earth's shadow

### shadow-step1-text
Earth’s shadow crosses {n:sh.plates_n} of the {n:cat.plates} plates in the catalogue. {n:sh.vasco_plates_cap} of these are among the plates in VASCO's tests, as Brian Doherty, a co-author of the 2026 VASCO paper, also found. The remaining plate, XE472, was taken on 18–19 November 1949, two nights before VASCO's list of nights begins.

Throughout this section, we use the central 6 by 6 degrees of each plate, matching the plate size VASCO use in their own analysis. VASCO’s random positions cover this square evenly on every plate, as shown above. The shadow's radius is {n:sh.deg} degrees, as in VASCO's 2025 paper.

Using our corrected dates, the areas of those plates falling in Earth's shadow hold {n:sh.transients} transients. Together they cover {n:sh.worth} plates’ worth of sky: all of plate XE422, 70% of two others, half of XE472, and smaller slices of the remaining six. An average plate holds {n:sh.count_avg} transients in its central 6 by 6 degrees, so if transients were spread evenly across the survey, the shadowed parts should hold about {n:sh.count_exp}.

VASCO’s archive dates place the shadow slightly differently. A date one day out moves the shadow by about a degree: this is small compared to the width of the shadow, but still changes which parts of which plates fall inside it. On the {n:sh.vasco_plates} plates in VASCO’s tests, the shadow covers {n:sh.vasco_worth} plates’ worth of sky, on plates averaging {n:sh.vasco_avg} transients. Counting each transient once, as their 2025 paper does, the shadow holds {n:sh.vasco_obs} transients against {n:sh.vasco_exp} expected, nearly a third fewer. VASCO’s binomial model, described in {2_circle_symbol}, puts this at {n:sh.vasco_z} standard deviations, close to their published 7.6.

Counted by rating, the {n:sh.transients} transients in shadow on our corrected dates give a rated count of {n:sh.observed}, against {n:sh.expected} for the same area of an average plate. From here on, counts in this section are rated counts unless specified otherwise.

### shadow-step1-caption
Of the {n:cat.plates} plates in the catalogue, {n:sh.plates_n} intersect with Earth's shadow. This includes one plate (XE422) which falls fully within Earth's shadow.

### shadow-step2-title
Measuring the probability of this deficit

### shadow-step2-text
VASCO score the shortfall with a binomial test, which treats every transient as an independent event. Their headline figure, 22 standard deviations, compared the count in shadow with the shadow’s share of the northern sky. Their figure of 7.6 standard deviations used the area of each plate the shadow actually covers, as we do here. Applied to our corrected dates and rated counts, the binomial model expects {n:sh.expected} transients in shadow, give or take about {n:sh.binom_sd}. The spread of {n:sh.binom_sd} is one standard deviation, and the observed count of {n:sh.observed} lies {n:sh.binom_z} standard deviations below {n:sh.expected}.

Transients are not independent, however: as shown above, transient counts on plates differ considerably more than counting noise permits. The count inside a given shadow therefore depends on whether the plate it falls on happens to be rich or poor in transients. To measure how much the count varies by chance, we laid each of the {n:sh.plates} shadow patches, with its real shape and position, onto a plate chosen at random. We call each such set of placements a trial.

A randomly chosen plate has no connection to Earth's shadow, so any shortfall in a trial comes from chance alone. Plates also differ in where they point, so we scale each trial count to what the pointing model expects at the shadow plate’s pointing. Here the pointing model is fitted to rated counts, leaving out the shadow plates. Across {n:sh.trials} trials, a count of {n:sh.observed} or fewer turned up {n:sh.p_point}% of the time, a statistical significance of {n:sh.sig_point} standard deviations. Without the pointing model, such counts turn up {n:sh.p_random}% of the time. With the {n:sh.ml_deg}-degree shadow radius used in the code for VASCO's 2026 paper, they turn up {n:sh.ml_p}% of the time.

The shadow plates hold few transients throughout, in sunlight as well as in shadow. On the {n:sh.partial} plates the shadow partly covers, the shadowed areas hold {n:sh.own_obs}. At each plate’s own average density, they would be expected to hold about {n:sh.own_exp}. Almost the whole shortfall therefore comes from the plates themselves. A shadow can hide glints only within its outline, yet these plates are just as sparse in their sunlit parts. The shadow is centred on the ecliptic, and during the survey it fell mostly on sparse parts of the catalogue: of the {n:sh.plates} shadow plates, {n:sh.taurus_lc} lie towards Taurus and Perseus, where the catalogue holds unusually few transients, and {n:sh.eq_plates} are on the equatorial row.

The pointing model allows for the equatorial row and for dust across each plate. Even so, {n:sh.empty_n} of the shadow plates, {n:sh.empty_plates}, hold under a tenth of what their pointing predicts, and none of our modelled factors explain why. Watters et al. (2026) also noted these sparse plates. With their own random-shadow test, on a larger published list of transients that overlaps VASCO’s, they found that a shortfall as large as the one in their list arises by chance about 30% of the time.

Every shadow plate was taken with its centre {n:sh.shadow_lo} to {n:sh.shadow_hi} degrees from the point opposite the Sun. We therefore repeated the trials drawing only from the {n:sh.near_hosts} plates within {n:sh.near_max} degrees of that point, scaling for pointing as before. Random plates then match the shortfall {n:sh.p_near}% of the time ({n:sh.sig_near} standard deviations).

VASCO also repeated their test using only transients within 2 degrees of each plate’s centre, away from the dense edges. In that region the shadow holds {n:sh.centre_obs} against {n:sh.centre_exp} expected, a shortfall random plates match {n:sh.p_centre}% of the time.

Seen from Palomar rather than from the centre of the Earth, the shadow at geostationary height shifts by several degrees and crosses {n:sh.topo_plates} plates instead of {n:sh.plates}. In this version, the shadow holds {n:sh.topo_obs} against {n:sh.topo_exp} expected, a shortfall random plates match {n:sh.p_topo}% of the time ({n:sh.sig_topo} standard deviations).

### shadow-step2-caption
Comparing the transients in Earth's shadow, each counted by its rating, to two measurements of chance. Use the arrows to step through how each result is built.

Purple: VASCO’s model, in which every transient is independent.

Blue: the real Earth shadows laid on random plates across {n:sh.trials} trials, each scaled to the shadow plate's pointing.

### shadow-plates-title
Taking a closer look at the plates in shadow

### shadow-plates-text
Earth’s shadow crosses {n:sh.plates} of the catalogue’s plates. One of them, XE422, lies fully inside Earth's shadow, so it has no sunlit area to compare with. The other {n:sh.partial} plates are shown below, with the position of every transient from VASCO’s catalogue indicated.

The plates hold very different numbers of transients. XE415 holds {n:pl.XE415_n} transients, while XE355 only holds {n:pl.XE355_n}. Within a plate, transients bunch into clumps and gather along some edges, while some plates are almost empty in places. For example, XE640’s transients crowd along its northern edge and are missing from two large rounded regions. The sky has no reason to be empty in two round patches of a single plate, so these gaps most likely come from the plate itself, or from how it was scanned and searched.

Earth’s shadow splits each plate into two sections. On XE617, the shadow covers {n:pl.XE617_area}% of the plate and holds {n:pl.XE617_O} of its {n:pl.XE617_n} transients. Its rated count in shadow is {n:pl.XE617_Ow} against {n:pl.XE617_E} expected, so this plate holds more in shadow than expected. On XE640, the shadow clips {n:pl.XE640_area}% of the plate, in a corner that holds none of its {n:pl.XE640_n} transients.

### shadow-plates-caption
The {n:sh.partial} plates that Earth’s shadow partly covers, showing the central 6 by 6 degrees of each. Each dot is a transient in VASCO’s catalogue, shaded by its rating: fainter dots are rated less likely to be genuine. Transients in Earth's shadow are shown in red, and those in sunlight blue. The red line marks the edge of Earth's shadow, {n:sh.deg} degrees from the point opposite the Sun when the plate was taken.

Under each plate are two numbers: the rated count in shadow, and the count expected if rated transients there were as dense as in the sunlit region. The expected count allows for the higher density towards plate edges, measured on plates the shadow never crosses.

The ratio is the observed count in shadow divided by the expected count. A ratio of 1 means no difference in density between the shadow and sunlit regions.

### shadow-fit-title
Fitting the shadow transient density

### shadow-fit-text
VASCO report that transients are rarer in Earth’s shadow. We express this as a ratio, <em>r</em>: the density of transients in the shadowed part of a plate divided by the density in its sunlit part. With no shadow effect, <em>r</em> should be 1. If some transients genuinely vanish in the shadow, <em>r</em> should be below 1. VASCO’s analyses pool transients across plates rather than comparing shadow with sunlight on each plate, so they do not measure <em>r</em>. Because <em>r</em> compares each plate with itself, it does not depend on which plates the shadow fell on.

On each plate, <em>r</em> compares the share of its rated transients that fall in the shadowed part with the share that part would hold if the shadow made no difference. That expected share is not simply the shadowed fraction of the area, because transients are denser towards the edges and corners of every plate (see the plate maps above), and the ratings themselves vary across a plate. We therefore weight each part of the plate by the typical density of rated transients there, measured on the {n:sh.never} plates the shadow never crosses, so that a real shortfall in shadow cannot affect the weights.

Take {n:fit.ex_plate}, where Earth's shadow covers {n:fit.ex_area}% of the plate. With no shadow effect, the shadowed part would hold {n:fit.ex_weighted}% of the plate’s rated count. It actually holds {n:fit.ex_share}%, a rated count of {n:fit.ex_O} out of {n:fit.ex_n}. If transients vanished in shadow, this share would fall below {n:fit.ex_weighted}%, and <em>r</em> measures by how much.

Because each plate is compared only with itself, <em>r</em> is unaffected by how many transients a plate holds or how highly transients are rated overall.

The single value of <em>r</em> that best fits all {n:sh.partial} plates together is {n:fit.r}: transients are as dense in shadow as in sunlight. If every transient were an independent event, counting statistics would give a 68% range for <em>r</em> of {n:fit.count_lo} to {n:fit.count_hi}. This range allows for the ratings: a plate of many low-rated transients carries less weight than its raw count suggests.

The plates themselves show that this uncertainty is too small: single-plate ratios run from 0 on XE640 to {n:fit.max_r} on {n:fit.max_plate}, a much wider spread than counting allows. In {5_circle_symbol} we measure the true uncertainty, using trials built from the catalogue’s other plates.

### shadow-fit-caption
Grey points: the density of rated transients in shadow divided by the density in sunlight, for each plate the shadow partly covers. Error bars are taken from counting statistics. Plates with nothing in shadow have only upper limits.

Blue: the best-fit ratio across all {n:sh.partial} plates, with a one-standard-deviation range taken from counting statistics.

### shadow-toys-title
Finding the true uncertainty

### shadow-toys-text
Counting statistics treat every transient as independent. To measure the true uncertainty, we ran the same fit on plates where no shadow effect is possible. Each trial lays the {n:sh.partial} real shadow outlines onto {n:sh.partial} other plates that Earth’s shadow never crossed, and counts the rated transients inside and outside the shadow. Each random plate has a rated count similar to that of the shadow plate it replaces.

Over {n:fit.trials} trials, the fit gives a typical <em>r</em> of {n:fit.null_median}, with the middle 68% of results between {n:fit.null_lo} and {n:fit.null_hi}. The fit to the {n:sh.partial} real shadow plates, <em>r</em> = {n:fit.r}, lies well within this range, so it shows no shortfall in shadow beyond chance. The same trials give the range of true ratios consistent with the data: {n:fit.ci68_lo} to {n:fit.ci68_hi} at 68% confidence, and {n:fit.ci90_lo} to {n:fit.ci90_hi} at 90%.

The same trials show that our fit can detect a genuine shortfall. In a second set of trials, we removed transients inside the shadow outlines until they were only {n:fit.inj_r} times as dense there, a shortfall as large as VASCO’s. The fit then returns {n:fit.r} or higher only {n:fit.inj_pct}% of the time. VASCO did not measure the ratio plate by plate, but a shortfall of the size their count implies is disfavoured at {n:fit.inj_cl}% confidence. The test is one-sided: it asks only whether <em>r</em> could be as low as {n:fit.inj_r}.

These results hold when the random plates are selected differently, when each shadow is instead moved around on its own plate, or when any one plate is left out.

### shadow-toys-caption
Blue: fit result to the real plates.

Grey: {n:fit.trials} trials using the real shadows applied to random plates.

Purple: trials with a shortfall of transients in Earth's shadow added to the data.

Use the arrows to step through how this comparison is built.

### shadow-bands-title
Transients rated most likely to be genuine

### shadow-bands-text
VASCO’s 2026 work splits the transients into ten bands by rating, each holding a tenth of the transients, and finds the largest shortfall in shadow in the top band. Their test counts each band against an expected share of the sky, so it has the same problem as the count in {2_circle_symbol}: the shadow plates hold fewer transients. Under VASCO’s binomial model, every band falls short, by {n:sh.band_z_lo} to {n:sh.band_z_hi} standard deviations, as a shortage across whole plates would predict.

The shortfall is largest in the top band because the rating model tends to give high ratings on plates that hold many transients, and the shadow plates hold few. Of the {n:sh.partial} shadow plates, {n:fit.top_none} hold no top-band transients at all, in sunlight or in shadow. Of those, {n:fit.top_none_sunlit} lie more than 90% in sunlight, so the shadow cannot explain why they lack top-band transients.

Each band holds an unknown mix of glints (which would vanish in shadow) and other transients such as flaws (which would not vanish). Inside the shadow, only the non-glint transients remain, so a band’s ratio equals one minus its share of glints. We fit all ten bands at once, each with its own share of glints as a free parameter. Allowing every band some share of glints improves the fit no more than chance does: in {n:fit.joint_p}% of trials with no glints, the fit improves at least as much (twice the log-likelihood, a standard measure of fit quality, rises by {n:fit.joint_q} against a typical {n:fit.joint_med}).

If glints were among the transients the rating model rates as genuine, their share should rise with rating. The fitted share of glints rises no more than chance allows. It scatters from band to band, with a slight upward drift that trials with no glints match or exceed {n:fit.trend_p}% of the time. The two lowest bands, which the rating model judges mostly flaws, show as large a fitted share of glints as most of the bands above them.

The top band, where VASCO find their largest shortfall, also gives our lowest ratio, {n:fit.top_r}. We tested this three ways, and none gives a significant shortfall:

- Plate by plate, the ratio rests on {n:fit.top_O} transients in shadow, out of {n:fit.top_n}. Trials with no effect give a ratio this low {n:fit.top_p}% of the time.
- Over all {n:sh.plates} shadow plates, as in {2_circle_symbol}, the top band holds {n:sh.top_obs} transients in shadow against {n:sh.top_exp} expected. Tested on its own, that shortfall would arise by chance {n:sh.top_p}% of the time. But the band was picked out as the most extreme of ten, and the most extreme of ten bands looks this extreme {n:sh.top_any}% of the time (look-elsewhere effect).
- We also fitted a signal spread across the bands, like the transients VASCO labelled as genuine. That fit finds a shortfall this large in {n:fit.labels_p}% of trials with no signal.

All three tests are limited by the small number of top-band transients on the shadow plates. As such, the top band can neither confirm nor rule out a signal.

### shadow-bands-caption
Grey points: for each band of VASCO’s ratings, the density of transients in shadow divided by the density in sunlight, fitted across the {n:sh.partial} plates the shadow partly covers. Each transient counts once. Thick and thin bars show the 68% and 90% ranges from trials on plates the shadow never crossed.

Purple: the top band, where VASCO find their largest shortfall.

Blue: the same fit to all transients together, each counted by its rating, as in {4_circle_symbol} and {5_circle_symbol}.

### shadow-verdict-title
No significant shadow shortfall

### shadow-verdict
VASCO’s binomial model puts the shortfall of transients in Earth’s shadow at {n:sh.binom_z} standard deviations. Measured against how transient counts vary between plates, the shortfall amounts to {n:sh.sig_point} standard deviations. This shortfall reflects the plates the shadow happened to fall on, which tend to hold fewer transients in sunlight as well as in shadow.

Comparing shadow with sunlight on the same plates removes the influence of which plates the shadow crossed. In this test, rated transients are as dense in shadow as in sunlight: the best-fit ratio is {n:fit.r}, with a 90% confidence range of {n:fit.ci90_lo} to {n:fit.ci90_hi}. No band of VASCO’s ratings, including the top band, shows a shortfall beyond chance.

Earth’s shadow also crosses {n:sh.south} plates south of the celestial equator, beyond the edge of VASCO’s catalogue. Adding these plates would allow the same test to be performed on around {n:sh.south_ratio} times as much shadowed sky as now.

### nuclear-h2
Do transients cluster around nuclear tests?

### nuclear-intro
During the sky survey, the United States tested nuclear weapons in Nevada and the Pacific. In two papers (Scientific Reports, 2025 and 2026), VASCO report that observing nights within a day of a nuclear test held more transients than other nights. The 2026 paper uses US tests and counts each transient by its machine-learning rating; it puts the chance of the excess at about 1 in 40 (p = 0.024). It keeps the window of one day either side from the 2025 paper, as a hypothesis set before the analysis. The authors take the excess being confined to so narrow a window as a sign that it is genuine.

Repeating their measurement with the full list of US tests (VASCO’s list leaves out the {n:nuc.missing} tests of spring 1953), we find a similar excess. It reduces to a level consistent with chance once two further factors are taken into account:

- Nights near tests tended to have more photographs, and each photograph adds its own transients.
- On nights close to tests, the telescope tended to point at parts of the sky where more transients survive the catalogue’s selection, mostly away from the Milky Way. Our pointing model (see the section on transients per plate) expects these plates to hold about {n:nc2.exp_pct}% more transients than others, counted by rating.

### nuc-vasco-title
Repeating VASCO’s test

### nuc-vasco-text
VASCO’s test counts the transients on each observing night and compares the nights within a day of a US test with all other nights. It then asks how often the same number of nights, selected at random, hold at least as many transients.

We first repeat VASCO’s test as they set it up in their 2026 paper, with the archive dates, their list of test nights, and each transient counted by its rating. The {n:nuc.vasco_nights} nights near tests hold {n:nuc.vasco_ratio} times as many transients as the average of other nights. Random nights do at least as well about one time in {n:nuc.vasco_odds} ({n:nuc.vasco_sig} standard deviations), close to VASCO’s 1 in 40.

Our version of the test applies two corrections. We use the corrected dates described above, and the full list of US tests from the Department of Energy (DOE/NV-209). VASCO’s 2026 paper cites the same report. The list in their code repository gives the paper’s {n:nuc.vasco_nights} nights near tests, but it leaves out all {n:nuc.missing} tests of Operation Upshot–Knothole, in spring 1953.

As in VASCO’s test, a night is dated by the local date on which it began. With both corrections, {n:nuc.nights} nights fall within a day of a test. We find that they hold {n:nuc.ratio} times as many transients as the average of other nights, and random nights do as well about one time in {n:nuc.odds} ({n:nuc.sig} standard deviations).

### nuc-vasco-caption
The excess of transients, counted by rating, on the {n:nuc.nights} nights within a day of a US test (red line), compared to the same number of nights drawn at random {n:nuc.trials} times from the survey (grey). The analysis is performed using corrected dates and the full list of US tests.

### nuc-corr-title
Counting photographs, and considering telescope pointing

### nuc-corr-text
VASCO’s test counts transients per night, but a night’s transients come from its photographs. The {n:nuc.nights} nights near nuclear tests averaged {n:nuc.plates_window} plates per night, compared to {n:nuc.plates_other} on other nights. Part of this gap comes from nights with no plates in the catalogue: {n:nc2.empty_near}% of the nights near tests, against {n:nc2.empty_other}% of other nights. A night with no plates contributes zero transients, even though nothing was observed. Watters et al. (2026) make the same point: counted over only the days with plates north of the celestial equator, the link between transients and test days is no longer significant. Counted per plate rather than per night, the excess falls from {n:nc2.night_r} to {n:nc2.photo_r}.

The sky above Palomar changes through the year as Earth orbits the Sun. On spring nights, the telescope pointed away from the Milky Way. There, fewer stars crowd the field, so fewer detections are rejected for lying near a known star, and more survive as transients. Plates within 5 degrees of the Milky Way hold a median of {n:crowd.inner} transients, compared to {n:crowd.outer} in plates more than 40 degrees away. Of the plates near nuclear tests, {n:nc2.spring}% were taken between February and June. Allowing for where every plate pointed, our pointing model expects the plates near tests to hold around {n:nc2.exp_pct}% more rated transients than others. The transient ratings widen this gap, because they rise with distance from the Milky Way. Counting each transient once, the pointing model expects {n:nc.exp_pct}% more.

A like-for-like comparison uses the pointing model detailed earlier: the number of transients each plate’s pointing alone would leave in the catalogue. Against these expected counts, the excess for plates near tests falls from {n:nc2.photo_r} to {n:nc2.point_r}. Fitting the pointing model without the plates near tests, so they cannot shape their own expected counts, gives {n:nc2.held_out}.

To judge how often chance alone gives a result like {n:nc2.point_r}, VASCO’s test draws the same number of nights at random from across the survey. Random nights do at least as well {n:nc2.point_pv}% of the time ({n:nc2.point_pv_sig} standard deviations).

Nights near tests are not a random set, however: the tests came in series, so those nights bunch into a few observing runs, which share batches of emulsion and similar conditions. As a stricter check that keeps this bunching, we shift the real test calendar by 1 to 5 years and up to two months, leave out plates taken within two weeks of a real test, and score each shifted calendar in the same way. Shifted calendars match or exceed the result {n:nc2.point_pc}% of the time ({n:nc2.sig} standard deviations). The shifted calendars centre slightly below 1 (see chart below), because plates from the same months of the year hold a little less than their pointing predicts.

Most of the change in statistical significance comes from the pointing model. Counted per photograph without it, shifted calendars match or exceed the excess only {n:nc2.photo_pc}% of the time, against {n:nc2.point_pc}% with it. We checked that the pointing model does not also remove a genuine excess. We add a fake excess, multiplying the number of transients on every plate near a test by {n:nc2.inj_f}, and repeat the whole measurement. The pointing model is refitted to all plates, fake excess included, so it has every chance to absorb the excess. The plates near tests still hold {n:nc2.inj_r} times their expected number, so the fake excess of {n:nc2.inj_f} is recovered.

With VASCO’s own list of US tests in place of the full one, the result is {n:nuc.vlist_r}, which shifted calendars match or exceed {n:nuc.vlist_pct}% of the time. With the full list of US tests, any excess of transients near tests is found to be below {n:nc2.lim}% at 95% confidence. If there were no excess, the test would typically set a limit of about {n:nc2.lim_exp}%. Our limit is looser because the measured {n:nc2.point_r} lies a little above the centre of the shifted calendars.

VASCO’s ratings depend on which examples were used to train the rating model. We retrained it with VASCO’s code on {n:ml.label_sets} different sets of their published labels. Our result of {n:nc2.point_r} then moves only between {n:ml.point_lo} and {n:ml.point_hi}, and shifted calendars match or exceed it {n:ml.point_plo}% to {n:ml.point_phi}% of the time. Counting each transient once, without ratings, gives {n:nc.point_r}. In VASCO’s own version of the test, without the pointing model, the same retraining moves the chance of random nights doing as well from one time in {n:ml.vasco_odds_hi} to one time in {n:ml.vasco_odds_lo}.

The choice of training labels matters most among the transients the rating model rates highest. The top band of ratings holds {n:ml.top_r} times its expected number near US tests. Two details of the training bear on this. VASCO’s labelling code picks part of each batch for review from low-rated transients on test dates, as examples of flaws. Several of the rating model’s inputs are also averaged over a whole night, so the model can learn to rate particular nights higher or lower than others.

Retrained without the {n:ml.window_labels} labels from nights near tests, the top band holds {n:ml.top_now} times its expected number. Retrained without the night-level inputs, it holds {n:ml.top_flat_lo} to {n:ml.top_flat_hi} times its expected number, depending on the labels used. Neither retrained result is significant, nor is the original {n:ml.top_r}: it is the most extreme of ten bands, and with ten bands to choose from, chance gives a band this extreme about one time in {n:ml.top_any_odds}.

### nuc-corr-caption
Once the number of photographs and telescope pointing are accounted for, the plates near nuclear tests hold {n:nc2.point_r} times their expected number of transients.

Compared with {n:nc2.slides} shifted test calendars, each catching about as many plates as the real one, about {n:nc2.frac} trials show at least the same level of excess.

### nuclear-verdict-title
No significant link to the nuclear tests

### nuclear-verdict
We reproduce VASCO’s excess of transients near nuclear tests, but find that it is primarily due to where the telescope pointed. Once pointing is allowed for, no significant effect of the nuclear test schedule on transient counts is found.

### align-h2
Do transients form straight lines?

### align-intro
A flat object tumbling in orbit would flash each time it caught the Sun, leaving a line of points during a long plate exposure. VASCO searched for such alignments and shortlisted five candidates. Their best candidate was reported at 3.9 standard deviations, corresponding to an estimated 1 in 10,000 probability of arising by chance.

A catalogue of 100,000 unevenly spread points may contain many chance alignments. Three questions naturally arise:

- Do the alignment quality checks used introduce biases in the directions of the alignments they can find?
- Does the real transient catalogue contain more alignments than expected by random chance?
- How many genuine glint tracks could the catalogue hold?

### align-step1-title
Searching for alignments in all directions

### align-step1-text
VASCO kept a set of transients only if its two sky coordinates, right ascension (east–west) and declination (north–south), rise and fall together almost perfectly: a Pearson correlation coefficient above 0.99. This approach works well for diagonal lines. However, a line running close to east–west barely changes in declination along its length. The small scatter in each measured position then hides that change, so the correlation falls below the cut and the line is rejected. The same happens to lines close to north–south.

In simulations (image below), the VASCO method keeps {n:align.pear_diag}% of diagonal lines, {n:align.pear_2}% of lines 2 degrees from east–west, and only {n:align.pear_ew}% of lines running exactly east–west. A satellite at geostationary height stays nearly fixed above the observatory while the telescope follows the stars, so on the plate it moves close to east–west. Seen from Palomar, the geostationary belt itself lies just south of the catalogue’s plates, so any such object on them is on a slightly tilted orbit, and its track runs within a few degrees of east–west. VASCO’s selection would reject almost all of the tracks such an object would leave.

The same alignment task often arises in particle physics, where the paths of charged particles must be reconstructed from hits they leave in a detector. The simplest tracking method fits a straight line through a set of points, resulting in a track candidate. The quality of each candidate can be assessed by checking how far the points lie from the line of best fit. Such a method works equally well in all directions, including east–west when applied to this setting.

### align-step1-caption
The share of simulated straight lines kept by each method, shown as a function of direction in the sky. VASCO’s selection (orange) collapses along east–west, where genuine geostationary satellite glints would lie, as well as along north–south. A straight-line fit (blue) keeps nearly all lines in every direction.

### align-step2-title
Fitting every track candidate

### align-step2-text
VASCO’s five shortlisted alignments came from an earlier list of 298,165 candidate transients (Solano et al. 2022). Of their {n:align.vasco_transients} transients, {n:align.vasco_in_cat} appears in the later catalogue of {n:cat.transients} transients that VASCO use for their other tests. We have applied our track-fit method to the full catalogue of {n:cat.transients} transients.

We build track candidates from pairs of transients up to 15 arcminutes apart. Each pair marks the two ends of a track, and the transients lying closest to the line between them are added, giving three to six transients in total. We then fit a straight line through each candidate, and measure how closely its transients follow it.

Two sources closer than about 4 arcseconds (the width of a star’s image on the scans) blur into one, so we ignore candidates with points closer than 8 arcseconds. The released catalogue itself holds almost no transients closer together than this. {n:align.vasco_close_n} of VASCO’s five shortlisted alignments contain points only {n:align.vasco_close} arcseconds apart, closer together than almost any pair in the released catalogue. Watters et al. (2026) also find that 8 of the 24 transients in the five alignments lie within 5 arcseconds of sources seen by the NEOWISE infrared survey, so may be stars rather than transients.

The position of each transient on a scan is uncertain by about an arcsecond. As such, even points from a genuine straight track can stray from a true straight line. Each track fit is scored with a χ² per degree of freedom (χ²/ndof), which compares the points’ actual distances from the line with a one-arcsecond uncertainty. A genuine track gives a value near 1 or below. Points scattered loosely about the line give a much larger value.

A good fit alone does not make a track stand out: in crowded parts of a plate, pure chance can cause transients to align frequently. We therefore assign each track candidate a track score, which measures how easily chance alone could produce it. The score uses four numbers, as explained in the diagram below: the candidate’s number of transients (<em>k</em>), its length (<em>L</em>), its width (<em>w</em>), and the number of transients in a box at least 10 arcminutes across around it (<em>N</em>).

The idea and the formula match those VASCO used to judge their own candidates (Edmunds and George, 1985). We apply our method in all directions, measure straightness from the fitted line, and score only track candidates that are straight to within the measurement error (a χ² probability of at least 1%). The lower the track score, the less often chance alone can produce such a track.

Transients often come in close pairs: the catalogue holds far more pairs of transients within 20 arcseconds of each other than chance would give. Two transients this close fix a direction only loosely, so a third lying roughly in line with them makes a short track that looks stronger than it is. Following Lezama and colleagues (2015), who developed a method for finding alignments among clustered points, we count a track’s transients in boxes 20 arcseconds long, the distance over which close pairs are unusually common: transients closer together than this along the track count once.

Unlike the shadow and nuclear tests, the track search counts every transient once. A track is a pattern made by several detections together, and a rating judges each detection on its own image, so it has no natural place in the track score.

Candidates built from a set of six transients also yield candidates composed of any combination of three, four, or five of that set. Where track candidates share a transient, we keep only the best-scoring candidate, such that each transient counts once.

### align-method-caption
How the track score rates a real track candidate on plate {n:align.ex_plate}. A candidate is scored only if it first passes the straightness test.

The track score, <em>μ</em>, is the number of tracks with <em>k</em> transients, at most <em>L</em> long and lying within <em>w</em> of a straight line, that we would expect if the <em>N</em> transients in the box were scattered at random, with <em>n</em> their density (Edmunds and George 1985):

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
Top: track score distributions for all 3-, 4-, 5-, and 6-transient candidates in the catalogue. The real data distribution is shown in red, and the expected random background distribution is shown in grey.

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
Using our track score method, we searched the full catalogue for glint tracks of 3, 4, 5, and 6 transients along straight lines. We find no excess above the background level expected from random catalogues.

With 95% probability, the catalogue holds at most {n:align.lim6} 6-transient, {n:align.lim5} 5-transient, and {n:align.lim4} 4-transient glint tracks. To the best of our knowledge, these are the first such limits to be set on glint tracks.

### end-h2
Summary

### card-4-question
Are the plates correctly dated?

### card-4-ours
Around {n:dates.early_pct}% of plates have an archive date that is one day early. Asteroid trails, Palomar’s own log, and the handwritten dates on plates independently corroborate this finding. We provide a correction to the dates in our code repository.

### card-1-question
Do transients avoid Earth’s shadow?

### card-1-ours
Counting each transient by VASCO’s rating and allowing for where each plate points, chance gives a shortfall in shadow as large as the one seen {n:sh.p_point}% of the time. The shortfall comes almost entirely from the plates the shadow fell on. Plate by plate, transients are as dense in shadow as in sunlight, and no band of the ratings shows a shortfall beyond chance. A shortfall as large as VASCO’s count implies is disfavoured at {n:fit.inj_cl}% confidence (one-sided).

### card-2-question
Do transients cluster around nuclear tests?

### card-2-ours
Nights near tests had more photographs, and pointed at parts of the sky where more transients survive the catalogue’s selection. Counted per photograph and allowing for pointing, with each transient counted by VASCO’s rating and the full list of US tests, plates near tests hold {n:nc2.point_r} times their expected number. Shifted test calendars match or exceed this {n:nc2.point_pc}% of the time. Retraining the ratings on different sets of VASCO’s labels moves the result only between {n:ml.point_lo} and {n:ml.point_hi}. The top band of ratings moves more, from {n:ml.top_now} to {n:ml.top_r}, but no choice of labels makes either result significant.

### card-3-question
Do transients form straight lines?

### card-3-ours
Our track score, which looks in every direction and counts each close group of transients once, finds as many track candidates in the released catalogue as chance predicts. We set the first limits on the number of glint tracks it holds; VASCO’s five shortlisted alignments come from an earlier list and are not tested here.

### end-text
Our corrected plate dates, track score, and background-only model are in our code repository, free for anyone to use. We would welcome the VASCO team’s response to this analysis, and thank them again for making their data and code available.

### footer-about
The analysis is our own code, run on the VASCO team’s public data from Brian Doherty’s VASCO-ML and VASCO-challenge repositories: their transient catalogue and its machine-learning ratings, random positions on each plate, and lists of observing nights and nuclear tests. Where we reproduce VASCO’s results, we follow their published method and code. Asteroid positions come from JPL Horizons, star positions from Gaia, plate times from Palomar’s log (VizieR VI/25), and test times from the US Department of Energy report DOE/NV-209. Every result of our own analysis can be reproduced from our repository, which also holds the table of corrected plate dates; figures quoted from VASCO and from other authors are taken from their papers. Watters et al. (2026, PASA) re-examined the shadow and nuclear-test results using other published lists of transients, and the VASCO team has responded (Villarroel et al. 2026, arXiv:2602.15171). Their response argues that Watters et al. used a reduced, filtered subset of transients, lacking complete timing information, too small for the shadow test, and without uncertainties. The analysis here uses VASCO’s full released catalogue, every plate’s corrected date and time, and uncertainties from random trials built on the real data. We thank Beatriz Villarroel and her colleagues for publishing their data and code; their openness made this work possible.
