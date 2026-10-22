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

VASCO considered objects in geostationary orbit (36,000 kilometres up), where a single orbit takes one day and a satellite stays above a point on the equator. At this height, Earth’s shadow covers a disc of sky {n:sh.diam} degrees across, centred on the point opposite the Sun. If some transients are glints from objects in orbit, the part of a plate inside the shadow should show fewer transients than the same area of sky elsewhere.

VASCO place the shadow as seen from the centre of the Earth. The telescope, however, sat on the Earth’s surface, 6,400 kilometres from its centre, and from there an object 36,000 kilometres up appears in a different direction. On the plates the shadow reaches, it falls {n:sh.shift_lo} to {n:sh.shift_hi} degrees from where the Earth-centred picture puts it. We therefore place the shadow as seen from Palomar.

Seen from Palomar, the shadow also moves during an exposure: as the Earth turns, the observatory travels more than 1,000 kilometres in 50 minutes, and the shadow’s edge slides {n:sh.sweep_lo} to {n:sh.sweep_hi} degrees across the plate. A glint lasts a fraction of a second, at an unknown moment during the exposure, so we count each transient by the share of the exposure that its place on the plate spent in shadow. The shadow’s edge is not sharp either: between {n:sh.umbra} and {n:sh.penumbra} degrees from its centre, the Earth hides only part of the Sun, and a glint there would be dimmer. Following the shadow through each exposure requires accurate plate times. On {n:dates.ast_mid_plates} plates, asteroid trails fix the middle of the exposure independently, and they typically agree with our times to within a minute, in which the shadow moves by less than a tenth of a degree.

Flaws in the photographic emulsion, dust on the glass, and marks from copying have no connection with where Earth’s shadow fell. If fewer transients fall inside the shadow than the differences between plates can explain, that would point to something in orbit.

Compared with VASCO’s test, our approach:

- Places the shadow where it fell on each plate, as seen from Palomar, and follows it through each exposure.
- Uses each plate’s corrected date and time.
- Measures chance by laying the real shadows on plates no shadow reaches, allowing for where each plate points, rather than by treating every transient as independent.
- Adds a test that compares shadow with sunlight on the same plates, which does not depend on how many transients a plate holds in total.
- Repeats both tests in each band of VASCO’s ratings.

### shadow-step1-title
Plates intersecting Earth's shadow

### shadow-step1-text
As seen from Palomar, Earth’s shadow reaches {n:sh.plates} of the {n:cat.plates} plates in the catalogue during their exposures. Throughout this section, we use the central 6 by 6 degrees of each plate, matching the plate size VASCO use in their own analysis. VASCO’s random positions cover this square evenly on every plate, as shown above.

We first reproduce VASCO’s own count, which uses their archive dates, counts each transient once, and places a shadow with a sharp edge, seen from the centre of the Earth, at a single moment for each plate. It holds {n:sh.vasco_obs_geo} transients against {n:sh.vasco_exp_geo} expected, which VASCO’s binomial model, described in {2_circle_symbol}, puts at {n:sh.vasco_z_geo} standard deviations (close to their published 7.6). Seen from Palomar instead, with nothing else changed, the shadow reaches {n:sh.vasco_plates} of the plates in their tests and holds {n:sh.vasco_obs} transients against {n:sh.vasco_exp} expected, about a fifth fewer: {n:sh.vasco_z} standard deviations.

VASCO’s count also misreads the plate times. Each plate’s header gives the start of the exposure in hours and hundredths of an hour, but the count reads it as hours and minutes, so 06:82 becomes 07:22 rather than 06:49, and many exposures are placed tens of minutes from their true time. From the centre of the Earth this barely matters, but as seen from Palomar it moves the shadow by up to about a degree.

Using our corrected dates and times, the parts of these plates that spent time in shadow hold {n:sh.transients} transients, each counted by its share of the exposure in shadow. Together they cover {n:sh.worth} plates’ worth of sky: all of plate XE422, most of XE472, about half of XE415, XE369, and XE532, and slivers of four others. An average plate holds {n:sh.count_avg} transients in its central 6 by 6 degrees, so if transients were spread evenly across the survey, the shadowed parts should hold about {n:sh.count_exp}.

Counted by rating, the transients in shadow on our corrected dates give a rated count of {n:sh.observed}, against {n:sh.expected} for the same area of an average plate. From here on, counts in this section are rated counts, with each transient also counted by its share of the exposure in shadow.

### shadow-step1-caption
Of the {n:cat.plates} plates in the catalogue, as seen from Palomar, Earth's shadow reaches {n:sh.plates} during their exposures. One plate, XE422, lies inside the shadow throughout. Each circle marks the shadow’s edge at the middle of a plate’s exposure.

### shadow-step2-title
Measuring the probability of this deficit

### shadow-step2-text
VASCO score their shortfall with a binomial test, which treats every transient as an independent event. Their headline figure, 22 standard deviations, compared the count in shadow with the shadow’s share of the northern sky; their figure of 7.6 standard deviations used the area of each plate the shadow actually covers, as we do here. Applied to our rated counts, the binomial model expects {n:sh.expected} transients in shadow, give or take {n:sh.binom_sd}, so the observed count of {n:sh.observed} lies {n:sh.binom_z} standard deviations below.

Transients are not independent, however: as shown above, counts on plates differ far more than counting noise permits, so the count inside a shadow depends on whether the plates it falls on happen to be rich or poor. To measure chance directly, we laid the {n:sh.plates} real shadows, with their shapes, positions, and movement during the exposure, onto plates chosen at random from those no shadow reaches. Each count is scaled to what the pointing model, fitted to rated counts on those plates, predicts at the shadow plate’s pointing. As VASCO recommend, the test works with total counts: no plate’s own count of transients sets what it is expected to hold.

Across {n:sh.trials} such trials, a count of {n:sh.observed} or fewer turned up {n:sh.p_point1}% of the time, a significance of {n:sh.sig_point} standard deviations; without the pointing model, {n:sh.p_random}%. Our trials measure chance correctly: on plates no shadow reaches, treated as if it had fallen there, a result below 5% came up {n:sh.cal_5}% of the time. This comparison also has limited power: had transients been as scarce in shadow as VASCO’s own count implies, it would have reached 95% significance only {n:sh.power}% of the time.

Most of the observed shortfall comes from two plates, XE472 and XE532, both of which exhibit transient levels well below expected given their pointing: each holds about {n:pl.XE472_qsu} of what its pointing predicts in sunlight, and {n:pl.XE472_qsh} in shadow. A shadow can hide glints only inside its outline, so the shadow did not empty these plates. They lie in a sparse patch of the catalogue towards Taurus and Perseus, where plates the shadow never reaches also hold far less than their pointing predicts, and which our pointing model does not model.

If each shadow is laid only on plates within 12 or 20 degrees of its own, the shortfall is found at random {n:sh.local_12}% or {n:sh.local_20}% of the time. Watters et al. (2026) also noted these sparse plates, and with their own random-shadow test, conducted on a larger published list of transients, found that a shortfall like the one in their list arises by chance about 30% of the time.

Other choices matter little. How the shadow’s edge is drawn, moving every plate’s time by up to five minutes, and sampling the exposure more finely keep the result between {n:sh.sys_p_lo}% and {n:sh.sys_p_hi}%. The {n:sh.ml_deg}-degree shadow radius in the code for VASCO’s 2026 paper gives {n:sh.ml_p1}%, and drawing random plates only from near the point opposite the Sun, where every shadow plate lies, gives {n:sh.p_near1}%.

Within 2 degrees of each plate’s centre, away from the dense edges, as VASCO also tested, the shadow holds {n:sh.centre_obs} against {n:sh.centre_exp} expected, which random plates match {n:sh.p_centre}% of the time.

The pooled count is therefore sensitive to which plates the shadow happened to fall on, and to how well the pointing model describes them. From {3_circle_symbol} onwards, we compare shadow with sunlight on the same plates, which removes both of these factors.

### shadow-step2-caption
Comparing the transients in Earth’s shadow, each counted by its rating and its share of the exposure in shadow, to two measurements of chance. Use the arrows to step through how each result is built.

Purple: VASCO’s model, in which every transient is independent.

Blue: the real Earth shadows laid on random plates across {n:sh.trials} trials, each scaled to the shadow plate’s pointing.

### shadow-plates-title
Taking a closer look at the plates in shadow

### shadow-plates-text
Earth’s shadow reaches {n:sh.plates} of the catalogue’s plates. One of them, XE422, lies inside the shadow throughout the exposure, so it has no sunlit area to compare with. Of the other {n:sh.partial}, seven are shown below, with the position of every transient from VASCO’s catalogue. The eighth, XE355, is left out: the shadow touches only a corner of it that holds no transients.

The plates hold very different numbers of transients. XE415 holds {n:pl.XE415_n} transients, while XE355 holds only {n:pl.XE355_n}. Within a plate, transients bunch into clumps and gather along some edges, while some plates are almost empty in places.

During each exposure, the shadow’s edge would have slid across the plate; in the diagrams below, the shadow’s edge at the start of each exposure is shown by a dashed line, and at the end by a solid line.

Between these two lines, the chance that a transient appeared while its location was in shadow rises smoothly from 0 to 1. On XE415, the shadow covers {n:pl.XE415_area}% of the plate over the exposure, and the rated count in shadow is {n:pl.XE415_Ow} against {n:pl.XE415_E} expected. On XE356, the shadow reaches only a sliver along one edge, late in the exposure, which holds {n:pl.XE356_Ow} transients against {n:pl.XE356_E} expected.

### shadow-plates-caption
Plates where Earth's shadow crossed during the exposure, as seen from Palomar.

Light areas were in sunlight for the whole exposure, and dark areas in shadow throughout. The band between them is where the shadow passed during the exposure: the darker the shading, the longer that area was in shadow. The dashed red line is the shadow’s edge at the start of the exposure, and the solid line is its edge at the end.

Each white dot is a transient. Below each plate, the rated count of transients within the shadow is given, along with the expected count if the shadow had no effect on transient levels. Transients falling in the band count in proportion to how long their position was in shadow.

### shadow-fit-title
Fitting the shadow transient density

### shadow-fit-text
VASCO report that transients are rarer in Earth’s shadow. We express this as a ratio, <em>r</em>: the density of transients where a plate is in shadow, divided by the density where it is in sunlight. With no shadow effect, <em>r</em> should be 1. If some transients genuinely vanish in the shadow, <em>r</em> should be below 1. VASCO’s analyses pool transients across plates rather than comparing shadow with sunlight on each plate, so they do not measure <em>r</em>. Because <em>r</em> compares each plate with itself, it does not depend on which plates Earth's shadow fell upon.

On each plate, <em>r</em> compares the share of its rated transients that fall in the shadowed part with the share that part would hold if the shadow made no difference. A transient in the band where the shadow passed during the exposure is split between the two: if its spot was in shadow for a third of the exposure, it counts one third towards shadow and two thirds towards sunlight.

Even with no shadow effect, the shadowed part of a plate will not hold a number of transients in proportion to its share of the plate's area. Transients are denser towards the edges and corners of every plate (see the plate maps above), and their ratings also vary across a plate. We therefore weight each part of the plate by the typical density of rated transients there, measured on the {n:sh.never} plates no shadow ever reaches, so that a real shortfall in shadow cannot affect the weights.

Take {n:fit.ex_plate} for example, where Earth’s shadow covers {n:fit.ex_area}% of the plate over the exposure. With no shadow effect, the shadowed part would hold {n:fit.ex_weighted}% of the plate’s rated count of transients. It actually holds {n:fit.ex_share}%, a rated count of {n:fit.ex_O} out of {n:fit.ex_n}. For this plate, that gives <em>r</em> = {n:pl.XE369_r}: transients are {n:pl.XE369_r} times as dense in shadow as in sunlight.

Because each plate is compared only with itself, <em>r</em> is unaffected by how many transients a plate holds, how highly they are rated overall, or what the pointing model predicts for it. XE472 and XE532, which drive the shortfall in {2_circle_symbol}, hold only {n:pl.XE472_n} and {n:pl.XE532_n} transients, so they carry little weight here. And since each plate is compared only with itself, whatever caused these two plates to be nearly empty (both in shadow and sunlight) does not enter the comparison.

The single value of <em>r</em> that best fits all {n:sh.partial} plates together is {n:fit.r}: transients are somewhat less dense in shadow than in sunlight. If every transient were an independent event, counting statistics would give a 68% range for <em>r</em> of {n:fit.count_lo} to {n:fit.count_hi}. This range takes into account the VASCO transient ratings: a plate of many low-rated transients carries less weight than its raw transient count suggests.

The plates themselves show that the uncertainty from counting statistics is too small, because single plates disagree far more than counting allows. XE415 gives a ratio of {n:pl.XE415_r}, with a counting range of {n:pl.XE415_lo} to {n:pl.XE415_hi}; XE356, where the shadow reaches only a sliver late in the exposure, gives {n:pl.XE356_r}, with a range of {n:pl.XE356_lo} to {n:pl.XE356_hi}. In {5_circle_symbol} we measure the true uncertainty, using trials built from the catalogue’s other plates.

### shadow-fit-caption
Grey points: the density of rated transients in shadow divided by the density in sunlight, for each plate the shadow partly covers. XE355, where Earth's shadow touches only an empty corner, is not considered. Error bars are taken from counting statistics.

Blue: the best-fit ratio across all {n:sh.partial} plates, with a one-standard-deviation range taken from naive counting statistics.

### shadow-toys-title
Finding the true uncertainty

### shadow-toys-text
Counting statistics treat every transient as independent. To measure the true uncertainty, we ran the same fit on plates where no shadow effect is possible. Each trial lays the {n:sh.partial} real shadows, with their movement during the exposure, onto {n:sh.partial} other plates that no shadow ever reached, and fits the transients there exactly as on the real plates. Each random plate has a rated count similar to that of the shadow plate it replaces. On simulated plates with a known ratio, the fit recovers the simulated ratio without bias, and its 68% counting range contains the simulated value {n:fit.closure_cov}% of the time, as expected.

Over {n:fit.trials} trials, the fit gives a typical <em>r</em> of {n:fit.null_median}, with the middle 68% of results between {n:fit.null_lo} and {n:fit.null_hi}. Compared to the fit to real data, which finds <em>r</em> = {n:fit.r}, trials find deficit values at least this extreme {n:fit.p_low}% of the time ({n:fit.sig_low} standard deviations). The same trials give the range of true ratios consistent with the data: {n:fit.ci68_lo} to {n:fit.ci68_hi} at 68% confidence, and {n:fit.ci90_lo} to {n:fit.ci90_hi} at 90%.

We find that the result rests predominantly on a single plate: without {n:fit.single_plate}, the fit measures <em>r</em> = {n:fit.drop_r}. On its own, {n:fit.single_plate} gives a ratio of {n:fit.single_r}. To see whether that is unusual, we moved its shadow to other positions on the same plate, leaving its transients where they are: {n:fit.single_p_moved}% of positions give a ratio at least this low. Laying the same shadow on other plates with a similar number of transients gives a ratio this low {n:fit.single_p_similar}% of the time. A ratio this low can arise by chance, based on the uneven distribution of transients across a plate.

However, our results hold when the random plates are selected differently, when each shadow is instead moved around on its own plate, and under multiple choices of shadow edge, plate times, and sampling of the exposure. All of these variations of the analysis return r</em> values between {n:fit.sys_r_lo} and {n:fit.sys_r_hi}.

### shadow-toys-caption
Blue line: fit result to the real plates.

Grey: {n:fit.trials} trials using the real shadows applied to random plates.

Shaded: trials with ratio values at least as low as the observed value in real data.

Use the arrows to step through how this comparison is built.

### shadow-bands-title
Transients rated most likely to be genuine

### shadow-bands-text
VASCO’s 2026 work splits the transients into ten bands by rating, each holding a tenth of the transients, and finds the largest Earth's shadow shortfall in the top band. Their test counts each band against an expected share of the sky, so it has the same problem as the count in {2_circle_symbol}: the shadow plates tend to hold fewer transients. Under VASCO’s binomial model, every band falls short, by {n:sh.band_z_lo} to {n:sh.band_z_hi} standard deviations, as a general shortage across entire plates would predict.

The shortfall is largest in the top band because the rating model tends to give high ratings on plates that hold many transients, and the shadow plates hold few. Of the {n:sh.partial} shadow plates, {n:fit.top_none} hold no top-band transients at all, in sunlight or in shadow. Of those, {n:fit.top_none_sunlit} lie more than 90% in sunlight, so the shadow cannot explain why they lack top-band transients.

Each band holds an unknown mix of glints (which would vanish in shadow) and other transients such as flaws (which would not vanish). Inside the shadow, only the non-glint transients remain, so a band’s ratio equals one minus its share of glints. We fit all ten bands at once, each with its own share of glints as a free parameter. Allowing every band some share of glints improves the fit no more than chance does: in {n:fit.joint_p}% of trials with no glints, the fit improves at least as much.

If the rating model picked out glints, the clearest sign would be a steady trend: the ratio falling further below 1 as the rating rises. The trend fitted through the ten bands (shown with the button above the chart) does fall, from {n:fit.trend_lo} in the lowest band to {n:fit.trend_hi} in the top band, but trials with no glints give a trend at least this steep {n:fit.trend_p}% of the time. The bands also scatter more than a steady trend would allow: the lowest band, which the rating model judges to contain mostly flaws, has one of the lowest ratios.

The top band, where VASCO find their largest shortfall, gives a standalone ratio of {n:fit.top_r}; the lowest ratio, {n:fit.band_lo}, comes from the {n:fit.band_low_ord} band. We tested the top band three ways, none of which indicate a statistically significant shortfall:

- Plate by plate, the ratio rests on {n:fit.top_O} transients in shadow out of {n:fit.top_n}. Trials with no effect give a ratio this low {n:fit.top_p}% of the time.
- Over all {n:sh.plates} shadow plates, as in {2_circle_symbol}, the top band holds {n:sh.top_obs} transients in shadow against {n:sh.top_exp} expected. Tested on its own, that shortfall would arise by chance {n:sh.top_p}% of the time. But the top band is one of ten, and here the {n:sh.top_lowest_ord} band’s shortfall is more significant. By chance alone, at least one of ten bands looks this extreme {n:sh.top_any}% of the time (the look-elsewhere effect).
- We also fitted a signal spread across the bands, like the transients VASCO labelled as genuine. That fit finds a shortfall this large in {n:fit.labels_p}% of trials with no signal.

All three tests are limited by the small number of top-band transients on the shadow plates. As such, the top band can neither confirm nor rule out a signal.

### shadow-bands-caption
The ratio as measured in each band of VASCO's transient ratings. Thick and thin bars mark the 68% and 90% ranges, found from trials on plates the shadow never reaches.

The blue point shows the ratio measured across all transients, each weighted by its rating, as in {4_circle_symbol} and {5_circle_symbol}.

Use the button to view the trend in the fitted ratio with rating score. The shaded bands show the 68% and 90% intervals.

### shadow-verdict-title
When measured at plate level, no significant shortfall in Earth's shadow

### shadow-verdict
VASCO’s binomial model puts the shortfall of transients in Earth’s shadow at {n:sh.binom_z} standard deviations. Measured against how transient counts vary between plates, we find the an effect amounting to {n:sh.sig_point} standard deviations. Most of the shortfall is driven by two plates, XE472 and XE532, which are low in transients both in sunlight and shadow. These plates fall in a sparse patch of the catalogue, which our pointing model applies no correction for.

We perform a direct measurement of the Earth's shadow effect, by measuring transients in shadow and sunlight at a per-plate level. We find that rated transients are somewhat less dense in shadow, with a best-fit ratio of {n:fit.r} ({n:fit.sig_low} standard deviations). The ratio predominantly rests on one plate, {n:fit.single_plate}, which exhibits an uneven spread of candidates capable of creating an effect by chance. Without it, the ratio measured is {n:fit.drop_r}. No single band of VASCO’s ratings shows a shortfall beyond chance.

Earth’s shadow also reaches {n:sh.south} plates south of the celestial equator, beyond the edge of VASCO’s catalogue. Adding these plates would allow the same tests to be performed on around {n:sh.south_ratio} times as much shadowed sky.

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
We place Earth’s shadow where it fell, as seen from Palomar, and count each transient by VASCO’s rating. Allowing for where each plate points, the shadow holds fewer transients than expected ({n:sh.sig_point} standard deviations). Most of this shortfall comes from two plates that hold few transients, including in their sunlit areas.

When comparing transients in sunlight and shadow at the per-plate level, a shortfall is also measured albeit less significant ({n:fit.sig_low} standard deviations). This result is largely driven by a single plate.

### card-2-question
Do transients cluster around nuclear tests?

### card-2-ours
Nights near tests had more photographs, and pointed at parts of the sky where more transients survive the catalogue’s selection. Counted per photograph and allowing for pointing, with each transient counted by VASCO’s rating and the full list of US tests, plates near tests hold {n:nc2.point_r} times their expected number. Shifted test calendars match or exceed this {n:nc2.point_pc}% of the time ({n:nc2.sig} standard deviations).

### card-3-question
Do transients form straight lines?

### card-3-ours
Our track score, which looks in every direction and counts each close group of transients once, finds as many track candidates in the released catalogue as chance predicts. We set the first limits on the number of glint tracks it holds; VASCO’s five shortlisted alignments come from an earlier list and are not tested here.

### end-text
Our corrected plate dates, track score, and background-only model are in our code repository, free for anyone to use. We would welcome the VASCO team’s response to this analysis, and thank them again for making their data and code available.

### footer-about
The analysis is our own code, run on the VASCO team’s public data from Brian Doherty’s VASCO-ML and VASCO-challenge repositories: their transient catalogue and its machine-learning ratings, random positions on each plate, and lists of observing nights and nuclear tests. Where we reproduce VASCO’s results, we follow their published method and code. Asteroid positions come from JPL Horizons, star positions from Gaia, plate times from Palomar’s log (VizieR VI/25), and test times from the US Department of Energy report DOE/NV-209. Every result of our own analysis can be reproduced from our repository, which also holds the table of corrected plate dates; figures quoted from VASCO and from other authors are taken from their papers. Watters et al. (2026, PASA) re-examined the shadow and nuclear-test results using other published lists of transients, and the VASCO team has responded (Villarroel et al. 2026, arXiv:2602.15171). Their response argues that Watters et al. used a reduced, filtered subset of transients, lacking complete timing information, too small for the shadow test, and without uncertainties. The analysis here uses VASCO’s full released catalogue, every plate’s corrected date and time, and uncertainties from random trials built on the real data. We thank Beatriz Villarroel and her colleagues for publishing their data and code; their openness made this work possible.
