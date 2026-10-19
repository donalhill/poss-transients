### hero-title
Examining the evidence for glints from orbit in 1950s sky photographs

### hero-lede
Between 1949 and 1957, Palomar Observatory in California photographed the northern sky on glass plates. The plates show faint points of light that have no counterpart in modern sky catalogues. The VASCO project, led by Dr Beatriz Villarroel, has catalogued more than 100,000 of these transients.

Villarroel’s team has proposed that some transients may be produced by sunlight glinting off flat, reflective objects in orbit around Earth in the pre-satellite era. In papers published in 2025 and 2026, the team reported three statistical patterns in support of this theory: fewer transients appear inside Earth’s shadow; more transients appear around the dates of nuclear weapons tests; some transients appear in straight patterns.

We are grateful to the team for publishing its catalogue and code. Building the catalogue took years of careful work, and releasing it openly made this re-examination possible. We repeat the shadow and nuclear-test analyses on their data, and apply both their alignment method and our own to the released catalogue.

### survey-h2
The Palomar Sky Survey

### survey-intro
Images were taken on plates using Palomar’s 48-inch Schmidt telescope, which was built to photograph wide areas of sky. Each plate is a sheet of glass 14 inches square, and covers a patch of sky 6.5 degrees across (about thirteen full Moons). Every patch was photographed through both a blue and a red filter, with red exposures lasting 40 to 60 minutes. Copies of the plates were scanned in the 1990s, and astronomers now work from those scans.

VASCO compared the red plates with three modern catalogues of the sky: Gaia, Pan-STARRS, and WISE. A point of light became a transient if no modern catalogue lists anything within 5 arcseconds of it. An arcsecond is a 3,600th of a degree, and 60 arcseconds make an arcminute. The public catalogue holds {n:cat.transients} such transients spanning {n:cat.plates} plates.

Many transients are likely to be artefacts, such as dust on the glass, specks in the photographic emulsion, or marks made in copying and scanning. The catalogue has not been checked by eye. The VASCO team has since trained a machine-learning model that rates how likely each transient is to be genuine, but the catalogue was built before this and includes every transient whatever its rating. The ratings are not yet public, so our analysis uses the whole catalogue.

### survey-map-caption
Each tile shown represents one plate, coloured by the number of transients on a logarithmic scale. The map is centred on the north celestial pole and keeps areas true, so equal patches of sky occupy equal space. Plates along the Milky Way hold fewer transients, as do plates close to a patch near Taurus and Perseus. Plates outlined with a dashed line are in the catalogue but were left out of VASCO’s published tests; we include them in our analysis.

### dates-h2
When was each plate taken?

### dates-intro
In the sky map above, {n:cat.left_out} plates are outlined with a dashed line. These are plates that were left out of VASCO’s tests, due to a problem with their dates.

The date of each plate is usually taken from the digitised archive at the Space Telescope Science Institute. This archive stores dates in a short block of text called a header, which are attached to each scan.

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
Asteroids move relative to stars, typically by around half an arcsecond a minute. During long image exposures, asteroids thus appear as short trails, while stars appear fixed since the telescope tracks them. Using NASA’s JPL Horizons service, we can determine the position of every known asteroid on any night in the 1950s with high precision.

For two example plates below, we show the trail a known asteroid would follow on the archive’s date (red dashed line) and where it would be a day later (red solid line). In one day, an asteroid moves around a quarter of a degree, so the two positions are well separated on the plate.

Using the archive’s date, we find no trail where the asteroid is predicted to be. One day later, however, we find trails left by asteroids at the expected position. Across the plates we checked, {n:dates.ast_late_trails} asteroid trails across {n:dates.ast_late_plates} plates appear one day after the archive date. On another {n:dates.ast_ok_plates} plates, {n:dates.ast_ok_trails} trails confirm that the archive date is correct, so the method is able to distinguish cases.

### dates-asteroid-1
Plate XE356, showing asteroid (210) Isabella. The archive date for this plate is 29 November 1951, but the asteroid trail actually appears on 30 November.

### dates-asteroid-2
Plate XE617, showing asteroid (163) Erigone. The archive date for this plate is 8 April 1956, but the asteroid trail appears on 9 April.

### dates-fix-title
Correcting the dates

### dates-fix-text
Palomar numbered the red plates in the order they were taken. The one-day error is not scattered at random: it affects long, unbroken runs of plate numbers, listed in the first table below. Palomar’s own log of the survey also survives. It records each night by its evening date, so the correct date is the log’s date plus one.

Where the log, the asteroids, and the handwritten labels cover the same plates, they agree on {n:dates.log_agree} of {n:dates.log_checked} cases. The two exceptions are typing errors in the log, which we correct by hand. For every other plate in a run with a one-day error, we add one day to the archive date.

Across the {n:dates.dated} plates we dated (which include all {n:cat.plates} in the catalogue), Palomar’s logs date {n:dates.by_log}, asteroid trails {n:dates.by_asteroids}, and the handwritten labels another {n:dates.by_label}. Runs of plate numbers settle a further {n:dates.by_runs}. That leaves {n:dates.unresolved} plates, XE723 and XE884, with no evidence for their date; neither is in the catalogue, so they play no part in the tests.

{n:dates.exceptions} plates don't follow this pattern, and are listed in the second table below. {n:dates.gross} of these plates are wrong by months or years, and we correct each of them by hand. All of the date corrections can be found in our code repository.

### dates-why
The mix of date conventions also explains most of the outlined plates on the sky map above. VASCO’s tests use only plates taken on one of {n:dates.vasco_nights} nights listed in Palomar’s log. The log names each night by its evening date, so a plate taken on the night of 8–9 April is listed under 8 April.

A plate whose archive date is one day early, 8 April, matches the list. A plate with the correct archive date, 9 April, matches only if Palomar also observed on the night of 9–10 April. If Palomar did not observe on the following night, the plate is dropped, which happened to {n:dates.dropped} of the {n:cat.left_out} outlined plates. If Palomar did observe that night, the plate is kept but counted as part of the following night, which happened to {n:dates.shifted} plates.

With the dates corrected, every one of the catalogue’s {n:cat.plates} plates falls on a night in Palomar’s log, so our analysis uses them all.

### cat-h2
Number of transients varies greatly across plates

### cat-intro
The {n:cat.plates} plates in the catalogue show a large degree of variation in how many transients they exhibit, ranging from {n:cat.min} to {n:cat.max}.

### cat-spread-title
More variation than expected by random chance

### cat-spread-text
VASCO’s Earth shadow test (see below) scores its result by treating every transient as an independent event: each has the same small chance of falling inside the shadow, regardless of other transients.

If transients were independent of one another and spread uniformly over the sky, plates would differ only by counting noise. An average plate holds {n:cat.avg} transients, and the counting noise on that is about its square root ({n:cat.noise}), so almost every plate would hold between {n:cat.band_lo} and {n:cat.band_hi} transients.

However, only {n:cat.in_band} of the {n:cat.plates} plates fall within this range. The actual counts spread smoothly from a handful of transients to several hundred: {n:cat.under30} plates hold fewer than 30 transients, while {n:cat.over500} hold more than 500. Some of this variation follows from where each plate points on the sky, as the sky map above shows.

How many transients a plate holds depends on the plate itself, so transients cannot be treated as independent events. A count drawn from only a few plates, such as the {n:sh.plates} crossed by Earth’s shadow, can therefore vary much more than counting noise suggests.

### cat-spread-caption
The blue histogram shows the distribution of transient numbers across the {n:cat.plates} plates. The final bar contains the {n:cat.over1000} plates with more than 1,000 transients.

The orange shaded band shows the range expected if transients were spread uniformly and independently. In this case, each plate would exhibit a number of transients close to the average ({n:cat.avg}), varying only by counting noise.

### cat-edge-title
More transients at plate edges

### cat-edge-text
Across the middle of a plate, transients are spread evenly. Towards the edges they grow denser, reaching up to about two and a half times the central density. Random positions that the VASCO team generated as a check show no such rise.

The edges of a plate are set by where the telescope happened to point, and the sky has no reason to be busier there. The rise most likely comes from the plates themselves, or from how they were copied, scanned, and matched. Near the edge, star images blur and can fail to match a modern catalogue.

### cat-edge-map-caption
Density of transients across the face of a plate, averaged over the {n:cat.plates} plates with each plate counted equally, and shown relative to the density near the centre (north up, east left). Left: transients. Right: the random positions VASCO generated as a check, which cover the plate evenly up to its edge.

### cat-edge-caption
Density of transients at each distance from the plate centre, shown relative to the centre density. Blue: transients, with shading showing how far the curve moves when the set of plates is resampled at random (90% of resamples). Grey: random positions. Of the {n:cat.plates} plates, {n:cat.edge_share}% show a higher transient density at the edge of the plate compared to the middle.

### cat-crowd-title
Fewer transients near the Milky Way

### cat-crowd-text
A point of light is kept as a transient only if no modern catalogue lists anything within 5 arcseconds. Near the Milky Way, stars are much more crowded, causing a greater proportion of transients to fail the 5 arcsecond threshold.

A typical plate within 5 degrees of the Milky Way holds {n:crowd.inner} transients. Counts rise steadily away from the galactic plane, and beyond about 40 degrees a typical plate holds around {n:crowd.outer}, {n:crowd.ratio} times as many.

Dust between the stars is also linked to fewer transients. Clouds of dust in our own Galaxy dim everything behind them, and the sparse patch near Taurus and Perseus on the sky map lies behind some of the nearest of these clouds.

To separate dust from crowding, we compared plates more than 20 degrees from the Milky Way, where stars are not crowded. The {n:dust.n_dusty} plates that look through the most dust hold a median of {n:dust.median_dusty} transients, around half the {n:dust.median_clear} held by the {n:dust.n_clear} plates that look through the least. Dust does not cover the whole of the Taurus–Perseus patch, however, and we have not found an explanation for the shortfall across the rest of it.

### cat-crowd-caption
Average number of transients per plate, plotted against distance from the plane of the Milky Way. At each distance, the line marks the median of all plates within 6 degrees of it, and the shading spans the middle 50% of plates.

### shadow-h2
Do transients avoid Earth’s shadow?

### shadow-intro
An object in orbit is seen only by the sunlight it reflects. When it passes into Earth’s shadow, it should therefore go dark.

VASCO considered objects in geostationary orbit (36,000 kilometres up), where a single orbit takes one day and a satellite stays above a spot on the equator. At this height, Earth’s shadow covers a disc of sky around 17 degrees across, centred on a point directly opposite the Sun. If some transients are due to glints from objects in orbit, the shadow disc should exhibit fewer transients than the same area of sky elsewhere.

Flaws in the photographic emulsion, dust on the glass, and marks from copying have no connection with where Earth’s shadow fell. A shortfall inside the disc, beyond what the differences between plates can produce, would therefore point to something in orbit. VASCO found such a shortfall, and we find it too. Where we differ is in the statistical interpretation of the deficit.

### shadow-step1-title
Nine plates intersect Earth's shadow

### shadow-step1-text
Across the {n:cat.plates} plates in the catalogue, Earth's shadow happened to intersect with {n:sh.plates_n} of them. Brian Doherty, a co-author of the 2026 paper, found {n:sh.vasco_plates} among the plates in VASCO's tests, as do we; the tenth, XE472, was taken on the night of 18–19 November 1949, two nights before VASCO's list of nights begins.

Throughout this section, we use the central 6 by 6 degrees of each plate, the plate size VASCO use in their own analysis. VASCO’s random positions cover this square evenly on every plate.

Using the corrected dates we describe above, the areas of those plates falling in Earth's shadow hold a total of {n:sh.observed} transients.

Together, the shadowed areas cover {n:sh.worth} plates’ worth of sky (all of plate XE422, two-thirds of two others, nearly half of XE472, and smaller slices of the remaining six). An average plate holds {n:sh.avg} transients in this square, so if transients were spread evenly across the survey, the shadowed parts would be expected to hold about {n:sh.expected} transients ({n:sh.worth} × {n:sh.avg}).

With the archive dates VASCO used, on the {n:sh.vasco_plates} plates in their tests, the shadow falls slightly differently, and the count is {n:sh.vasco_obs} against {n:sh.vasco_exp} expected: a similar shortfall of about a third.

### shadow-step1-caption
Of the {n:cat.plates} plates in the catalogue, {n:sh.plates_n} intersect with Earth's shadow. This includes one plate (XE422) which falls fully within Earth's shadow.

### shadow-step2-title
Measuring the probability of this deficit

### shadow-step2-text
VASCO scores this shortfall with a binomial test, which treats every transient as an independent event. Their headline figure of 22 standard deviations compared the shadow with its share of the whole sky; their 7.6 standard deviations used each plate’s actual coverage, as we do here. With that model and our corrected dates, the expected {n:sh.expected} would vary by about its square root ({n:sh.binom_sd}), which puts the observed {n:sh.observed} about {n:sh.binom_z} standard deviations below expectation.

Transients are not independent, however: as shown above, transient counts on plates differ far more than counting noise permits. The count inside a given shadow therefore depends on whether the plate it falls on happens to be rich or poor in transients. To measure how far that count can move, we took the {n:sh.plates} Earth shadow patches, each with its real shape and position, and laid them on {n:sh.plates} plates chosen at random.

A randomly chosen plate has no connection to Earth's shadow, so any shortfall in such trials comes from chance alone. Across 20,000 such trials, a shadow transient count of {n:sh.observed} or fewer turned up {n:sh.p_random}% of the time. This corresponds to a statistical significance of {n:sh.sig_random} standard deviations.

Where the shadow fell on plates also matters. {n:sh.taurus_cap} of the {n:sh.plates} shadow plates lie towards Taurus and Perseus, where the catalogue holds unusually few transients, so an average plate overstates what they should hold. Wesley Watters and colleagues made the same point, and found with their own random-shadow test, on a smaller public part of the catalogue, that a deficit this large arises by chance about 30% of the time. Drawing our random plates only from those near each shadow plate on the sky gives a probability of around one in {n:sh.odds_near} ({n:sh.sig_near} standard deviations).

VASCO also repeated their test using only transients within 2 degrees of each plate’s centre, away from the dense edges. In that region we count {n:sh.centre_obs} transients in shadow against {n:sh.centre_exp} expected, a shortfall random plates match about one time in {n:sh.odds_centre}.

Seen from Palomar rather than from the centre of the Earth, the shadow at geostationary height shifts by several degrees and crosses {n:sh.topo_plates} plates instead of {n:sh.plates}. This version of our analysis gives {n:sh.topo_obs} transients against {n:sh.topo_exp} expected, a shortfall random plates match about one time in {n:sh.odds_topo} ({n:sh.sig_topo} standard deviations).

### shadow-step2-caption
Orange: how the count would vary if every transient were independent (an average of {n:sh.expected}, with a standard deviation of {n:sh.binom_sd}). Blue: how the count varies when the real Earth shadow patches are laid on random plates; this represents the data-driven null distribution. The observed count of {n:sh.observed} is well below the orange distribution. However, when using the data-driven null distribution, no significant deviation is found.

### shadow-plates-title
Taking a closer look at the plates in shadow

### shadow-plates-text
Earth’s shadow crosses {n:sh.plates} of the catalogue’s plates. One of them, XE422, lies wholly inside Earth's shadow, so it has no sunlit area to compare with. The other {n:sh.partial} plates are shown below, with the position of every transient from VASCO’s catalogue indicated.

The plates differ significantly in terms of total transient counts. XE415 holds {n:pl.XE415_n} transients, while XE355 only holds {n:pl.XE355_n}. Within a plate, transients bunch into clumps and gather along some edges, while some plates are almost empty in places. For example, XE640’s transients crowd along its northern edge and are missing from two large rounded regions. The sky has no reason to be empty in two round patches of a single plate, so these gaps most likely come from the plate itself, or from how it was scanned and searched.

Earth’s shadow splits each plate in two. On XE617 the shadow covers {n:pl.XE617_area}% of the plate, and holds {n:pl.XE617_O} of its {n:pl.XE617_n} transients. On XE640 it clips {n:pl.XE640_area}% of the plate, in a corner that holds none of its {n:pl.XE640_n} transients.

### shadow-plates-caption
The {n:sh.partial} plates that Earth’s shadow partly covers, showing the central 6 by 6 degrees of each. Each dot is a transient in VASCO’s catalogue. Transients in Earth's shadow are shown in red, and those in sunlight blue. The red line marks the edge of Earth's shadow, 8.5 degrees from the point opposite the Sun when the plate was taken. Under each plate, “expected” shows how many transients the shadowed region would hold if transients there were as dense as in the sunlit region, allowing for the higher density towards plate edges. The ratio is the observed number of shadow transients divided by expected - a ratio of 1 means no difference in transient density between the shadow and sunlit regions.

### shadow-fit-title
Fitting the shadow transient density

### shadow-fit-text
The claim is that transients are rarer in Earth’s shadow. Expressed as a measurement: in the shadowed part of a plate, transients are r times as dense as in the sunlit part. With no Earth's shadow effect, the ratio r should be 1. If some transients genuinely vanish in the shadow, r should be below 1.

In practice, r compares two shares on each plate: the share of its transients that lie in the shadowed part, and the share the shadowed part would hold with no effect. Because transients are denser towards the edges and corners of every plate (see the plate maps above), this second share is not simply the shadowed area. We weight each part of the plate by the average density of transients there, measured only on the {n:sh.never} plates that Earth’s shadow never crossed, so a real shortfall in the shadow cannot affect the weights.

{n:fit.ex_plate}, for example, has {n:fit.ex_area}% of its area in shadow, which with no effect would hold {n:fit.ex_weighted}% of its transients. It actually holds {n:fit.ex_share}% ({n:fit.ex_O} of {n:fit.ex_n}). If transients vanish in the shadow, their share should fall below the expected share, and the ratio r measures by how much.

Because each plate is compared only with itself, it does not matter how many transients a plate holds in total. The variation between plates, large as it is, plays no part here.

The single value of r that best fits all {n:sh.partial} plates together is {n:fit.r}: transients are as dense in shadow as in sunlight. If every transient were an independent event, counting statistics would put r between {n:fit.count_lo} and {n:fit.count_hi}.

The plates themselves show that this uncertainty is too optimistic: single-plate ratios run from 0 on XE640 to {n:fit.max_r} on {n:fit.max_plate}, a much wider spread than counting allows. In {5_circle_symbol} we measure the real uncertainty, using trials built from the catalogue’s other plates.

### shadow-fit-caption
Grey points: the density of transients in shadow divided by the density in sunlight, for each plate that intersects Earth's shadow. Error bars are taken from counting statistics. Plates with no transients in shadow have only upper limits. Blue: the best-fit ratio across all {n:sh.partial} plates, with a one-standard-deviation range taken from counting statistics.

### shadow-toys-title
Finding the real uncertainty

### shadow-toys-text
Counting statistics treat every transient as independent. To measure the real uncertainty, we ran the same fit on plates where no shadow effect is possible. Each trial lays the {n:sh.partial} real shadow outlines onto {n:sh.partial} other plates that Earth’s shadow never crossed, and counts the number of transients inside and outside the shadow. In these trials, random plates are selected to contain a similar number of transients to the corresponding plates in the real analysis.

Over 20,000 such trials, the fit gives a typical r value of {n:fit.null_median}, with the middle 68% of results falling between {n:fit.null_lo} and {n:fit.null_hi}. Our best fit result on the {n:sh.partial} genuine plates (r = {n:fit.r}) falls well within this range, indicating that there is no evidence of reduced transient signals in Earth's shadow beyond random chance.

The same trials show that the fit can find a real shortfall. When we remove transients from inside the shadow outlines on the random plates, so that they are only {n:fit.inj_r} times as dense there (a shortfall as large as VASCO’s), the fit returns a value as high as {n:fit.r} only about one time in {n:fit.inj_odds}.

These results hold when the random plates are selected differently, when each shadow is instead moved around on its own plate, or when any one plate is left out.

### shadow-toys-caption
The best-fit ratio in 20,000 random trials, using plates that Earth’s shadow never crossed (grey). The best-fit result on the actual plates is shown in blue. Given that the best-fit result falls within the range of results expected purely by chance, there is no evidence of a deficit of transients in Earth's shadow.

### shadow-verdict-title
No significant shadow deficit

### shadow-verdict
Judged against a data-driven null distribution rather than counting noise, the pooled shortfall falls from {n:sh.binom_z} standard deviations to {n:sh.sig_random}. Within each plate, transients are also found to be no rarer in shadow than in sunlight: the best-fit ratio is {n:fit.r}, and at 90% confidence level the true ratio lies between {n:fit.ci90_lo} and {n:fit.ci90_hi}.

Earth’s shadow also crosses {n:sh.south} plates south of the celestial equator, beyond the edge of VASCO’s catalogue. Adding these plates would allow the same test to be performed with around {n:sh.south_ratio} times as much shadowed sky.

### nuclear-h2
Do transients cluster around nuclear tests?

### nuclear-intro
During the sky survey, the United States tested nuclear weapons in Nevada and the Pacific. VASCO’s two papers (Scientific Reports, 2025 and 2026) report that observing nights within a day of a US test held more transients than other nights, at a chance of about one in forty (p = 0.024) in the 2026 paper. The authors fixed this window before looking at the data, and read the sharpness of the association in time as a sign that it is genuine.

We find the same excess. Most of it, however, comes from the photographs themselves: the nights near tests held more photographs, pointed at parts of the sky that yield more transients. Once the comparison allows for this, the remaining excess is not statistically significant.

### nuc-vasco-title
VASCO’s test, as they ran it

### nuc-vasco-text
VASCO’s test counts the transients on each observing night, compares the nights within a day of a US test with all other nights, and asks how often the same number of nights picked at random do as well. We ran it counting each transient once, since VASCO’s machine-learning scores are not public.

On VASCO’s own construction, with the archive dates and their list of nights, the {n:nuc.vasco_nights} nights near tests hold {n:nuc.vasco_ratio} times the usual number of transients, and random nights do as well about one time in {n:nuc.vasco_odds} ({n:nuc.vasco_sig} standard deviations), close to VASCO’s one in forty. With every plate in the catalogue on its corrected date, {n:nuc.nights} nights fall within a day of a test. They hold {n:nuc.ratio} times the usual number, and random nights do as well about one time in {n:nuc.odds} ({n:nuc.sig} standard deviations).

### nuc-vasco-caption
The excess of transients on the {n:nuc.nights} nights within a day of a US test (red line), against the same number of nights drawn at random from the survey, 20,000 times (grey). Every catalogue plate, on corrected dates.

### nuc-corr-title
Counting photographs, and where they pointed

### nuc-corr-text
VASCO’s test counts transients per night, but a night’s transients come from its photographs. The {n:nuc.nights} nights near nuclear tests averaged {n:nuc.plates_window} plates a night, against {n:nuc.plates_other} on other nights. Only {n:nc2.empty_near} of them have no plates in the catalogue, against {n:nc2.empty_other} of the other {n:nc2.nights_other}; a night with no plates counts as zero transients, although nothing was observed. Counted per plate rather than per night, the excess falls from {n:nc2.night_r} to {n:nc2.photo_r}.

The sky above Palomar changes through the year as Earth orbits the Sun. On spring nights, the telescope faced away from the Milky Way, where there are fewer stars to rule transients out. Plates within 10 degrees of the Milky Way hold a median of {n:nc2.med_low} transients, against {n:nc2.med_high} more than 40 degrees away. Of the plates near nuclear tests, {n:nc2.spring}% were taken between February and June, which on its own would give them about {n:nc2.exp_pct}% more transients.

A like-for-like comparison is the number of transients each plate should hold, given where the telescope pointed. A curve fitted to all {n:cat.plates} plates gives every plate an expected count from its distance from the Milky Way, with the dust in front of it and its distance from the ecliptic. Measured against expected counts, plates near nuclear tests hold {n:nc2.point_r} times their expected number, rather than {n:nc2.night_r}.

Random nights drawn from across the survey are not a like-for-like comparison either. The tests came in series, so the nights near them bunch into a few seasons, and a bunched set of nights varies more than a scattered one. Sliding the real test calendar by 1 to 5 years and up to two months keeps each test’s time of year and its spacing from the others. Slid calendars do as well as the real one about one time in {n:nc2.odds} ({n:nc2.sig} standard deviations), which is not statistically significant.

### nuc-corr-caption
Once the number of photographs and the telescope pointing are accounted for, the plates near nuclear tests hold {n:nc2.point_r} times their expected number of transients. Of {n:nc2.slides} slid test calendars, about one in {n:nc2.odds} does as well or better.

### nuc-checks-title
Nothing else moves the result

### nuc-checks-text
We also allowed for every other property we can measure on each plate: how deep it goes (the share of faint Gaia stars it records), how sharp its stars are, exposure time, height above the horizon, the Moon, and VASCO’s plate-quality flag. None moves the result by more than a few per cent; with all of them together, it is {n:nc3.all}.

Nor does what counts as a transient. Counting 2-arcminute cells rather than transients, keeping only the central part of each plate, or keeping only transients with no other within 5 arcminutes all give similar results. Neither does leaving out any one series of tests. Every version lies between {n:nc3.lo} and {n:nc3.hi} times, and slid calendars do as well between {n:nc3.p_lo}% and {n:nc3.p_hi}% of the time.

### nuc-checks-caption
How far each further check moves the ratio for plates near nuclear tests, from {n:nc2.point_r} times per photograph allowing for pointing. Each version is measured against its own slid test calendars.

### nuc-window-title
Is VASCO’s window special?

### nuc-window-text
VASCO fixed their window, a test’s date and a day either side, before looking at the data, so its chance needs no allowance for other windows. The 2026 paper also reads the sharpness of the association in time as a sign that it is genuine. That claim compares windows, so we scored every window of consecutive days within two weeks of a test in the same way: per photograph, against pointing, and against the same slid calendars.

VASCO’s window ranks {n:nc4.rank} of {n:nc4.windows}, so about half of all windows are more striking. The most striking are single days a week or two from a test, led by {n:nc4.best} at {n:nc4.best_r} times, which slid calendars match {n:nc4.best_p}% of the time. Slid calendars produce a best window at least that striking {n:nc4.g14}% of the time, so none of the windows stands out.

### nuc-window-caption
Every window of consecutive days within two weeks of a US test that holds at least 5 plates, by its first and last day, coloured by how rarely slid test calendars match it. VASCO’s window is outlined.

### align-h2
Do transients form geometric patterns?

### align-intro
A flat object tumbling in orbit would flash each time it caught the Sun, leaving a line of points during a long plate exposure. VASCO searched for such alignments and shortlisted five candidates. Their best candidate was quoted at 3.9 standard deviations, corresponding to an estimated 1 in 10,000 chance of occurring by chance.

A catalogue of 100,000 unevenly spread points could contain many chance alignments. Three questions naturally arise:

- Do the alignment quality checks used introduce biases in the directions of the alignments they can find?
- Does the real transient catalogue contain more alignments than expected by random chance?
- What do the strongest alignments look like?

### align-step1-title
Searching for alignments in all directions

### align-step1-text
VASCO kept a set of transients only if its two sky coordinates, right ascension (east–west) and declination (north–south), rise and fall together almost perfectly: a correlation coefficient (Pearson’s r) above 0.99. This approach works for diagonal lines. However, along an east–west line, declination does not change, so there is nothing to correlate, and the line is rejected; the same happens to north–south lines.

In simulations (image below), the VASCO method keeps 99% of diagonal lines, 24% of lines 2 degrees from east–west, and only 1% of lines running exactly east–west. A geostationary satellite stays fixed in the sky above the observatory while the telescope follows the stars, so on the plate it moves exactly east–west. Objects at the same height on slightly tilted orbits move close to east–west too. VASCO’s selection rejects almost all of the tracks such an object would leave.

The same alignment task arises in particle physics, where the paths of charged particles must be reconstructed from hits they leave in a detector. The simplest tracking method fits a straight line through a set of points that could form a track, resulting in a track candidate. The quality of each candidate can be assessed, by checking how far the points lie from the line of best fit. Such a method can be applied in all spatial directions without losing signal, including along the east-west direction in this case.

### align-step1-caption
The share of simulated straight lines kept by each method, shown as a function of direction in the sky. VASCO’s selection (orange) collapses along east–west, where genuine geostationary satellite flashes would lie, as well as along north–south. A straight-line fit (blue) keeps nearly all lines in every direction.

### align-step2-title
Fitting every track candidate

### align-step2-text
VASCO’s five shortlisted alignments came from an earlier list of 298,165 candidate transients (Solano et al. 2022). Only {n:align.vasco_in_cat} of their {n:align.vasco_transients} transients appears in the later, cleaner catalogue of {n:cat.transients} transients that VASCO uses for its other tests. We have applied our track-fit method to the full catalogue of {n:cat.transients} transients.

We treated every set of three to six transients up to 15 arcminutes long as a track candidate, fitted a straight line through it, and recorded how closely its points follow the line. Two sources closer than about 4 arcseconds (the width of a star’s image on these scans) blur into one, so we ignored candidates with points closer than 8 arcseconds. The released catalogue itself holds almost no transients closer together than that. Three of VASCO’s five shortlisted alignments contain points 1 to 5 arcseconds apart, so they could not appear in the catalogue, regardless of the method used to search for alignments.

The position of each transient on a scan is uncertain by about an arcsecond. As such, even points from a genuine straight track can stray from the line by about that much. Each track fit is scored with a χ² per degree of freedom (χ²/ndof), which compares the points’ actual distances from the line with a one-arcsecond uncertainty. A genuine track gives a value near 1 or below. Points scattered loosely about the line give a much larger value.

A good fit alone does not make a track stand out: in crowded parts of a plate, pure chance can cause transients to align frequently. We therefore give each track candidate a track score: the number of tracks of the same length that would arise by chance in a 10 arcminute box around the candidate, given how many transients fall in the box.

This concept, and the formula we use, aligns with the method VASCO used to judge its own candidates (Edmunds and George, 1985). We apply our method in all directions, measure straightness from the fitted line, and score only track candidates that are straight to within the measurement error (a χ² probability of at least 1%). The lower the score, the less often chance alone can produce such a track.

A track candidate built from six transients also contains all shorter possible track candidates inside it. Where track candidates share a transient, we keep only the best-scoring one, so that each transient counts in at most one track candidate of each size.

### align-step3-title
Comparing track scores to random chance

### align-step3-text
We repeated the track score process for {n:align.copies} copies of the catalogue. In each copy, every transient is moved at random by up to {n:align.shift} arcseconds. Transients within {n:align.shift} arcseconds of one another are moved together, turned and shifted as one piece. This breaks genuine lines but retains the bunching of transients, down to the tightest clumps. Such copies enable us to measure how many track candidates are found purely due to chance.

Genuine glints in the real data coming from objects in orbit would show up as an excess of track candidates with low track scores, which chance alone can rarely produce. Across most of the track score distributions, the real catalogue matches the random data for all sizes of track candidate (see chart below).

The strongest-scoring tracks, however, sit somewhat above what chance alone would imply. {n:align.n5} track candidates are found with scores of 1 in 100,000 or better, against {n:align.exp5} expected by chance. {n:align.n4_cap} score 1 in 10,000 or better, against {n:align.exp4} expected due to chance. Allowing for the four sizes and two score levels we looked at, a difference at least this unusual turns up somewhere in around {n:align.p_global}% of random datasets.

It is important to assess such candidates in detail, to understand their physical characteristics in the plate imagery - we perform such an assessment in {5_circle_symbol}.

### align-step3-caption
Track score distributions for all 3-, 4-, 5-, and 6-track candidates in the catalogue. The real data distribution is shown in red, and the expected random background distribution is shown in grey. A small excess of track candidates is found in the low tails, beyond the rate expected purely by chance.

### align-vasco-title
Comparing our track score and VASCO’s score

### align-vasco-text
Our implementation of VASCO’s formula reproduces the scores in their published table of candidates. Where both methods find a track candidate, they agree closely (rank correlation of {n:align.rho}), and {n:align.top50} of the 50 top-scoring candidates found by the VASCO method are among the 50 ranked highest by the track score method. However, the methods diverge in two places.

The track score sets aside {n:align.vasco_bent} of the {n:align.vasco_strong} candidates VASCO’s score rates as 1 in 1,000 or better, because their transients do not fall along a straight line within the measurement error. The track score also rates candidates in all spatial directions, including those close to north-south and east-west which are rejected by VASCO’s Pearson correlation screening.

To measure the candidate-finding efficiency of each method, we added 100 artificial glint tracks to the catalogue in a set direction, and counted how many each method found. For a fair comparison, each score requires a pass mark that lets through the same number of false alarms, namely track candidates created purely by chance. We set both thresholds such that around {n:align.false_alarms} fake candidates pass across the whole survey, measured using copies of the catalogue where real alignments are broken (see {4_circle_symbol}).

For tracks running diagonally, both methods find about nine in ten alignments. Within two degrees of east–west, where a geostationary object’s track would lie, VASCO’s score finds fewer than half, and along east–west exactly it finds about {n:align.vasco_ew} in 100. In this same region, the track score still finds {n:align.track_ew}.

### align-vasco-caption
Left: VASCO’s score against the track score for the {n:align.both} track candidates in the real catalogue found by both methods. Darker cells hold more candidates, and the dashed line represents equal scores across both methods. Right: the proportion of artificial tracks found by each method, as a function of spatial direction. The higher efficiency of the track score method along the east-west and north-south directions is evident.

### align-tail-title
Taking a closer look at the highest-scoring tracks

### align-tail-text
The {n:align.n4} tracks with track scores of 1 in 10,000 or better are shown below. Chance alone would produce around {n:align.exp4} tracks with scores this strong, and less than one track with a score of 1 in 100,000 or better (the candidates on plates XE017 and XE492 both score in this range).

Using the toggle, you can switch between three photographs: the red plate the transients were found on, the blue plate taken on the same night, and the red plate from a second survey taken decades later. By definition, transients are missing from the later survey; all {n:align.n4} are also missing from the blue plates taken on the same night.

The strongest track, appearing on plate XE017, exhibits six star-like points along a line 15 arcminutes long. At high contrast in the red plate image, a pale band runs between the central four transients. The candidate on plate XE243 also exhibits similar pale banding structure. Light reaching the plate darkens it, so a pale band more likely points to the plate itself: damage to the emulsion, or something resting on it during the exposure.

On plates XE409 and XE493, the track candidates sit close to blots that exist only on the red plate, potentially caused by plate damage. The 5-transient track found on plate XE492 is the subject of further investigation in {7_circle_symbol}.

Four of the 3-transient track candidates exhibit a similar spatial structure: two points falling close together, with a third further along the line. We find that such structures are also found in the randomly varied data used to form the null distribution, and appear at around the same rate as in the real data.

A second, independent scan of the same plates, made in Edinburgh from a different set of copies, confirms that these features are also present on their scans. These features are thus present on the original plates, or on a master both sets of copies were created from, rather than appearing as damage to one copy.

### align-tail-findings
XE017: Pale band along the track | At high contrast, a pale band runs between the central transients on the red plate.

XE492: Curved dotted trail | The transients sit on a dotted trail that bends across the red plate and is absent from the blue.

XE493: Beside a blot | The points sit beside a large blot that appears only on the red plate.

XE294: Two close points and a third | Chance produces this layout about as often; other specks nearby appear only on the red plate.

XE409: Among blots | The points sit among blots that appear only on the red plate, one of them star-like.

XE516: Two close points and a third | Only loosely in line; other specks nearby appear only on the red plate.

XE243: Pale band along the track | At high contrast, a pale band joins the three points on the red plate.

### align-tail-caption


### align-widths-title
Comparing the transients with the stars around them

### align-widths-text
Every star on a plate was photographed through the same telescope during the same exposure. The air above the telescope blurs and shifts each star’s light, and a star that shines for the whole exposure records all of that movement. Stars of similar brightness on the same plate therefore share a typical shape, and averaging them provides a star model for that plate. Comparing each transient with a corresponding star model tells us something about how the image of that transient was produced.

Relative to a star, the image of a transient can either be wider, similar, or narrower:

- <em>Wider.</em> Light from anything farther than about {n:widths.defocus_km} kilometres, including an object in orbit, comes to the same focus as the stars, and cannot spread wider than them in every direction. A brief glint from a moving object would leave a short streak instead, whereas the wide transients observed appear round like the stars. A wider image is most likely due to something on the plate itself.
- <em>About as wide.</em> The light may have come through the telescope from a short-lived source, such as a glint from orbit, which would explain its absence from the blue plate and later surveys. But a speck on the plate the size of a star would look the same, and width alone cannot tell them apart.
- <em>Narrower.</em> A flash much shorter than the exposure records the air’s blurring for only an instant, and can be somewhat sharper than a star. A small, sharp flaw in the emulsion would also be narrow.

Measuring how much of each transient’s light lies away from its centre, compared with the star model, shows which of these scenarios applies. An image shaped like the stars gives a flat line at 1. A wider image has relatively more light away from its centre, and its line climbs above 1; a narrower one falls below. To judge each track, we compare it with groups of real stars from the same plate, matched to its transients in both number and brightness.

{n:widths.wider} tracks are flagged as wider than the profiles of similar stars: XE017 and XE243, which also feature pale bands around the transients, and XE409, which sits among blots. The other {n:widths.like} tracks feature star-like transients; none of the transients appear narrower than a typical star.

In general across the data, transients appear wider than typical stars. On average across {n:widths.control_n} transients selected at random, transients appear around {n:widths.control_pct}% wider than stars of similar brightness. Wide transient images are therefore common across the catalogue, which may point to them being features of the plates rather than light focused by the telescope. Hambly and Blair (2024) found the opposite: the nine transients from a single 1950 plate that they examined, on Edinburgh scans, had images sharper than stars. The two studies look at different samples. The catalogue we use also comes from a search that removed detections much sharper or broader than the stars around them (Solano et al. 2022), so it holds few sharp transients by construction. Both results, however, find transients unlike the stars around them. On plates from a different telescope, Hamburg Observatory’s 0.6-metre reflector, Busko (2026) examined 11 transients by eye and found that their images share the slight off-axis blur (coma) of the stars around them, which would point to light focused by a telescope; three of them were narrower than stars.

### align-widths-caption
For each of the {n:align.n4} strongest track candidates, we show the component transients alongside a star model, built from stars of similar brightness on the same plate. The graphs on the right show how each objects brightness falls away from its centre, relative to the star model; by definition, star-shaped images product a flat line at 1. The red line is the average over a track’s transients, with its 68% range shaded; the grey band shows the range for groups of real stars of similar brightness. Transients wider than stars produce a line that rises with distance from the centre. Tracks whose transients are star-like stay within the grey band. No track contains transients that appear narrower than the stars. The bottom row shows the same analysis for {n:widths.control_n} transients selected at random from the catalogue, where each is compared to stars of similar brightness from its own plate. On average, transient images are found to be around {n:widths.control_pct}% wider than stars.

### align-trail-title
The trail on plate XE492

### align-trail-text
The second-strongest track candidate, found on plate XE492, consists of five transients along a dotted trail. The transients appear similar in shape to stars on the same plate, and the trail is missing from the blue plate taken on the same night. The trail also appears on the Edinburgh copy of the same plate.

Extending the view across the plate, the trail runs for at least {n:trail.length} arcminutes (see image below). However, the full trail is not straight: along its path, it bends away from the track fit line by around {n:trail.sag} arcseconds. This corresponds to an arc with a radius of around {n:trail.radius} degrees.

Anything coasting through space, such as a meteor or an object in orbit, crosses a single plate in a straight line, so a bend this tight rules both out. An aircraft could turn this tightly, but it would be far too close to the telescope to be in focus. To be in focus, an object would need to be farther than about {n:widths.defocus_km} kilometres (see {6_circle_symbol}), and to trace this arc it would have to turn on a radius of at least {n:trail.turn_km} kilometres under its own power during the exposure.

The dots are also not evenly spaced. Their spacing runs from about {n:trail.spacing_lo} to {n:trail.spacing_hi} arcseconds, and differs from one part of the trail to the next. As such, they do not indicate a source flashing at a steady rate. Over its last 14 arcminutes, the trail also thickens into a band two to three dots wide, with dots sitting side by side across it. A single point of light crossing the sky would leave a line one dot wide.

A trail of this nature may be caused by an object coming into contact with the original plate, or the master that copies were printed from. For instance, a fibre lying against the emulsion would leave a sharp, gently curving line, and its twists could break the line into uneven dots as observed.

### align-trail-caption
Top: plate XE492 showing the 5 VASCO transients (red circles) and the straight line track fit (blue line).Bottom: the full extent of the trail, with an arc fitted to points extracted along the trail (red curve). The trail is consistent with an arc of {n:trail.radius} degrees curvature.

### align-verdict-title
Strong track candidates found, but context is important

### align-verdict
The track score we propose builds on VASCO’s own measure but can see in every direction, including along the east–west direction where genuine geostationary objects would leave tracks.

Applied to the full catalogue of transients, our method finds track candidates at the rate expected by chance, with the exception of a handful of strong track candidates. These candidates appear only on the red plate, and not in the blue plate or subsequent survey imagery. Close inspection of their imagery can help understand the context of each candidate, for example whether plate damage or a foreign object is present.

### end-h2
Summary

### card-4-question
Are the plates correctly dated?

### card-4-ours
Around {n:dates.early_pct}% of plates have an archive date that is one day early. Asteroid trails, Palomar’s own logs, and the handwritten dates on plates independently corroborate this finding. We provide a correction to the dates in our code repository.

### card-1-question
Do transients avoid Earth’s shadow?

### card-1-ours
Plates differ so much that the shortfall happens by chance about one time in {n:sh.odds_random}. Within each plate, allowing for the higher density at plate edges, transients are as dense in shadow as in sunlight; a shortfall as large as VASCO’s would give this result only about one time in {n:fit.inj_odds}.

### card-2-question
Do transients cluster around nuclear tests?

### card-2-ours
The nights near tests held more photographs, pointed at parts of the sky that yield more transients. Counted per photograph, and allowing for where the telescope pointed, the excess is {n:nc2.point_r} times, which slid test calendars match about one time in {n:nc2.odds}.

### card-3-question
Do transients form linear patterns?

### card-3-ours
A track score that sees in every direction reveals a catalogue mostly consistent with chance alignments. Careful assessment of the plate imagery for strongest candidates can help provide additional context, for example the presence of nearby plate damage.

### end-text
There are three things we would like to see next:

- Public release of the machine-learning scores. The 2026 VASCO paper reports stronger results for the transients its model rates most likely to be genuine, but the score column in the released catalogue assigns the same value to 99.8% of transients. We would be glad to put the scores through the plate-by-plate tests above.
- Extension of the analyses to the plates south of the celestial equator. Earth’s shadow alone crosses {n:sh.south} of them, which would give around {n:sh.south_ratio} times as much shadowed sky data.
- Use of our corrected plate dates and track score in future analyses, to improve their accuracy and to find tracks in every direction.

We would welcome the VASCO team’s response to this analysis, and reiterate our gratitude for making their data and code available.

### footer-about
Built from the VASCO team’s public catalogue, VASCO’s analysis code in Brian Doherty’s VASCO-ML repository, and our own scripts. Asteroid positions come from JPL Horizons, star positions from Gaia, plate times from Palomar’s log (VizieR VI/25), and test times from the US Department of Energy report DOE/NV-209. Every result of our own analysis can be reproduced from our repository, which also holds the table of corrected plate dates; figures quoted from VASCO and from other authors are taken from their papers. We thank Beatriz Villarroel and her colleagues for publishing their data and code; their openness made this work possible.
