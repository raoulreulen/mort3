/*==============================================================================
FIGURE: LIFE EXPECTANCY BY YEAR OF CHILDHOOD CANCER DIAGNOSIS

This do-file creates a figure showing:
- Observed life expectancy by year of diagnosis (overall, all cancers combined)
- Expected life expectancy (under general population mortality) by year of diagnosis
- Uses stpm2 (flexible parametric survival models) for smooth estimates

==============================================================================*/

*-------------------------------------------------------------------------------
* LOAD STSET DATA FOR ALL-CAUSE MORTALITY
*-------------------------------------------------------------------------------
use "$temp/x-mort3-stset-allcauses", clear

*-------------------------------------------------------------------------------
* CREATE YEAR AND DECADE OF DIAGNOSIS VARIABLES
*-------------------------------------------------------------------------------
// Extract year from doe (date of entry = date of diagnosis)
gen year_diag = year(doe)

// Create decade categories for easier modeling
gen decade_diag = floor(year_diag/10)*10

// Check range
summ year_diag decade_diag, detail

// Label decade
label define ldecade 1940 "1940s" 1950 "1950s" 1960 "1960s" 1970 "1970s" ///
	1980 "1980s" 1990 "1990s" 2000 "2000s" 2010 "2010s"
label values decade_diag ldecade

*-------------------------------------------------------------------------------
* CALCULATE OBSERVED AND EXPECTED LE BY DECADE (EMPIRICAL FIRST)
*-------------------------------------------------------------------------------
// Calculate person-time and deaths
gen person_years = _t - _t0
gen expected_deaths = person_years * rate_allcauses

// Collapse by decade to get totals
preserve
collapse (sum) deaths=_d (sum) py=person_years (sum) exp_d=expected_deaths, ///
	by(decade_diag)

// Calculate SMR
gen smr = deaths / exp_d

// Simple LE approximation: assume constant hazard within decade
// Observed hazard rate
gen obs_rate = deaths / py

// Expected hazard rate
gen exp_rate = exp_d / py

// LE = 1/rate (under constant hazard assumption, truncated at 80)
gen le_observed_simple = min(1/obs_rate, 80)
gen le_expected_simple = min(1/exp_rate, 80)
gen le_lost_simple = le_expected_simple - le_observed_simple

list decade_diag deaths py smr le_observed_simple le_expected_simple le_lost_simple

tempfile decade_summary
save `decade_summary'
restore

*-------------------------------------------------------------------------------
* FIT STPM2 MODEL FOR SMOOTH ESTIMATES
*-------------------------------------------------------------------------------
// Center year for modeling
summ year_diag, meanonly
local mean_year = r(mean)
gen year_diag_c = year_diag - `mean_year'

di as text _newline "Fitting flexible parametric survival model..."
di as text "Model includes year of diagnosis with time-varying effects"
di as text "This may take a few minutes..."

// Fit stpm2 model
stpm2 year_diag_c, scale(hazard) df(4) tvc(year_diag_c) dftvc(2) nolog

// Store estimates
estimates store obs_model

*-------------------------------------------------------------------------------
* PREDICT LIFE EXPECTANCY CURVES BY YEAR OF DIAGNOSIS
*-------------------------------------------------------------------------------
// Get year range
summ year_diag, meanonly
local min_year = r(min)
local max_year = r(max)

// Determine how many years to predict for
local year_range = `max_year' - `min_year' + 1

// Create results dataset
preserve
clear
set obs `year_range'
gen year_diag = `min_year' + _n - 1
gen year_diag_c = year_diag - `mean_year'

// Restore estimates
estimates restore obs_model

// Create time variable for prediction
expand 200
bysort year_diag: gen _temptime = (_n - 1) * 80 / 199

// Generate the year_diag_c values for each time point
// (already exists in the data)

// Predict survival for each year-time combination
predict _tempsurv, survival timevar(_temptime) zeros

// Calculate RMST by integrating for each year
bysort year_diag (_temptime): gen _area = (_tempsurv + _tempsurv[_n-1])/2 * (_temptime - _temptime[_n-1]) if _n > 1
by year_diag: egen le_observed = total(_area)
by year_diag: keep if _n == 1

// Clean up
drop _temptime _tempsurv _area

di as text "Life expectancy prediction completed"

// Use the general population mortality rates to calculate expected LE
// Merge back with original data to get average background rates by year
tempfile pred_le
save `pred_le'
restore

// Calculate expected LE by integrating general population survival
// Group original data by year and calculate average rate at each age
preserve
gen age_year = floor((_t0 + _t) / 2)
collapse (mean) avg_bg_rate=rate_allcauses, by(year_diag age_year)

// For each year of diagnosis, calculate expected survival
gen exp_surv = .
bysort year_diag (age_year): replace exp_surv = exp(-sum(avg_bg_rate)) if age_year <= 80

// Integrate to get life expectancy
bysort year_diag (age_year): gen _area_exp = (exp_surv + exp_surv[_n-1])/2 * (age_year - age_year[_n-1]) if _n > 1 & age_year <= 80
by year_diag: egen le_expected = total(_area_exp)
by year_diag: keep if _n == 1

keep year_diag le_expected
tempfile expected_by_year
save `expected_by_year'
restore

// Merge with predictions
use `pred_le', clear
merge 1:1 year_diag using `expected_by_year', nogen

// Calculate life years lost
gen le_lost = le_expected - le_observed

// Save dataset
save "$temp/x-mort3-le-by-year-diagnosis", replace

// Display sample
list year_diag le_observed le_expected le_lost in 1/10

*-------------------------------------------------------------------------------
* CREATE SMOOTHED FIGURES
*-------------------------------------------------------------------------------
// Determine y-axis range
summ le_expected, meanonly
local ymax = ceil(r(max)/5)*5
summ le_observed, meanonly
local ymin = floor(r(min)/5)*5

// Sort for graphing
sort year_diag

// Create figure showing observed and expected life expectancy by year
twoway ///
	(line le_expected year_diag, ///
		lcolor(blue) lwidth(medium) lpattern(dash)) ///
	(line le_observed year_diag, ///
		lcolor(red) lwidth(medium) lpattern(solid)) ///
	, ///
	legend(order(1 "Expected LE (general population)" ///
		2 "Observed LE (cancer survivors)") ///
		position(6) ring(0) cols(1) size(small)) ///
	ytitle("Life Expectancy (years)", size(medium)) ///
	xtitle("Year of Childhood Cancer Diagnosis", size(medium)) ///
	ylabel(`ymin'(5)`ymax', angle(horizontal) format(%9.0f) labsize(small)) ///
	yscale(range(`ymin' `ymax')) ///
	xlabel(, labsize(small)) ///
	title("Life Expectancy by Year of Diagnosis", size(medium)) ///
	subtitle("All childhood cancers combined (stpm2 smoothed)", size(small)) ///
	graphregion(color(white)) plotregion(color(white)) ///
	scheme(s2color)

// Save graph
graph export "$temp/fig_le_by_year_diagnosis.png", replace width(2400) height(1800)
graph export "$temp/fig_le_by_year_diagnosis.pdf", replace
graph save "$temp/fig_le_by_year_diagnosis.gph", replace

*-------------------------------------------------------------------------------
* CREATE ALTERNATIVE FIGURE: LIFE YEARS LOST BY YEAR OF DIAGNOSIS
*-------------------------------------------------------------------------------
twoway ///
	(line le_lost year_diag, ///
		lcolor(red) lwidth(medium) lpattern(solid)) ///
	, ///
	ytitle("Life Years Lost (years)", size(medium)) ///
	xtitle("Year of Childhood Cancer Diagnosis", size(medium)) ///
	ylabel(, angle(horizontal) format(%9.1f) labsize(small)) ///
	xlabel(, labsize(small)) ///
	yline(0, lcolor(black) lpattern(solid)) ///
	title("Life Years Lost by Year of Diagnosis", size(medium)) ///
	subtitle("All childhood cancers combined (stpm2 smoothed)", size(small)) ///
	graphregion(color(white)) plotregion(color(white)) ///
	scheme(s2color)

// Save graph
graph export "$temp/fig_le_lost_by_year_diagnosis.png", replace width(2400) height(1800)
graph export "$temp/fig_le_lost_by_year_diagnosis.pdf", replace
graph save "$temp/fig_le_lost_by_year_diagnosis.gph", replace

*-------------------------------------------------------------------------------
* EXPORT TO WORD DOCUMENT
*-------------------------------------------------------------------------------
cap: putdocx clear
putdocx begin, font(arial, "10")

// Add first figure
putdocx paragraph, style(Heading1)
putdocx text ("Life Expectancy by Year of Childhood Cancer Diagnosis")
putdocx paragraph
putdocx image "$temp/fig_le_by_year_diagnosis.png", width(6.5)

putdocx pagebreak

// Add second figure
putdocx paragraph, style(Heading1)
putdocx text ("Life Years Lost by Year of Diagnosis")
putdocx paragraph
putdocx image "$temp/fig_le_lost_by_year_diagnosis.png", width(6.5)

putdocx save "$temp/figure_le_by_year_diagnosis.docx", replace

di as text ""
di as text "Figures exported to:"
di as text "  Word document: " as result "$temp/figure_le_by_year_diagnosis.docx"
di as text "  Individual PNGs and PDFs in: " as result "$temp/"
di as text "  Dataset saved: " as result "$temp/x-mort3-le-by-year-diagnosis.dta"

exit
