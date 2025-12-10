/*==============================================================================
 GRAPH: ABSOLUTE EXCESS RISK (AER) BY ATTAINED AGE - SMOOTHED

 Purpose: Calculate AER by attained age with finer age intervals and produce
          a smoothed graph using lowess

 Input:   Uses stset data from mort3-02a-rr-stsetdata.do
 Output:  Graph showing observed AER points and smoothed curve with CIs
==============================================================================*/

*-------------------------------------------------------------------------------
* LOAD DATA FOR ALL-CAUSE MORTALITY
*-------------------------------------------------------------------------------
use "$temp/x-mort3-stset-allcauses", clear

*-------------------------------------------------------------------------------
* CALCULATE AER BY FINER AGE INTERVALS
*-------------------------------------------------------------------------------
// Create finer age groups (e.g., single year or 2-year intervals)
// Using the midpoint of follow-up time (_t + _t0)/2 for attained age
gen attained_age = (_t + _t0)/2

// Round to create age groups (adjust interval as needed)
// Using 1-year intervals for finer granularity
gen age_group = floor(attained_age)

// Calculate person-years, observed deaths, and expected deaths by age group
gen person_years = _t - _t0
gen expected = person_years * rate_allcauses

// Collapse by age group to get totals
preserve
collapse (sum) obs=_d exp=expected py=person_years (mean) age_mid=attained_age, by(age_group)

// Calculate AER per 10,000 person-years
// AER = (Observed - Expected) / Person-years * 10,000
gen aer = ((obs - exp) / py) * 10000

// Calculate 95% CI for AER using Poisson-based approach
// Variance of (O-E)/PY = (O/PY^2) approximately
// SE of AER
gen aer_se = (sqrt(obs) / py) * 10000

// 95% CI
gen aer_lower = aer - 1.96 * aer_se
gen aer_upper = aer + 1.96 * aer_se

// Keep only age groups with sufficient follow-up
drop if exp < 0.5  // Exclude age groups with very small expected deaths

// Label variables
label variable age_group "Attained age (years)"
label variable aer "Absolute Excess Risk (per 10,000 PY)"
label variable aer_lower "AER 95% CI Lower"
label variable aer_upper "AER 95% CI Upper"
label variable obs "Observed deaths"
label variable exp "Expected deaths"

*-------------------------------------------------------------------------------
* FIT LOWESS SMOOTHING
*-------------------------------------------------------------------------------
// Sort by age for proper line drawing
sort age_group

di as text "Fitting lowess smooth..."
lowess aer age_group, bwidth(0.5) gen(aer_lowess) nograph

// Approximate confidence intervals using local standard deviation
// Calculate residuals from lowess fit
gen lowess_resid = aer - aer_lowess

// Estimate local variability using moving window
gen aer_lowess_se = .
levelsof age_group, local(ages)
foreach age of local ages {
	// Calculate SE based on nearby residuals (within ±5 years)
	quietly summarize lowess_resid if abs(age_group - `age') <= 5 & lowess_resid != ., detail
	if r(N) > 0 {
		local se = r(sd) / sqrt(r(N))
		quietly replace aer_lowess_se = `se' in `=_N'
		quietly replace aer_lowess_se = `se' if age_group == `age'
	}
}

// Create 95% confidence intervals
gen aer_lowess_lower = aer_lowess - 1.96 * aer_lowess_se
gen aer_lowess_upper = aer_lowess + 1.96 * aer_lowess_se

drop lowess_resid

*-------------------------------------------------------------------------------
* CREATE GRAPH
*-------------------------------------------------------------------------------
// Determine maximum and minimum y-axis values from data in display range
summ aer_lowess_upper if age_group >= 20 & age_group <= 70, meanonly
local ymax = ceil(r(max)/100)*100  // Round up to nearest 100

summ aer_lowess_lower if age_group >= 20 & age_group <= 70, meanonly
local ymin = floor(r(min)/100)*100  // Round down to nearest 100
if `ymin' > 0 {
	local ymin = 0  // Start at 0 if all values are positive
}

// Create graph with lowess smoothing and confidence intervals
twoway ///
	(rarea aer_lowess_lower aer_lowess_upper age_group if age_group >= 20 & age_group <= 70, ///
		color(gs10) fintensity(50) lwidth(none)) ///
	(line aer_lowess age_group if age_group >= 20 & age_group <= 70, ///
		lcolor(blue) lwidth(medium) lpattern(solid)) ///
	(scatter aer age_group if age_group >= 20 & age_group <= 70, ///
		mcolor(red) msize(vsmall) msymbol(circle)), ///
	legend(order(2 "Smoothed AER (lowess)" 1 "95% CI (smoothed)" ///
		3 "Observed AER") ///
		position(11) ring(0) cols(1) size(*.5)) ///
	ytitle("Absolute Excess Risk (per 10,000 person-years)", size(medium)) ///
	xtitle("Attained age (years)", size(medium)) ///
	ylabel(`ymin'(100)`ymax', angle(horizontal) format(%9.0f) labsize(small)) ///
	yscale(range(`ymin' `ymax')) ///
	xlabel(20(5)70, labsize(small)) ///
	xscale(range(20 70)) ///
	yline(0, lcolor(black) lpattern(dash)) ///
	title("Absolute Excess Risk by Attained Age", size(medium)) ///
	subtitle("All-cause mortality with lowess smoothing", size(small)) ///
	graphregion(color(white)) plotregion(color(white)) ///
	scheme(s2color)

// Save graph
graph export "$temp/aer_by_age_smooth.png", replace width(2400) height(1800)
graph export "$temp/aer_by_age_smooth.pdf", replace
graph save "$temp/aer_by_age_smooth.gph", replace

*-------------------------------------------------------------------------------
* EXPORT GRAPH TO WORD AND PDF DOCUMENTS
*-------------------------------------------------------------------------------
// Create Word document with AER graph
cap: putdocx clear
putdocx begin, font(arial, "10")

// Add AER graph
putdocx paragraph, style(Heading1)
putdocx text ("Absolute Excess Risk by Attained Age")
putdocx paragraph
putdocx image "$temp/aer_by_age_smooth.png", width(6.5)

putdocx save "$temp/figure_aer_by_age.docx", replace

// Export high-res version for PDF
graph use "$temp/aer_by_age_smooth.gph", clear
graph export "$temp/aer_by_age_smooth_for_pdf.png", replace width(3200) height(2400)

di as text ""
di as text "Graph exported to:"
di as text "  Word document: " as result "$temp/figure_aer_by_age.docx"
di as text "  PNG and PDF in: " as result "$temp/"

*-------------------------------------------------------------------------------
* SAVE DATASET WITH AER BY AGE
*-------------------------------------------------------------------------------
save "$temp/x-mort3-aer-by-age-allcauses", replace

// Display summary statistics
list age_group obs exp py aer aer_lower aer_upper aer_lowess in 1/20

restore

exit
