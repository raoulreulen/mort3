/*==============================================================================
 GRAPH: STANDARDIZED MORTALITY RATIO (SMR) BY ATTAINED AGE - SMOOTHED

 Purpose: Calculate SMR by attained age with finer age intervals and produce
          a smoothed graph using fractional polynomials

 Input:   Uses stset data from mort3-02a-rr-stsetdata.do
 Output:  Graph showing observed SMR points and smoothed curve with CIs
==============================================================================*/

*-------------------------------------------------------------------------------
* LOAD DATA FOR ALL-CAUSE MORTALITY
*-------------------------------------------------------------------------------
use "$temp/x-mort3-stset-allcauses", clear

*-------------------------------------------------------------------------------
* CALCULATE EXPECTED DEATHS AND OBSERVED DEATHS BY FINER AGE INTERVALS
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

// Calculate SMR and 95% CI
// SMR = Observed / Expected
gen smr = obs / exp

// 95% CI for SMR using exact Poisson method
// Lower CI: χ²(2*obs, 0.025) / (2*expected)
// Upper CI: χ²(2*(obs+1), 0.975) / (2*expected)
gen smr_lower = (invchi2(2*obs, 0.025) / 2) / exp if obs > 0
gen smr_upper = (invchi2(2*(obs+1), 0.975) / 2) / exp

// For obs=0, use rule of three: upper limit = -ln(0.05)/expected
replace smr_lower = 0 if obs == 0
replace smr_upper = -log(0.05) / exp if obs == 0

// Keep only age groups with sufficient follow-up
drop if exp < 0.5  // Exclude age groups with very small expected deaths

// Label variables
label variable age_group "Attained age (years)"
label variable smr "Standardized Mortality Ratio"
label variable smr_lower "SMR 95% CI Lower"
label variable smr_upper "SMR 95% CI Upper"
label variable obs "Observed deaths"
label variable exp "Expected deaths"

*-------------------------------------------------------------------------------
* FIT MULTIPLE SMOOTHING METHODS FOR COMPARISON
*-------------------------------------------------------------------------------

// METHOD 1: Fractional polynomial (degree 2)
di as text "Fitting fractional polynomial model..."
fracpoly, degree(2): regress smr age_group [aweight=exp]
predict smr_fp2, xb
predict smr_fp2_se, stdp
gen smr_fp2_lower = smr_fp2 - 1.96 * smr_fp2_se
gen smr_fp2_upper = smr_fp2 + 1.96 * smr_fp2_se
di as text "FP powers: " as result e(fp_k1) " and " e(fp_k2)
di as text "R-squared: " as result e(r2)

// METHOD 2: Restricted cubic spline (3 knots - less flexible)
di as text _newline "Fitting restricted cubic spline (3 knots)..."
mkspline age_rcs3 = age_group, cubic nknots(3)
regress smr age_rcs3* [aweight=exp]
predict smr_rcs3, xb
di as text "R-squared: " as result e(r2)
drop age_rcs3*

// METHOD 3: Restricted cubic spline (5 knots - more flexible)
di as text _newline "Fitting restricted cubic spline (5 knots)..."
mkspline age_rcs5 = age_group, cubic nknots(5)
regress smr age_rcs5* [aweight=exp]
predict smr_rcs5, xb
di as text "R-squared: " as result e(r2)
drop age_rcs5*

// METHOD 4: Lowess smoothing (bandwidth 0.5)
di as text _newline "Fitting lowess smooth..."
lowess smr age_group, bwidth(0.5) gen(smr_lowess) nograph

// Approximate confidence intervals using local standard deviation
// Calculate residuals from lowess fit
gen lowess_resid = smr - smr_lowess

// Estimate local variability using moving window
gen smr_lowess_se = .
levelsof age_group, local(ages)
foreach age of local ages {
	// Calculate SE based on nearby residuals (within ±5 years)
	quietly summarize lowess_resid if abs(age_group - `age') <= 5 & lowess_resid != ., detail
	if r(N) > 0 {
		local se = r(sd) / sqrt(r(N))
		quietly replace smr_lowess_se = `se' in `=_N'
		quietly replace smr_lowess_se = `se' if age_group == `age'
	}
}

// Create 95% confidence intervals
gen smr_lowess_lower = smr_lowess - 1.96 * smr_lowess_se
gen smr_lowess_upper = smr_lowess + 1.96 * smr_lowess_se

drop lowess_resid

// METHOD 5: Simple polynomial (degree 3)
di as text _newline "Fitting cubic polynomial..."
gen age_sq = age_group^2
gen age_cube = age_group^3
regress smr age_group age_sq age_cube [aweight=exp]
predict smr_poly3, xb
di as text "R-squared: " as result e(r2)

*-------------------------------------------------------------------------------
* CREATE COMPARISON GRAPH WITH ALL SMOOTHING METHODS
*-------------------------------------------------------------------------------
// Sort by age for proper line drawing
sort age_group

// Determine maximum y-axis value from data in display range
summ smr_fp2_upper if age_group >= 20 & age_group <= 75, meanonly
local ymax = ceil(r(max)/5)*5  // Round up to nearest 5

// Create graph comparing all smoothing methods
twoway ///
	(rarea smr_fp2_lower smr_fp2_upper age_group if age_group >= 20 & age_group <= 75, ///
		color(blue%15) lwidth(none)) ///
	(line smr_fp2 age_group if age_group >= 20 & age_group <= 75, ///
		lcolor(blue) lwidth(medium) lpattern(solid)) ///
	(line smr_rcs3 age_group if age_group >= 20 & age_group <= 75, ///
		lcolor(red) lwidth(medium) lpattern(dash)) ///
	(line smr_rcs5 age_group if age_group >= 20 & age_group <= 75, ///
		lcolor(green) lwidth(medium) lpattern(dash_dot)) ///
	(line smr_lowess age_group if age_group >= 20 & age_group <= 75, ///
		lcolor(orange) lwidth(medium) lpattern(shortdash)) ///
	(line smr_poly3 age_group if age_group >= 20 & age_group <= 75, ///
		lcolor(purple) lwidth(medium) lpattern(longdash)) ///
	(scatter smr age_group if age_group >= 20 & age_group <= 75, ///
		mcolor(gs8) msize(vsmall) msymbol(circle)), ///
	legend(order(2 "Fractional polynomial" 3 "Cubic spline (3 knots)" ///
		4 "Cubic spline (5 knots)" 5 "Lowess" 6 "Cubic polynomial" ///
		7 "Observed SMR" 1 "95% CI (FP)") ///
		position(11) ring(0) cols(1) size(*.55)) ///
	ytitle("Standardised Mortality Ratio (SMR)", size(medium)) ///
	xtitle("Attained age (years)", size(medium)) ///
	ylabel(0 1 5(5)`ymax', angle(horizontal) format(%9.0f) labsize(small)) ///
	yscale(range(0 `ymax')) ///
	xlabel(20(5)75, labsize(small)) ///
	xscale(range(20 75)) ///
	yline(1, lcolor(black) lpattern(dash)) ///
	title("Standardised Mortality Ratio by Attained Age", size(medium)) ///
	subtitle("Comparison of smoothing methods", size(small)) ///
	graphregion(color(white)) plotregion(color(white)) ///
	scheme(s2color)

// Save comparison graph
graph export "$temp/smr_by_age_smooth_comparison.png", replace width(2400) height(1800)
graph export "$temp/smr_by_age_smooth_comparison.pdf", replace
graph save "$temp/smr_by_age_smooth_comparison.gph", replace

*-------------------------------------------------------------------------------
* CREATE INDIVIDUAL GRAPH WITH LOWESS (PREFERRED)
*-------------------------------------------------------------------------------
// Create graph with lowess smoothing and confidence intervals
twoway ///
	(rarea smr_lowess_lower smr_lowess_upper age_group if age_group >= 20 & age_group <= 75, ///
		color(gs10) fintensity(50) lwidth(none)) ///
	(line smr_lowess age_group if age_group >= 20 & age_group <= 75, ///
		lcolor(blue) lwidth(medium) lpattern(solid)) ///
	(scatter smr age_group if age_group >= 20 & age_group <= 75, ///
		mcolor(red) msize(vsmall) msymbol(circle)), ///
	legend(order(2 "Smoothed SMR (lowess)" 1 "95% CI (smoothed)" ///
		3 "Observed SMR") ///
		position(1) ring(0) cols(1) size(*.5)) ///
	ytitle("Standardised Mortality Ratio (SMR)", size(medium)) ///
	xtitle("Attained age (years)", size(medium)) ///
	ylabel(0 1 5(5)`ymax', angle(horizontal) format(%9.0f) labsize(small)) ///
	yscale(range(0 `ymax')) ///
	xlabel(20(5)75, labsize(small)) ///
	xscale(range(20 75)) ///
	yline(1, lcolor(black) lpattern(dash)) ///
	title("Standardised Mortality Ratio by Attained Age", size(medium)) ///
	subtitle("All-cause mortality with lowess smoothing", size(small)) ///
	graphregion(color(white)) plotregion(color(white)) ///
	scheme(s2color)

// Save individual graph
graph export "$temp/smr_by_age_smooth.png", replace width(2400) height(1800)
graph export "$temp/smr_by_age_smooth.pdf", replace
graph save "$temp/smr_by_age_smooth.gph", replace

*-------------------------------------------------------------------------------
* EXPORT GRAPHS TO WORD AND PDF DOCUMENTS
*-------------------------------------------------------------------------------
// Create Word document with all three graphs
cap: putdocx clear
putdocx begin, font(arial, "10")

// Add comparison graph
putdocx paragraph, style(Heading1)
putdocx text ("Comparison of Smoothing Methods")
putdocx paragraph
putdocx image "$temp/smr_by_age_smooth_comparison.png", width(6.5)

putdocx pagebreak

// Add individual SMR graph
putdocx paragraph, style(Heading1)
putdocx text ("Standardised Mortality Ratio by Attained Age")
putdocx paragraph
putdocx image "$temp/smr_by_age_smooth.png", width(6.5)

putdocx save "$temp/figure_smr_by_age.docx", replace

// Export to PDF using graph combine
graph use "$temp/smr_by_age_smooth_comparison.gph", clear
graph export "$temp/smr_by_age_smooth_comparison_for_pdf.png", replace width(3200) height(2400)

graph use "$temp/smr_by_age_smooth.gph", clear
graph export "$temp/smr_by_age_smooth_for_pdf.png", replace width(3200) height(2400)

di as text ""
di as text "Graphs exported to:"
di as text "  Word document: " as result "$temp/figure_smr_by_age.docx"
di as text "  Individual PNGs and PDFs in: " as result "$temp/"

*-------------------------------------------------------------------------------
* SAVE DATASET WITH SMR BY AGE
*-------------------------------------------------------------------------------
save "$temp/x-mort3-smr-by-age-allcauses", replace

// Display summary statistics
list age_group obs exp smr smr_fp2 smr_rcs3 smr_rcs5 smr_lowess smr_poly3 in 1/20

restore

exit
