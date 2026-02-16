/*==============================================================================
TABLE 3: ALL-CAUSE SMRs AND AERs BY CHILDHOOD CANCER TYPE (DIAG)

This do-file creates a Word table showing:
- All-cause SMR and AER by diagnosis (cancer type)
- Observed life expectancy for each cancer type
- Expected life expectancy (under general population mortality)

==============================================================================*/

*-------------------------------------------------------------------------------
* LOAD STSET DATA FOR ALL-CAUSE MORTALITY
*-------------------------------------------------------------------------------
use "$temp/x-mort3-stset-allcauses", clear

*-------------------------------------------------------------------------------
* CALCULATE SMR AND AER BY DIAGNOSIS GROUP
*-------------------------------------------------------------------------------
// Calculate observed deaths, expected deaths, and person-years
gen _e  = (_t - _t0) * rate_allcauses
gen _y  = (_t - _t0) / 100000

// Expand dataset to create overall row
expand 2, gen(tag)
replace diag = 999 if tag == 1  // Use 999 for overall category
label define ldiag 999 "Overall", modify

// Collapse by diagnosis to get SMR and AER
preserve
collapse (sum) _d _e _y, by(diag)

// Calculate SMR and AER using custom smraer command
smraer _d _e _y, desaer(0) dessmr(1)

// Save SMR/AER results
tempfile smraer_results
save `smraer_results'
restore

*-------------------------------------------------------------------------------
* CALCULATE LIFE EXPECTANCY BY DIAGNOSIS (LIFE TABLE APPROACH)
*-------------------------------------------------------------------------------
preserve

// Sort by individual and time
sort indexno _t0

// Calculate conditional survival probabilities by diagnosis and age
// We'll use discrete age intervals for the life table

// Create age groups for life table (single year intervals)
gen age_int = floor((_t0 + _t) / 2)  // Mid-point age of interval

// Create person-time variable before collapse
gen person_time = _t - _t0

// Collapse to get deaths and person-time by diagnosis and age
collapse (sum) deaths = _d (sum) person_time, by(diag age_int)

// Calculate mortality rate (hazard) for each age-diagnosis combination
gen hazard = deaths / person_time

// Calculate survival probabilities
// For each diagnosis, we need cumulative survival from birth
bysort diag (age_int): gen cum_survival = exp(-sum(hazard))

// Replace with 1 for first age group (survival at birth = 1)
by diag: replace cum_survival = 1 if _n == 1
by diag: replace cum_survival = cum_survival[_n-1] * exp(-hazard) if _n > 1

// Calculate life expectancy as sum of survival probabilities
// LE = sum of person-years lived across all ages
collapse (sum) le_observed = cum_survival, by(diag)

// Save life expectancy results
tempfile le_results
save `le_results'
restore

*-------------------------------------------------------------------------------
* CALCULATE EXPECTED LIFE EXPECTANCY (GENERAL POPULATION RATES)
*-------------------------------------------------------------------------------
preserve

// Use general population rates to calculate expected survival
gen age_int = floor((_t0 + _t) / 2)

// Calculate expected deaths using general population rates
gen exp_deaths_age = (_t - _t0) * rate_allcauses
gen person_time = _t - _t0

// Collapse to get expected deaths and person-time by diagnosis and age
collapse (sum) exp_deaths = exp_deaths_age (sum) person_time, by(diag age_int)

// Calculate expected hazard based on general population
gen exp_hazard = exp_deaths / person_time

// Calculate expected cumulative survival
bysort diag (age_int): gen exp_cum_survival = exp(-sum(exp_hazard))
by diag: replace exp_cum_survival = 1 if _n == 1
by diag: replace exp_cum_survival = exp_cum_survival[_n-1] * exp(-exp_hazard) if _n > 1

// Calculate expected life expectancy
collapse (sum) le_expected = exp_cum_survival, by(diag)

// Save expected life expectancy results
tempfile le_expected_results
save `le_expected_results'
restore

*-------------------------------------------------------------------------------
* MERGE ALL RESULTS TOGETHER
*-------------------------------------------------------------------------------
use `smraer_results', clear

// Merge life expectancy results
merge 1:1 diag using `le_results', nogen
merge 1:1 diag using `le_expected_results', nogen

// Create formatted strings for life expectancies
gen str_le_obs = string(le_observed, "%9.1f")
gen str_le_exp = string(le_expected, "%9.1f")

// Calculate life years lost
gen le_lost = le_expected - le_observed
gen str_le_lost = string(le_lost, "%9.1f")

*-------------------------------------------------------------------------------
* PREPARE TABLE VARIABLES AND FORMATTING
*-------------------------------------------------------------------------------
// Sort by diagnosis (put Overall at top)
gsort -diag

// Create diagnosis label variable
decode diag, gen(strdiag)

// Create observed deaths with percentage
gen strobsperc = string(_d) + " (" + string((_d/_d[1])*100, "%9.1f") + "%)"

// Label variables for table headers
label variable strdiag "Cancer type"
label variable strobsperc "Observed deaths (%)"
label variable strexp "Expected deaths"
label variable smrstr "SMR (95% CI)"
label variable aerstr "AER* (95% CI)"
label variable str_le_obs "Observed LE"
label variable str_le_exp "Expected LE"
label variable str_le_lost "LE lost"

*-------------------------------------------------------------------------------
* EXPORT TO WORD DOCUMENT
*-------------------------------------------------------------------------------
cap: putdocx clear
putdocx begin, font(arial narrow, "10")
putdocx paragraph

// Create table with all relevant columns
putdocx table tbl1 = data("strdiag strobsperc strexp smrstr aerstr str_le_obs str_le_exp str_le_lost"), ///
varnames ///
border(start, nil) border(insideV, nil) border(end, nil) width(100%) ///
layout(autofitcontents) border(all, nil) ///
title("Table 3. All-cause mortality by childhood cancer type: SMR, AER, and life expectancy") ///
note("SMR = Standardized Mortality Ratio; AER = Absolute Excess Risk per 100,000 person-years; LE = Life Expectancy in years" , font(arial narrow, "8"))

// Format table headers and cells
putdocx table tbl1(1,.), bold
putdocx table tbl1(2,.), bold

// Alignment
putdocx table tbl1(.,1), halign(left) bold
putdocx table tbl1(.,2), halign(center)
putdocx table tbl1(.,3), halign(center)
putdocx table tbl1(.,4), halign(center)
putdocx table tbl1(.,5), halign(center)
putdocx table tbl1(.,6), halign(center)
putdocx table tbl1(.,7), halign(center)
putdocx table tbl1(.,8), halign(center)

// Save and open document
putdocx save "$temp/table3_smraers_bydiag.docx", replace
shell "$temp/table3_smraers_bydiag.docx"

exit
