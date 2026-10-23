/*==============================================================================
 WORD TABLE: MORTALITY AGE GAP (from CSVs made by mort3-11c-rr-mortality-age-gap.do)
 Output: $output/table_mortality_age_gap.docx    RCR (Claude-assisted) 23 Oct 2026
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"

import delimited "$output/table_mortality_age_gap.csv", clear varnames(1) stringcols(_all)
keep sex age surv_rate_per1000 pop_rate_per1000 eqstr gapstr
destring surv_rate_per1000 pop_rate_per1000, replace
format surv_rate_per1000 pop_rate_per1000 %5.1f
gen str_surv = string(surv_rate_per1000, "%5.1f")
gen str_pop  = string(pop_rate_per1000,  "%5.2f")
keep sex age str_surv str_pop eqstr gapstr
label var sex "Sex"
label var age "Attained age of survivor"
label var str_surv "Survivor all-cause death rate per 1,000 PY (smoothed)"
label var str_pop "General population rate per 1,000 at same age (E&W 2022)"
label var eqstr "Equivalent general-population age (95% CI)"
label var gapstr "Mortality age gap, years (95% CI)"

putdocx clear
putdocx begin, font("arial narrow", 10) landscape
putdocx paragraph
putdocx table t1 = data("sex age str_surv str_pop eqstr gapstr"), varnames ///
	border(start, nil) border(insideV, nil) border(end, nil) width(100%) ///
	layout(autofitcontents) border(all, nil) ///
	title("Table. Mortality age gap: the general-population age with the same all-cause death rate as a survivor of a given attained age") ///
	note("Survivor rates smoothed with a Poisson model with restricted cubic splines of attained age (ages 20-70), by sex. Population rates: England & Wales 2022, single ages by log-linear interpolation of 5-year bands. CI from 200 Poisson bootstrap resamples; population rates treated as fixed.", font("arial narrow", 8))
putdocx table t1(1,.), bold

import delimited "$output/table_mortality_age_gap_bydiag.csv", clear varnames(1) stringcols(_all)
keep group sex age deaths eqstr gapstr
label var group "Childhood cancer"
label var sex "Sex"
label var age "Attained age"
label var deaths "Deaths (all ages)"
label var eqstr "Equivalent population age (95% CI)"
label var gapstr "Mortality age gap, years (95% CI)"
putdocx paragraph
putdocx table t2 = data("group sex age deaths eqstr gapstr"), varnames ///
	border(start, nil) border(insideV, nil) border(end, nil) width(100%) ///
	layout(autofitcontents) border(all, nil) ///
	title("Table. Mortality age gap at attained ages 30 and 50 for the four largest cancer groups") 
putdocx table t2(1,.), bold
putdocx save "$output/table_mortality_age_gap.docx", replace
exit
