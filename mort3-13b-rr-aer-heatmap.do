/*==============================================================================
MORT3-13b: DATA FOR LEXIS-STYLE HEAT MAP OF AER (and SMR) BY ATTAINED AGE x
           DECADE OF DIAGNOSIS

PURPOSE
	For all causes, recurrence/progression, SPN, circulatory, respiratory and
	external causes: observed deaths, expected deaths, person-years, AER
	(per 100,000 person-years) and SMR with 95% CIs (smraer.ado), by attained
	age group (5-19, 20-29, 30-39, 40-49, 50-59, 60+) x decade of diagnosis
	(<1970, 1970-79, 1980-89, 1990-99, 2000-06). For recurrence there is no
	population comparator (newrate==0): the crude rate per 100,000 PY is
	reported in the aer column and SMR is left missing. Cells with <1000
	person-years are blanked (aer/smr set to missing; obs/py retained).
	The heat map itself is drawn in R (mort3-13c-rr-aer-heatmap.R).

INPUTS
	$temp/x-mort3-stset-<cod>.dta  (cod = allcauses recur spn circulation
	                                respiratory external), from do-file 02a
OUTPUTS
	$output/aer_heatmap_data.csv  (cause agegrp decade obs exp py aer aerll
	                               aerul smr smrll smrul blank)

Date: 23 Oct 2026    Author: RCR (Claude-assisted)
==============================================================================*/

do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"

local cods allcauses recur spn circulation respiratory external
tempfile stacked
local i = 0
local n = 0
foreach x of local cods {
	local ++n
	use "$temp/x-mort3-stset-`x'", clear
	keep if _st==1
	* attained-age group from the stsplit ageband (0 1 5 10 ... 85)
	gen agegrp = .
	replace agegrp = 1 if ageband<20
	replace agegrp = 2 if inrange(ageband,20,29)
	replace agegrp = 3 if inrange(ageband,30,39)
	replace agegrp = 4 if inrange(ageband,40,49)
	replace agegrp = 5 if inrange(ageband,50,59)
	replace agegrp = 6 if ageband>=60
	assert !mi(agegrp)
	assert !mi(newrate)
	gen _e = newrate*(_t-_t0)
	gen _y = (_t-_t0)/100000          // PY in units of 100,000 -> AER per 100,000 PY
	collapse (sum) _d _e _y, by(agegrp decdxcat)
	gen str cause = "`x'"
	gen _order = `n'
	if `i'>0 append using `stacked'
	save `stacked', replace
	local i = 1
}
use `stacked', clear
smraer _d _e _y, dessmr(1) desaer(0)

* recurrence: no comparator -> crude rate per 100,000 PY, no SMR
replace smr   = . if cause=="recur"
replace smrll = . if cause=="recur"
replace smrul = . if cause=="recur"
replace smrstr = "" if cause=="recur"
* with _e==0 smraer's aerll/aerul equal the exact Poisson CI of the crude rate: keep

gen py  = _y*100000
gen obs = _d
gen exp = _e
replace exp = . if cause=="recur"
gen blank = py<1000
foreach v in aer aerll aerul smr smrll smrul {
	replace `v' = . if blank==1
}
label define lage 1 "5-19" 2 "20-29" 3 "30-39" 4 "40-49" 5 "50-59" 6 "60+"
label values agegrp lage
decode agegrp, gen(agegrp_str)
gen str10 decade = ""
replace decade = "<1970"   if decdxcat==1900
replace decade = "1970-79" if decdxcat==1970
replace decade = "1980-89" if decdxcat==1980
replace decade = "1990-99" if decdxcat==1990
replace decade = "2000-06" if decdxcat==2000
assert decade!=""
replace cause = "all causes"             if cause=="allcauses"
replace cause = "recurrence/progression" if cause=="recur"
replace cause = "SPN"                    if cause=="spn"
replace cause = "circulatory"            if cause=="circulation"
replace cause = "respiratory"            if cause=="respiratory"
replace cause = "external"               if cause=="external"
sort _order agegrp decdxcat
drop agegrp
rename agegrp_str agegrp
keep cause agegrp decade obs exp py aer aerll aerul smr smrll smrul blank _order
order cause agegrp decade obs exp py aer aerll aerul smr smrll smrul blank
format exp py %12.1f
format aer aerll aerul %9.1f
format smr smrll smrul %9.2f
list cause agegrp decade obs exp py aer smr blank if cause=="all causes", noobs clean
drop _order
export delimited using "$output/aer_heatmap_data.csv", replace
exit
