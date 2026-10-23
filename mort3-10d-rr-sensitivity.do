/*==============================================================================
 mort3-10d-rr-sensitivity.do
 PURPOSE
   Sensitivity analyses for the supplement: SMR and AER (per 100,000 PY) for
   all causes and the main specific causes (recurrence, SPN, circulatory,
   respiratory, external) under:
     (a) main analysis (follow-up and rates to 31 Dec 2022)
     (b) follow-up truncated at 31 Dec 2019 (person-time with yeargrp>=2020
         dropped) to exclude the COVID years
     (c) excluding survivors who embarked (vitstat==1; all their records)
     (d) all causes: deaths of unknown cause (codcat==0, n=28) removed from
         the observed count (person-time retained)
     (e) deaths certified under ICD-10 vs ICD-6/7/8/9: ICD-10 was used for
         all deaths from 1 Jan 2001 (checked: icdver==10 <=> year of death
         >=2001), so person-time is split at calendar year 2001.  Crude
         check of coding comparability; confounded with calendar period and
         attained age.
 INPUTS   $temp/x-mort3-stset-<cod>.dta  (do-file 02a); smraer.ado
 OUTPUTS  $output/table_sensitivity.docx ; $output/table_sensitivity.csv
 NOTES    No flag2020 variable exists in x-mort3-prepforstset; embarkation is
          vitstat==1 (758 survivors; vitstat 2 = dead, 3 = alive).
          recur has no population comparator: crude rate shown as AER.
 Date 23 Oct 2026   Author RCR (Claude-assisted)
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
version 18
set more off

local cod allcauses recur spn circulation respiratory external
local lab_allcauses   "all causes"
local lab_recur       "recurrence/progression"
local lab_spn         "SPN"
local lab_circulation "circulatory"
local lab_respiratory "respiratory"
local lab_external    "external causes"

* check the embarkation / ICD assumptions once
use indexno vitstat icdver dod allcauses codcat using "$temp/x-mort3-prepforstset", clear
tab vitstat, m
assert allcauses==1 if vitstat==2 & dod!=.
count if vitstat==1
di "embarked survivors: " r(N)
gen yod = year(dod)
tab icdver if allcauses==1 & yod>=2001, m
tab icdver if allcauses==1 & yod<2001, m
assert icdver==10 if allcauses==1 & yod>=2001
assert icdver<10  if allcauses==1 & yod<2001
count if allcauses==1 & codcat==0
di "deaths of unknown cause: " r(N)

*-------------------------------------------------------------------------------
* 1. Build scenario-specific collapsed rows, stack
*-------------------------------------------------------------------------------
tempfile stacked
local i = 0
local n = 0
foreach x of local cod {
	local ++n
	use _d _t _t0 newrate yeargrp vitstat codcat icdver using "$temp/x-mort3-stset-`x'", clear
	gen double _e = (_t-_t0)*newrate
	gen double _y = (_t-_t0)/100000
	* scenario indicators
	gen byte s1 = 1                                  // (a) main
	gen byte s2 = (yeargrp<2020)                     // (b) to 31 Dec 2019
	gen byte s3 = (vitstat!=1)                       // (c) excl. embarked
	gen byte s4 = 1                                  // (d) excl. unknown COD (allcauses only)
	gen byte s5 = (yeargrp<2001)                     // (e) ICD-6/7/8/9 period
	gen byte s6 = (yeargrp>=2001)                    // (e) ICD-10 period
	gen double _d4 = _d
	if "`x'"=="allcauses" replace _d4 = 0 if _d==1 & codcat==0
	tempfile base
	save `base'
	forvalues s = 1/6 {
		use `base', clear
		if `s'==4 & "`x'"!="allcauses" continue
		keep if s`s'==1
		if `s'==4 replace _d = _d4
		collapse (sum) _d _e _y
		gen scen = `s'
		gen str30 cod = "`lab_`x''"
		gen str20 codvar = "`x'"
		gen _order = `n'
		if `i'>0 append using `stacked'
		save `stacked', replace
		local i = 1
	}
}
use `stacked', clear
sort _order scen

*-------------------------------------------------------------------------------
* 2. SMR / AER
*-------------------------------------------------------------------------------
smraer _d _e _y, dessmr(1) desaer(0)
gen double crude = _d/_y
replace smrstr = "" if codvar=="recur"
replace aerstr = string(crude,"%9.0f") + "†" if codvar=="recur"
replace obsexp = strobs + "/-" if codvar=="recur"
gen str8 pystr = string(_y*100000/1000,"%9.1f")   // person-years in thousands

label define scen 1 "(a) main: follow-up to 2022" 2 "(b) follow-up to 31 Dec 2019" ///
	3 "(c) excluding embarked survivors" 4 "(d) excluding deaths of unknown cause" ///
	5 "(e) deaths before 2001 (ICD-6/7/8/9)" 6 "(e) deaths 2001-2022 (ICD-10)"
label values scen scen
decode scen, gen(scenlab)
gen str30 cod_disp = cod
bys _order (scen): replace cod_disp = "" if _n>1
rename (_d _e _y) (obs exp py)
replace py = py*100000
order cod scenlab obsexp pystr smrstr aerstr
label var cod     "Cause of death"
label var scenlab "Scenario"
label var obsexp  "O/E"
label var pystr   "PY (1000s)"
label var smrstr  "SMR (95% CI)"
label var aerstr  "AER (95% CI)*"
export delimited cod codvar scen scenlab obs exp py obsexp pystr smrstr aerstr using "$output/table_sensitivity.csv", replace
local nrows = _N

*-------------------------------------------------------------------------------
* 3. Word table
*-------------------------------------------------------------------------------
cap putdocx clear
putdocx begin, font("arial narrow", 10)
putdocx paragraph
putdocx text ("Table. "), bold
putdocx text ("Sensitivity analyses: SMRs and AERs for all causes and the main specific causes of death under alternative analysis assumptions")
putdocx table tbl1 = data("cod_disp scenlab obsexp pystr smrstr aerstr"), varnames ///
	border(start, nil) border(insideV, nil) border(end, nil) border(insideH, nil) border(all, nil) ///
	width(100%) layout(autofitcontents)
putdocx table tbl1(1,.), bold border(bottom, single)
putdocx table tbl1(.,1), halign(left) bold
putdocx table tbl1(.,2), halign(left)
forvalues c = 3/6 {
	putdocx table tbl1(.,`c'), halign(center)
}
local last = `nrows'+1
putdocx table tbl1(`last',.), border(bottom, single)
putdocx paragraph
putdocx text ("O/E = observed/expected deaths; PY = person-years; SMR with exact Poisson 95% CI; *AER per 100,000 person-years. (a) Main analysis: follow-up to 31 Dec 2022 with England & Wales rates to 2022. (b) Person-time and deaths in 2020-2022 excluded. (c) 758 survivors known to have embarked excluded entirely. (d) 28 deaths of unknown cause removed from the observed all-cause count (person-time retained). (e) Person-time and deaths split at 1 Jan 2001, when ICD-10 replaced ICD-9 for death certification: the two rows give SMRs for deaths coded under ICD-6/7/8/9 and ICD-10 respectively; this is confounded with calendar period and attained age and is only a crude check of coding comparability. †Recurrence/progression has no population comparator: crude death rate per 100,000 person-years shown instead of the AER."), font("arial narrow", 8)
putdocx save "$output/table_sensitivity.docx", replace
di "saved $output/table_sensitivity.docx"
exit
