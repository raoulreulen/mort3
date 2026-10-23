/*==============================================================================
 mort3-10c-rr-period-analysis.do
 PURPOSE
   "Contemporary" estimates: SMRs and AERs restricted to person-time (and
   deaths) in calendar years 2010-2022 (yeargrp>=2010), side by side with the
   full follow-up estimates.
   Part A: all causes in the list, overall.
   Part B: all causes, recurrence, SPN and circulatory by attained-age group
           (<20, 20-29, 30-39, 40-49, 50-59, 60+).
 INPUTS   $temp/x-mort3-stset-<cod>.dta  (do-file 02a); smraer.ado
 OUTPUTS  $output/table_period_2010_2022.docx ; $output/table_period_2010_2022.csv
 NOTES    recur has no population comparator: Obs and person-years only
          (crude rate per 100,000 PY given in place of the AER column).
 Date 23 Oct 2026   Author RCR (Claude-assisted)
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
version 18
set more off

local cod allcauses allcodexcrecur recur spn infection blood endocrine mental nervous ///
	circulation cardiac cerebrovascular othercirculatory respiratory digestive ///
	muscoskeletal genitourinary perinatal external suicide other
local codage allcauses recur spn circulation

local lab_allcauses       "all causes"
local lab_allcodexcrecur  "all causes excl. recurrence"
local lab_recur           "recurrence/progression"
local lab_spn             "SPN"
local lab_infection       "infection"
local lab_blood           "blood"
local lab_endocrine       "endocrine"
local lab_mental          "mental"
local lab_nervous         "nervous system"
local lab_circulation     "circulatory"
local lab_cardiac         "cardiac"
local lab_cerebrovascular "cerebrovascular"
local lab_othercirculatory "other circulatory"
local lab_respiratory     "respiratory"
local lab_digestive       "digestive"
local lab_muscoskeletal   "musculoskeletal"
local lab_genitourinary   "genitourinary"
local lab_perinatal       "perinatal/congenital"
local lab_external        "external causes"
local lab_suicide         "suicide"
local lab_other           "other"

*-------------------------------------------------------------------------------
* 1. Collapse: period (0 = full, 1 = 2010-22) x attained-age group, per cause
*-------------------------------------------------------------------------------
tempfile stacked
local i = 0
local n = 0
foreach x of local cod {
	local ++n
	use _d _t _t0 newrate ageband yeargrp using "$temp/x-mort3-stset-`x'", clear
	gen double _e = (_t-_t0)*newrate
	gen double _y = (_t-_t0)/100000
	gen agegrp = 5
	replace agegrp = 20 if ageband>=20 & ageband<30
	replace agegrp = 30 if ageband>=30 & ageband<40
	replace agegrp = 40 if ageband>=40 & ageband<50
	replace agegrp = 50 if ageband>=50 & ageband<60
	replace agegrp = 60 if ageband>=60
	gen byte p2010 = (yeargrp>=2010)
	* duplicate so that period 0 = all follow-up, period 1 = 2010-22 only
	expand 2, gen(dup)
	gen period = cond(dup==0, 0, 1)
	drop if period==1 & p2010==0
	collapse (sum) _d _e _y, by(period agegrp)
	gen str30 cod = "`lab_`x''"
	gen str20 codvar = "`x'"
	gen _order = `n'
	* overall rows (agegrp = 0)
	preserve
		collapse (sum) _d _e _y, by(period cod codvar _order)
		gen agegrp = 0
		tempfile ov
		save `ov'
	restore
	append using `ov'
	if "`x'"!="allcauses" & "`x'"!="recur" & "`x'"!="spn" & "`x'"!="circulation" keep if agegrp==0
	if `i'>0 append using `stacked'
	save `stacked', replace
	local i = 1
}
use `stacked', clear

*-------------------------------------------------------------------------------
* 2. SMR / AER
*-------------------------------------------------------------------------------
smraer _d _e _y, dessmr(1) desaer(0)
* recur: no comparator -> show crude rate per 100,000 PY in AER column
gen double crude = _d/_y
replace smrstr = "" if codvar=="recur"
replace aerstr = string(crude,"%9.0f") + "†" if codvar=="recur"
replace obsexp = strobs + "/-" if codvar=="recur"

label define agegrp 0 "all ages" 5 "<20" 20 "20-29" 30 "30-39" 40 "40-49" 50 "50-59" 60 "60+"
label values agegrp agegrp
decode agegrp, gen(agelab)
drop smr smrll smrul aer aerll aerul strobs strexp crude
rename (_d _e _y) (obs exp py)
reshape wide obs exp py obsexp smrstr aerstr, i(cod codvar _order agegrp agelab) j(period)
sort _order agegrp
order cod agelab obsexp0 smrstr0 aerstr0 obsexp1 smrstr1 aerstr1
label var cod     "Cause of death"
label var agelab  "Attained age"
label var obsexp0 "Full follow-up O/E"
label var smrstr0 "Full follow-up SMR (95% CI)"
label var aerstr0 "Full follow-up AER (95% CI)*"
label var obsexp1 "2010-2022 O/E"
label var smrstr1 "2010-2022 SMR (95% CI)"
label var aerstr1 "2010-2022 AER (95% CI)*"
export delimited cod codvar agelab obs0 exp0 py0 obsexp0 smrstr0 aerstr0 obs1 exp1 py1 obsexp1 smrstr1 aerstr1 ///
	using "$output/table_period_2010_2022.csv", replace

*-------------------------------------------------------------------------------
* 3. Word tables
*-------------------------------------------------------------------------------
cap putdocx clear
putdocx begin, font("arial narrow", 10) landscape

* Part A: overall, all causes
preserve
	keep if agegrp==0
	local nA = _N
	putdocx paragraph
	putdocx text ("Table. "), bold
	putdocx text ("SMRs and AERs for specific causes of death: full follow-up (to 31 Dec 2022) compared with person-time in calendar years 2010-2022 only")
	putdocx table tblA = data("cod obsexp0 smrstr0 aerstr0 obsexp1 smrstr1 aerstr1"), varnames ///
		border(start, nil) border(insideV, nil) border(end, nil) border(insideH, nil) border(all, nil) ///
		width(100%) layout(autofitcontents)
	putdocx table tblA(1,.), bold border(bottom, single)
	putdocx table tblA(.,1), halign(left) bold
	forvalues c = 2/7 {
		putdocx table tblA(.,`c'), halign(center)
	}
	local last = `nA'+1
	putdocx table tblA(`last',.), border(bottom, single)
restore

* Part B: by attained age for 4 causes
preserve
	keep if inlist(codvar,"allcauses","recur","spn","circulation")
	sort _order agegrp
	gen str30 cod_disp = cod
	bys _order (agegrp): replace cod_disp = "" if _n>1
	local nB = _N
	putdocx paragraph
	putdocx text ("Table (continued). "), bold
	putdocx text ("By attained age: all causes, recurrence/progression, SPN and circulatory deaths")
	putdocx table tblB = data("cod_disp agelab obsexp0 smrstr0 aerstr0 obsexp1 smrstr1 aerstr1"), varnames ///
		border(start, nil) border(insideV, nil) border(end, nil) border(insideH, nil) border(all, nil) ///
		width(100%) layout(autofitcontents)
	putdocx table tblB(1,.), bold border(bottom, single)
	putdocx table tblB(.,1), halign(left) bold
	putdocx table tblB(.,2), halign(left)
	forvalues c = 3/8 {
		putdocx table tblB(.,`c'), halign(center)
	}
	local last = `nB'+1
	putdocx table tblB(`last',.), border(bottom, single)
restore
putdocx paragraph
putdocx text ("O/E = observed/expected deaths; SMR with exact Poisson 95% CI; *AER per 100,000 person-years. 2010-2022 columns use only person-time and deaths occurring in calendar years 2010-2022 (expected deaths from contemporaneous England & Wales rates). †Recurrence/progression has no population comparator: crude death rate per 100,000 person-years shown instead of the AER."), font("arial narrow", 8)
putdocx save "$output/table_period_2010_2022.docx", replace
di "saved $output/table_period_2010_2022.docx"
exit
