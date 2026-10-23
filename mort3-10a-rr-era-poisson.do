/*==============================================================================
 mort3-10a-rr-era-poisson.do
 PURPOSE
   Has excess mortality fallen in later treatment eras after adjusting for
   attained age?  Multiplicative Poisson models of the SMR (observed deaths,
   offset log(expected)) with decade of diagnosis as exposure, unadjusted and
   adjusted for attained-age group, sex, diagnosis and age-at-diagnosis group.
   Relative SMR per decade (ref <1970) with 95% CI and p for trend (decade as
   a 1..5 score, likelihood-ratio test, adjusted model).
   Table A: all causes in the list, rows = cause (observed deaths / unadjusted
            RR / adjusted RR), columns = decade + p-trend.
   Table B: adjusted RR per decade restricted to attained age 20-39 and 40+
            for all causes, recurrence, SPN and circulatory.
 INPUTS   $temp/x-mort3-stset-<cod>.dta   (do-file 02a, rates to 2024)
 OUTPUTS  $output/table_era_poisson.docx ; $output/table_era_poisson.csv ;
          $output/table_era_poisson_agerestricted.csv
 NOTES    recur has no population comparator (newrate==0) so the model uses
          person-years as the exposure: RRs for recurrence are ratios of the
          crude death rate, not of the SMR.  Cells with E==0 are dropped.
 Date 23 Oct 2026   Author RCR (Claude-assisted)
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
version 18
set more off

local cod allcauses allcodexcrecur recur spn circulation cardiac cerebrovascular ///
	respiratory infection nervous digestive external
local decs 1900 1970 1980 1990 2000
local declab1900 "<1970"
local declab1970 "1970-79"
local declab1980 "1980-89"
local declab1990 "1990-99"
local declab2000 "2000-06"

* cause labels as used in the manuscript
local lab_allcauses      "all causes"
local lab_allcodexcrecur "all causes excl. recurrence"
local lab_recur          "recurrence/progression"
local lab_spn            "SPN"
local lab_circulation    "circulatory"
local lab_cardiac        "cardiac"
local lab_cerebrovascular "cerebrovascular"
local lab_respiratory    "respiratory"
local lab_infection      "infection"
local lab_nervous        "nervous system"
local lab_digestive      "digestive"
local lab_external       "external causes"

*-------------------------------------------------------------------------------
* helper: after a poisson fit, return RR (95% CI) strings per decade in locals
*         rr1900..rr2000 (ref = "1.0 (ref)"; inestimable = "-")
*-------------------------------------------------------------------------------
cap program drop _getrr
program define _getrr, rclass
	syntax, decs(numlist)
	foreach d of numlist `decs' {
		if `d'==1900 {
			return local rr`d' "1.0 (ref)"
			continue
		}
		cap local b  = _b[`d'.decdxcat]
		cap local se = _se[`d'.decdxcat]
		if _rc | `se'==0 | `se'==. {
			return local rr`d' "-"
		}
		else {
			local rr = exp(`b')
			local ll = exp(`b'-1.96*`se')
			local ul = exp(`b'+1.96*`se')
			return local rr`d' = string(`rr',"%9.2f") + " (" + string(`ll',"%9.2f") + "," + string(`ul',"%9.2f") + ")"
		}
	}
end

*-------------------------------------------------------------------------------
* 1. Build collapsed cell datasets per cause and fit models
*-------------------------------------------------------------------------------
tempname ph
tempfile res
postfile `ph' str30 cod str10 model str12 d1900 str22 r1900 str12 d1970 str22 r1970 ///
	str12 d1980 str22 r1980 str12 d1990 str22 r1990 str12 d2000 str22 r2000 ///
	str10 ptrend using `res', replace

tempname ph2
tempfile res2
postfile `ph2' str30 cod str10 agerest str22 r1900 str22 r1970 str22 r1980 str22 r1990 str22 r2000 ///
	str10 ptrend using `res2', replace

foreach x of local cod {
	di _n(2) "=================== `x' ==================="
	use _d _t _t0 newrate ageband decdxcat sex diag agedxcat using "$temp/x-mort3-stset-`x'", clear
	gen double E  = (_t-_t0)*newrate
	gen double PY = (_t-_t0)
	* attained-age group from the 5-y ageband (record lies wholly within band)
	gen agegrp = 5
	replace agegrp = 20 if ageband>=20 & ageband<30
	replace agegrp = 30 if ageband>=30 & ageband<40
	replace agegrp = 40 if ageband>=40 & ageband<50
	replace agegrp = 50 if ageband>=50 & ageband<60
	replace agegrp = 60 if ageband>=60
	collapse (sum) _d E PY, by(decdxcat agegrp sex diag agedxcat)
	gen decscore = (decdxcat==1900) + 2*(decdxcat==1970) + 3*(decdxcat==1980) + 4*(decdxcat==1990) + 5*(decdxcat==2000)

	if "`x'"=="recur" {
		local expo PY
		di "recur: no population comparator -> exposure = person-years (crude rate ratio)"
	}
	else {
		local expo E
		count if E==0
		di "cells with E==0 dropped: " r(N)
		drop if E==0
	}
	drop if PY==0

	* observed deaths per decade
	foreach d of local decs {
		qui sum _d if decdxcat==`d'
		local obs`d' = string(r(sum),"%9.0f")
	}

	* unadjusted
	poisson _d i.decdxcat, exposure(`expo') irr iterate(100)
	_getrr, decs(`decs')
	foreach d of local decs {
		local u`d' "`r(rr`d')'"
	}
	* adjusted
	poisson _d i.decdxcat i.agegrp i.sex i.diag i.agedxcat, exposure(`expo') irr iterate(100)
	di "converged = " e(converged)
	_getrr, decs(`decs')
	foreach d of local decs {
		local a`d' "`r(rr`d')'"
	}
	* p for trend (adjusted): LR test, decade as continuous score vs no decade
	qui poisson _d i.agegrp i.sex i.diag i.agedxcat, exposure(`expo') iterate(100)
	est store m0
	qui poisson _d decscore i.agegrp i.sex i.diag i.agedxcat, exposure(`expo') iterate(100)
	est store m1
	lrtest m0 m1
	local ptr = r(p)
	local ptrs = cond(`ptr'<0.001, "<0.001", string(`ptr',"%9.3f"))
	di "adjusted trend (per decade) RR = " exp(_b[decscore]) "  p-trend = `ptrs'"

	post `ph' ("`lab_`x''") ("unadjusted") ("`obs1900'") ("`u1900'") ("`obs1970'") ("`u1970'") ///
		("`obs1980'") ("`u1980'") ("`obs1990'") ("`u1990'") ("`obs2000'") ("`u2000'") ("")
	post `ph' ("`lab_`x''") ("adjusted") ("`obs1900'") ("`a1900'") ("`obs1970'") ("`a1970'") ///
		("`obs1980'") ("`a1980'") ("`obs1990'") ("`a1990'") ("`obs2000'") ("`a2000'") ("`ptrs'")

	*---------------------------------------------------------------------------
	* Table B: adjusted model restricted to attained age 20-39 and 40+
	*---------------------------------------------------------------------------
	if inlist("`x'","allcauses","recur","spn","circulation") {
		tempfile cells
		save `cells'
		foreach rest in 2039 40 {
			use `cells', clear
			if `rest'==2039 {
				keep if agegrp==20 | agegrp==30
				local restlab "20-39"
			}
			else {
				keep if agegrp>=40
				local restlab "40+"
			}
			di _n "---- `x' attained age `restlab' ----"
			tab decdxcat, sum(_d)
			poisson _d i.decdxcat i.agegrp i.sex i.diag i.agedxcat, exposure(`expo') irr iterate(100)
			_getrr, decs(`decs')
			foreach d of local decs {
				local b`d' "`r(rr`d')'"
			}
			* mark decades with zero observed deaths as inestimable
			foreach d of local decs {
				qui count if decdxcat==`d'
				if r(N)==0 local b`d' "-"
			}
			qui poisson _d i.agegrp i.sex i.diag i.agedxcat, exposure(`expo') iterate(100)
			est store m0
			qui poisson _d decscore i.agegrp i.sex i.diag i.agedxcat, exposure(`expo') iterate(100)
			est store m1
			lrtest m0 m1
			local ptr = r(p)
			local ptrs = cond(`ptr'<0.001, "<0.001", string(`ptr',"%9.3f"))
			post `ph2' ("`lab_`x''") ("`restlab'") ("`b1900'") ("`b1970'") ("`b1980'") ("`b1990'") ("`b2000'") ("`ptrs'")
		}
	}
}
postclose `ph'
postclose `ph2'

*-------------------------------------------------------------------------------
* 2. Table A: long layout, 3 rows per cause (deaths, unadjusted RR, adjusted RR)
*-------------------------------------------------------------------------------
use `res', clear
gen seq = _n
* build a deaths row from the unadjusted row
preserve
	keep if model=="unadjusted"
	replace model = "deaths (n)"
	foreach d of local decs {
		replace r`d' = d`d'
	}
	replace seq = seq - 0.5
	tempfile deaths
	save `deaths'
restore
append using `deaths'
sort seq
keep seq cod model r1900 r1970 r1980 r1990 r2000 ptrend
order cod model r1900 r1970 r1980 r1990 r2000 ptrend
label var cod    "Cause of death"
label var model  "Estimate"
foreach d of local decs {
	label var r`d' "`declab`d''"
}
label var ptrend "p-trend*"
* blank repeated cause labels for display
gen str30 cod_disp = cod
bys cod (seq): replace cod_disp = "" if _n>1
sort seq
export delimited cod model r1900 r1970 r1980 r1990 r2000 ptrend using "$output/table_era_poisson.csv", replace
local nA = _N

cap putdocx clear
putdocx begin, font("arial narrow", 10) landscape
putdocx paragraph
putdocx text ("Table. "), bold
putdocx text ("Relative SMRs by decade of childhood cancer diagnosis (reference <1970), unadjusted and adjusted for attained age, sex, type of childhood cancer and age at diagnosis (multiplicative Poisson models)")
putdocx table tblA = data("cod_disp model r1900 r1970 r1980 r1990 r2000 ptrend"), varnames ///
	border(start, nil) border(insideV, nil) border(end, nil) border(insideH, nil) border(all, nil) ///
	width(100%) layout(autofitcontents)
putdocx table tblA(1,.), bold border(bottom, single)
putdocx table tblA(.,1), halign(left) bold
putdocx table tblA(.,2), halign(left)
forvalues c = 3/8 {
	putdocx table tblA(.,`c'), halign(center)
}
local last = `nA'+1
putdocx table tblA(`last',.), border(bottom, single)
putdocx paragraph
putdocx text ("RR (95% CI) = ratio of SMRs relative to survivors diagnosed before 1970, from Poisson regression of observed deaths with log(expected deaths) as offset; adjusted models include attained-age group (5-19, 20-29, 30-39, 40-49, 50-59, 60+), sex, type of childhood cancer and age at diagnosis (0-3, 4-7, 8-11, 12-14). Cells with zero expected deaths were excluded. For recurrence/progression there is no population comparator, so the RRs are ratios of the crude death rate (offset log(person-years)). *p for trend from a likelihood-ratio test of decade fitted as a linear score (1-5) in the adjusted model."), font("arial narrow", 8)

*-------------------------------------------------------------------------------
* 3. Table B: adjusted RR restricted to attained age 20-39 and 40+
*-------------------------------------------------------------------------------
use `res2', clear
label var cod    "Cause of death"
label var agerest "Attained age"
foreach d of local decs {
	label var r`d' "`declab`d''"
}
label var ptrend "p-trend*"
export delimited using "$output/table_era_poisson_agerestricted.csv", replace
local nB = _N
putdocx paragraph
putdocx text ("Table (continued). "), bold
putdocx text ("Adjusted relative SMRs by decade of diagnosis, restricted to person-time at attained ages 20-39 and 40 or over")
putdocx table tblB = data("cod agerest r1900 r1970 r1980 r1990 r2000 ptrend"), varnames ///
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
putdocx paragraph
putdocx text ("Models as above, fitted within the age restriction (attained-age groups within the restriction retained as covariates). '-' = no person-time or no deaths in that decade within the age restriction (survivors diagnosed in 2000-06 had not reached age 40 by end of follow-up). Recurrence/progression: ratios of crude death rates."), font("arial narrow", 8)
putdocx save "$output/table_era_poisson.docx", replace
di "saved $output/table_era_poisson.docx"
exit
