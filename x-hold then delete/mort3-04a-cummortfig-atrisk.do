/*==============================================================================
FIGURE 1 - CUMULATIVE MORTALITY WITH NUMBER AT RISK TABLE
==============================================================================*/
local cod spn recur circulation respiratory external

local j = 0
local xtime 70

*-------------------------------------------------------------------------------
* CALCULATE EXPECTED OVERALL (ALLCAUSES)
*-------------------------------------------------------------------------------
use  "$temp/x-mort3-stset-allcauses"  , clear

/*
stexpect conditional , ratevar(newrate) 		///
out($temp/x-mort3-expected-allcause , replace) 	///
method(2) at(5(1)`xtime') npoints(10)
*/

*-------------------------------------------------------------------------------
* CALCULATE NUMBER AT RISK BY 10-YEAR AGE BANDS
*-------------------------------------------------------------------------------
use  "$temp/x-mort3-stset-allcauses"  , clear

* Calculate number at risk at the start of each 10-year age band
foreach age in 10 20 30 40 50 60 70 {
	count if  _t > `age'
	local n`age' = r(N)
	di "Age `age': `n`age''"
}

*-------------------------------------------------------------------------------
* STCOMPET PER COD
*-------------------------------------------------------------------------------
local i = 0
foreach x in `cod' {
	use  "$temp/x-mort3-prepforstset" , clear

	recode `x' (0=2) if allcauses==1
	stset dox, fail(`x'==1) id(indexno) origin(dob) entry(doe) scale(365.25)
	assert _st==1
	stcompet cuminc = ci , compet1(2)
	gen ci_`x' = cuminc*100 if _d==1
	keep if _d==1
	keep ci_`x' _t
	sort _t

	tempfile `x'
	save ``x''
}

*-------------------------------------------------------------------------------
* APPEND ALL FILES
*-------------------------------------------------------------------------------
use "$temp/x-mort3-expected-allcause" , clear
append using  `spn' `recur' `circulation' `respiratory' `external'
gen expected = (1-conditional)*100

*-------------------------------------------------------------------------------
* CREATE MAIN CUMULATIVE INCIDENCE GRAPH (WITH X-AXIS)
*-------------------------------------------------------------------------------
#delimit ;

twoway
(line ci_recur ci_spn ci_circulation ci_respiratory ci_external _t if (_t>=5 & _t<= 70),
 sort connect(J J J J) clp(solid solid solid solid)
 clw(med med med med) clc(black blue red green))
(lowess expected t_exp if (t_exp>=5 & t_exp<= 70), sort bw(0.3) clp(shortdash) clc(cyan))
	,
	xtitle("Attained age, years", size(small))
	ytitle("Cumulative mortality, %", size(small))
	ylabel(, angle(0) labsize(small))
	xlabel(5(5)70, labsize(small) nogrid)

	legend(on order(1 "recurrence" 2 "spn" 3 "circulation" 4 "respiratory" 5 "external" 6 "expected mortality")
	ring(0) position(11) size(small) rowgap(0.1) cols(1))

	plotregion(margin(b=1))
	graphregion(margin(b=1))

	name(main, replace)
;
#delimit cr



*-------------------------------------------------------------------------------
* CREATE RISK TABLE GRAPH
*-------------------------------------------------------------------------------
#delimit ;

twoway
	(scatteri 1 5 )
	,
	yscale(range(0 1) off)
	ylabel(none)
	xscale(range(5 70))
	xlabel(none)
	xtitle("")

	text(1 5 "Number at risk", place(e) size(small))
	text(0.95 15 "`n10'", size(small) place(c))
	text(0.95 25 "`n20'", size(small) place(c))
	text(0.95 35 "`n30'", size(small) place(c))
	text(0.95 45 "`n40'", size(small) place(c))
	text(0.95 55 "`n50'", size(small) place(c))
	text(0.95 65 "`n60'", size(small) place(c))
	text(0.95 75 "`n70'", size(small) place(c))
	
	plotregion(margin(t=0 b=0))
	graphregion(margin(t=0 b=0))

	name(risk, replace)
;
#delimit cr

*-------------------------------------------------------------------------------
* COMBINE GRAPHS
*-------------------------------------------------------------------------------
graph combine main risk, ///
	cols(1) ///
	imargin(0 0 0 0) ///
	graphregion(margin(zero)) ///
	ysize(6) xsize(5)
