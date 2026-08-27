/*==============================================================================
FIGURE 1 like in JAMA PAPER - CUM MORT ALL CAUSES 

MODIFIED BY HAMA 26/08/26
- local cod already uses recur/allcodexcrecur from earlier restructuring
- had to fix stexpect's output() option - it doesn't like full paths, needed 
  cd "$temp" first then just the bare filename
- had to install stexpect and stcompet (SSC packages, not on this machine before)
- added a small table with exact values at each 5-year age mark from the graph, 
  doesn't change the graph
- expected mortality now reflects the updated 2024 rates

Original figure 1 do-file by Raoul untouched, this is a modified copy 
==============================================================================*/
==============================================================================*/
local cod recur allcodexcrecur

/*
blood endocrine mental nervous circulation ///
respiratory digestive muscoskeletal genitourinary perinatal other external 
*/

local j = 0
local xtime 70

*-------------------------------------------------------------------------------
* CALCULATE EXPECTED OVERALL (ALLCAUES) - VERY SLOW!
*-------------------------------------------------------------------------------
use  "$temp/x-mort3-stset-allcauses"  , clear

cd "$temp"
stexpect conditional , ratevar(newrate) 		///
output(x-mort3-expected-allcause, replace) 	///
method(2) at(5(1)`xtime') npoints(10)


*-------------------------------------------------------------------------------
* STCOMPET PER COD
*-------------------------------------------------------------------------------
local i = 0
foreach x in `cod' {
	use  "$temp/x-mort3-prepforstset" , clear //read the original file
	
	//other cod to 2 as competing risk
	recode `x' (0=2) if allcauses==1
	
	//stset on each cause of death
	stset dox, fail(`x'==1) id(indexno) origin(dob) entry(doe) scale(365.25)
	assert _st==1
	stcompet cuminc = ci , compet1(2) //put other causes as competinng
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
append using  `recur' `allcodexcrecur'
gen expected = (1-conditional)*100

*-------------------------------------------------------------------------------
* ADDED BY HAMA: table of exact cumulative incidence/mortality values at each 
* 5-year age mark shown on the graph, for accurate reporting in text/tables 
* alongside the figure (doesn't affect the graph itself)
*-------------------------------------------------------------------------------
preserve
tempname results
postfile `results' age recur_ci allcodexcrecur_ci exp_mort using "$temp/figure1_keyages.dta", replace

foreach age in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 {
    quietly sum ci_recur if _t<=`age' & !missing(ci_recur), meanonly
    local rec = r(max)
    
    quietly sum ci_allcodexcrecur if _t<=`age' & !missing(ci_allcodexcrecur), meanonly
    local acer = r(max)
    
    quietly sum expected if t_exp<=`age' & !missing(expected), meanonly
    local exp = r(max)
    
    post `results' (`age') (`rec') (`acer') (`exp')
}
postclose `results'

use "$temp/figure1_keyages.dta", clear
format recur_ci allcodexcrecur_ci exp_mort %5.2f
list, clean noobs
restore

*-------------------------------------------------------------------------------
* FIGURE 
*-------------------------------------------------------------------------------
#delimit ;

twoway

/* CAUSES OF DEATH*/
(line ci_recur ci_allcodexcrecur _t if (_t>=5 & _t<= 70) 
, sort connect(J J J J) clp(solid solid solid solid) 
clw(med med med med) clc(black blue red green))
	
/* EXPECTED DEATHS*/
(lowess expected t_exp if (t_exp>=5 & t_exp<= 70), sort bw(0.3) clp(l) clc(cyan) clp(shortdash))	
	
	/* OPTIONS*/
	,
	
	/*LABELS*/
	xtitle("Attained age, years" ,size(small))
	ytitle("Cumulative mortality, %" , size(small))
	ylabel(0(5)30, angle(0)) xlabel(5(5)70)
	
	/*LEGEND*/
	legend(on) 
	legend(order(1 "recurrence" 2 "all causes except recurrence" 3 "expected mortality")
	ring(0) position(11) size(small) rowgap(0.1) cols(1))

;
#delimit cr
