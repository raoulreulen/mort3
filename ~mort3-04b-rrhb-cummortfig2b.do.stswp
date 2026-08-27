/*==============================================================================
FIGURE 2 like in JAMA PAPER - CUM MORT ALL CAUSES 

SPLIT BY HAMA INTO FIGURE 2A AND FIGURE 2B (26/08/26) - too many categories 
in one figure was too crowded to read.

FIGURE 2A: recurrence, SPN, circulatory (all), external, respiratory 
FIGURE 2B (THIS DO-FILE): circulatory subcategories (cardiac, cerebrovascular, other circulatory)

MODIFIED BY HAMA 27/08/26
- had to fix stexpect's output() option, needed cd "$temp" first then just 
  the filename, no full path
- modified labels within the legend
- changed y axis from 0-16 to 0-20 since recurrence line was close to the top
- added a table with exact values at 5 year age marks for reporting alongside the figure

Original figure 2 do-file by Raoul untouched, this is a modified copy 
==============================================================================*/
local cod cardiac cerebrovascular othercirculatory


local j = 0
local xtime 70

*-------------------------------------------------------------------------------
* CALCULATE EXPECTED OVERALL (ALLCAUES) - VERY SLOW!
*-------------------------------------------------------------------------------
use  "$temp/x-mort3-stset-allcauses"  , clear

cd "$temp"
stexpect conditional , ratevar(newrate) 		///
output(x-mort3-expected-allcause , replace) 	///
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
foreach x of local cod {
	append using ``x''
}
gen expected = (1-conditional)*100

*-------------------------------------------------------------------------------
* ADDED BY HAMA: table of exact cumulative incidence/mortality values at each 
* 5-year age mark shown on the graph, for accurate reporting in text/tables 
* alongside the figure (doesn't affect the graph itself)
*-------------------------------------------------------------------------------
preserve
tempname results
postfile `results' age cardiac_ci cerebrovascular_ci othercirculatory_ci exp_mort using "$temp/figure2b_keyages.dta", replace

foreach age in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 {
    quietly sum ci_cardiac if _t<=`age' & !missing(ci_cardiac), meanonly
    local v_cardiac = r(max)
    
    quietly sum ci_cerebrovascular if _t<=`age' & !missing(ci_cerebrovascular), meanonly
    local v_cerebrovascular = r(max)
    
    quietly sum ci_othercirculatory if _t<=`age' & !missing(ci_othercirculatory), meanonly
    local v_othercirculatory = r(max)
    
    quietly sum expected if t_exp<=`age' & !missing(expected), meanonly
    local v_exp = r(max)
    
    post `results' (`age') (`v_cardiac') (`v_cerebrovascular') (`v_othercirculatory') (`v_exp')
}
postclose `results'

use "$temp/figure2b_keyages.dta", clear
format cardiac_ci cerebrovascular_ci othercirculatory_ci exp_mort %5.2f
list, clean noobs
restore

*-------------------------------------------------------------------------------
* FIGURE — build line variables, styles and legend from cod local
*-------------------------------------------------------------------------------
// Colours to cycle through (add more if needed)
local colours black blue red

// Readable legend labels for each category
local lbl_cardiac "Cardiac"
local lbl_cerebrovascular "Cerebrovascular"
local lbl_othercirculatory "Other circulatory"

// Build: ci variable list, connect/clp/clw/clc options, legend order
local civars
local conn
local clps
local clws
local clcs
local legorder
local k = 0
foreach x of local cod {
	local ++k
	local civars `civars' ci_`x'
	local conn   `conn' J
	local clps   `clps' solid
	local clws   `clws' med
	local col : word `k' of `colours'
	local clcs   `clcs' `col'
	local legorder `legorder' `k' "`lbl_`x''"
}
// Expected mortality is the next series after the COD lines
local ++k
local legorder `legorder' `k' "expected mortality"

#delimit ;

twoway

/* CAUSES OF DEATH */
(line `civars' _t if (_t>=5 & _t<=`xtime')
, sort connect(`conn') clp(`clps') clw(`clws') clc(`clcs'))

/* EXPECTED DEATHS */
(lowess expected t_exp if (t_exp>=5 & t_exp<=`xtime')
, sort bw(0.3) clp(shortdash) clc(cyan))

	,

	/*LABELS*/
	xtitle("Attained age, years", size(small))
	ytitle("Cumulative mortality, %", size(small))
	ylabel(0(2)20, angle(0)) xlabel(5(5)`xtime')

	/*LEGEND*/
	legend(on order(`legorder')
	ring(0) position(11) size(small) rowgap(0.1) cols(1))

;
#delimit cr
