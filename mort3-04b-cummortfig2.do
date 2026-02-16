/*==============================================================================
FIGURE 2 like in JAMA PAPER - CUM MORT ALL CAUSES 

==============================================================================*/
local cod spn recur circulation respiratory 


local j = 0
local xtime 70

*-------------------------------------------------------------------------------
* CALCULATE EXPECTED OVERALL (ALLCAUES) - VERY SLOW!
*-------------------------------------------------------------------------------
use  "$temp/x-mort3-stset-allcauses"  , clear


stexpect conditional , ratevar(newrate) 		///
out($temp/x-mort3-expected-allcause , replace) 	///
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
* FIGURE — build line variables, styles and legend from cod local
*-------------------------------------------------------------------------------
// Colours to cycle through (add more if needed)
local colours black blue red green orange purple

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
	local legorder `legorder' `k' "`x'"
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
	ylabel(0(1)16, angle(0)) xlabel(5(5)`xtime')

	/*LEGEND*/
	legend(on order(`legorder')
	ring(0) position(11) size(small) rowgap(0.1) cols(1))

;
#delimit cr
