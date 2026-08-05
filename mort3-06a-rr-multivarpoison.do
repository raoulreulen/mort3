*=========================================================================================
* MULTIVARIABLE ANALYSIS WITH CHECK FOR OVERDISPERSION - MULTI-SHEET VERSION
*=========================================================================================
estimates clear

* Define sites to analyze
local cod allcauses spn circulation external respiratory nervous digestive 

/*
infection blood endocrine mental nervous circulation ///
respiratory digestive muscoskeletal genitourinary perinatal othersuicide
*/


* Define output Excel file
local excelfile "$temp/x-mvpcsf-allsites.xlsx"

* Delete existing file if it exists
cap: erase "`excelfile'"

*=========================================================================================
* LOOP THROUGH EACH SITE
*=========================================================================================
foreach site of local cod {

	di as text _newline(2) "===== PROCESSING SITE: `site' =====" _newline

	use "$temp/x-mort3-stset-`site'"  , clear
	cap: drop y
	cap: drop pyrs

	*-------------------------------------------------------------------------------
	* ADD ATTAINED AGE WHEN NECESSARY
	*-------------------------------------------------------------------------------
	stsplit attage , at(0 20(10)50 110) after(time=dob)

	/*
	stsplit fu , at(0 10(10)40 100) after(time=survdate) trim
	label define fu 0 "0-9 yrs" 10 "10-19 yrs" 20 "20-29 yrs" 30 "30-39 yrs" 40 "40+ yrs"
	label values fu fu
	*/

	*-------------------------------------------------------------------------------
	* CALCULATE EXPECTED & PERSON-YEARS
	*-------------------------------------------------------------------------------
	gen pyrs = (_t-_t0)
	gen y    = pyrs/100000
	gen E 	 = (rate_`x')*pyrs

	*-------------------------------------------------------------------------------
	* COLLAPSE
	*-------------------------------------------------------------------------------
	local collap sex diag attage decdxcat agedxcat  /*rt*/ /* ct */
	collapse (sum) _d E y ,  by(`collap')

	*===============================================================================
	* CHECK FOR OVERDISPERSION,see page 178 Hilbe book; boundary LR-test
	*===============================================================================

	*-------------------------------------------------------------------------------
	* RR/SIR MODEL
	*-------------------------------------------------------------------------------
	glm _d i.(`collap')  if E!=0  , fam(pois) eform link(log)  lnoffset(E) //Poisson
	scalar lrpois =  e(ll)

	glm _d i.(`collap')  if E!=0  , fam(nb ml) eform link(log) lnoffset(E) //NBREG
	scalar lrnbreg =  e(ll)

	*
	nbreg _d i.(`collap')  if E!=0  , exposure(E)   irr

	scalar lr = -2*(lrpois-lrnbreg)
	di chiprob(1,lr)/2
	cap: assert (chiprob(1,lr)/2) < 0.05 //if true, then NB

	if !_rc {
		di "NEGATIVE BINOMIAL"
		local options1  fam(nb ml) eform lnoffset(E)  //NB2
		local model1 negbin
	}

	else {
		di "POISSON"
		local options1 fam(pois) eform link(log)  lnoffset(E)
		local model1 poisson
		}

	*-------------------------------------------------------------------------------
	* AER/RER MODEL
	*-------------------------------------------------------------------------------
	glm _d i.(`collap')  if E!=0  , fam(pois) eform link(rs E)  lnoffset(y) //Poisson
	scalar lrpois =  e(ll)
	
	glm _d i.(`collap')  if E!=0  , fam(nb ml) eform link(rs E)  lnoffset(y) //NBREG
	scalar lrnbreg =  e(ll)

	scalar lr = -2*(lrpois-lrnbreg)
	di chiprob(1,lr)/2
	cap: assert (chiprob(1,lr)/2) < 0.05 //if true, then NB


	//RUN NEG BINOMIAL
	if !_rc {
		di "NEGATIVE BINOMIAL"
		local options2 fam(nb ml) link(rs E) lnoffset(y)  eform
		local model2 negbin
	}

	//RUN POISSON
	else {
		di "POISSON"
		local options2 fam(pois) eform link(rs E) lnoffset(y)
		local model2 poisson
	}

	*-------------------------------------------------------------------------------
	* EXPORT TO SITE-SPECIFIC FILE
	*-------------------------------------------------------------------------------
	local append replace

	*-------------------------------------------------------------------------------
	* LOOP OVER COLLAPSED VARIABLES
	*-------------------------------------------------------------------------------
	foreach x of local collap {

		*---------------------------------------------------------------------------
		* CREATE DUMMY VARIABLES FOR VARIABLE OF INTEREST
		*---------------------------------------------------------------------------
		local depvar `x'
		local omit`depvar' : subinstr local collap "`depvar'" ""  //remove local
		qui: tabulate `depvar' , gen(`depvar'_)
		local max = r(r)

		*---------------------------------------------------------------------------
		* MULTIVARIABLE MODEL (1=RR, 2=RER)
		*---------------------------------------------------------------------------
		forvalues i=1/2 {
		_eststo  `x'`i': glm _d i.(`omit`depvar'') `depvar'_2-`depvar'_`max' if E!=0 , `options`i''

		//p-hetero
		_eststo m1: glm _d i.(`omit`depvar'')  if E!=0 & !mi(`depvar') , `options`i''
		lrtest m1 `x'`i'
		estadd scalar p_het = r(p) :`x'`i' //save p-value in scalar

		//p-trend
		glm _d  i.(`omit`depvar'') `depvar'  if E!=0 & !mi(`depvar') , `options`i''
		lrtest m1
		estadd scalar p_trend = r(p) :`x'`i' //save p-value in scalar
		estimates drop m1
		}

		*---------------------------------------------------------------------------
		* ESTOUT TO SITE-SPECIFIC FILE
		*---------------------------------------------------------------------------
		estout `x'* 															///
		using "$temp/temp_`site'_results.txt" , 								///
			cells("b(fmt(%9.1f)) & ci(par( ( - ) ) fmt(%9.1f)) ") eform 		///
			label collabels(, none) eqlabels(none) mlabel(, none) 			 	///
			stats(p_trend p_het , fmt(%9.2f)) 								 	///
			refcat(`x'_2  ,  label(1.0 (ref.))) substitute("0.0 (0.0-.)" "-")	///
			keep(`depvar'*) `append'

		local append append
		estimates clear

	} // End collapsed variables loop

	*-------------------------------------------------------------------------------
	* IMPORT TEXT FILE AND EXPORT TO SHEET IN MAIN WORKBOOK
	*-------------------------------------------------------------------------------
	import delimited using "$temp/temp_`site'_results.txt", clear delimiter(tab) varnames(nonames)

	* Replace p-values of 0.00 with <0.001
	foreach v of varlist _all {
		cap: replace `v' = "<0.001" if strtrim(`v') == "0.00"
	}

	* Save copy for combined sheet & store model info
	save "$temp/temp_combined_`site'", replace
	local model1_`site' `model1'
	local model2_`site' `model2'

	* Add footnote with model selection
	local N = _N + 2
	set obs `N'
	replace v1 = "RR/SIR model: `model1'; AER/RER model: `model2'" in `N'

	local lastrow = _N
	export excel using "`excelfile'", sheet("`site'") sheetreplace

	* Format font: Aptos Narrow size 8
	putexcel set "`excelfile'", sheet("`site'") modify
	putexcel A1:Z`lastrow', overwritefmt font("Aptos Narrow", 8)

	* Clean up temp file
	cap: erase "$temp/temp_`site'_results.txt"

} // End site loop

*=========================================================================================
* COMBINED SHEET - ALL SITES SIDE BY SIDE
*=========================================================================================
di as text _newline(2) "===== BUILDING COMBINED SHEET =====" _newline

* Merge all sites side by side
local snum 0
foreach site of local cod {
	local snum = `snum' + 1

	use "$temp/temp_combined_`site'", clear
	gen _rownum = _n
	rename v2 `site'_rr
	rename v3 `site'_aer

	if `snum' == 1 {
		save "$temp/temp_combined_all", replace
	}
	else {
		drop v1
		save "$temp/temp_combined_merge", replace
		use "$temp/temp_combined_all", clear
		merge 1:1 _rownum using "$temp/temp_combined_merge", nogen
		save "$temp/temp_combined_all", replace
		cap: erase "$temp/temp_combined_merge.dta"
	}
	cap: erase "$temp/temp_combined_`site'.dta"
}

use "$temp/temp_combined_all", clear
drop _rownum

* Insert header row at top (site names)
set obs `=_N + 1'
gen _sort = _n
replace _sort = 0 in `=_N'
sort _sort
drop _sort

replace v1 = "" in 1
foreach site of local cod {
	replace `site'_rr  = "`site'" in 1
	replace `site'_aer = "" in 1
}

* Add footnote rows with model info per site
local N = _N + 2
set obs `N'
replace v1 = "Model selection per site (poisson/negbin):" in `N'

foreach site of local cod {
	local N = _N + 1
	set obs `N'
	replace v1 = "  `site': RR/SIR=`model1_`site'', AER/RER=`model2_`site''" in `N'
}

* Export combined sheet
local lastrow = _N
export excel using "`excelfile'", sheet("Combined") sheetreplace

* Format font: Aptos Narrow size 8
putexcel set "`excelfile'", sheet("Combined") modify
putexcel A1:ZZ`lastrow', overwritefmt font("Aptos Narrow", 8)

cap: erase "$temp/temp_combined_all.dta"

di as result _newline(2) "===== ANALYSIS COMPLETE =====" _newline
di as text "Results saved to: `excelfile'"

exit
