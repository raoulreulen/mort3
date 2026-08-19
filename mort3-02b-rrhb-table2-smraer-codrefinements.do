/*========================================================================================
TABLE 2: SIRs by SPN - creates table with SIRs for all separate SPNs in one 
small table 

MODIFIED BY HAMA 
- local cod updated to match circulatory subcategories and neoplasm/recur split - 
  see stset do-file for full rationale (mort3-02a-rrhb-stsetdata-codrefinements)
- added label formatting: nested cardiac/cerebrovascular/othercirculatory under 
  circulatory, suicide under external, and all cause-specific categories under 
  "all causes except recurrence"
- changed AER denominator from 100,000 to 10,000 person-years
- NOTE (19/08/2026): total person-years and expected numbers higher than Raoul's table 2 
from 2025, observed numbers remain identical. likely because this version uses  BCCSS_CoVariates_Jan26 (Jan 2026). 
will be flagged and checked with Raoul 
*=======================================================================================*/

local cod allcauses recur allcodexcrecur spn ///
circulation cardiac cerebrovascular othercirculatory ///
external suicide ///
respiratory nervous digestive infection perinatal endocrine genitourinary mental muscoskeletal blood other


*=========================================================================================
* SIRs check CALCULATE SIRs MANUALLY FOR EACH SITE (RCR 11 Dec 2020)
*=========================================================================================
estimates clear

local i = 0 

foreach x of local cod { 
use "$temp/x-mort3-stset-`x'"  , clear

*stsplit age60 , at(60) after(dob)
*keep if age60>=60

	gen _e  = (_t- _t0)*rate_`x'
	gen _y  = (_t-_t0)/10000
	
	collapse (sum) _d _e _y 
	gen strcod = "`x'"

	if `i'>0 append using `table2'
	tempfile table2
	save `table2'
	local i =1
	}
	
//tidy up labels for display, indent nested subcategories under their parent
replace strcod = "All causes" 				if strcod=="allcauses"
replace strcod = "Recurrence or progression" if strcod=="recur"
replace strcod = "All causes except recurrence" if strcod=="allcodexcrecur"
replace strcod = " Subsequent primary neoplasm" 	if strcod=="spn"
replace strcod = " Circulatory (all)" 			if strcod=="circulation"
replace strcod = "   Cardiac" 				if strcod=="cardiac"
replace strcod = "   Cerebrovascular" 		if strcod=="cerebrovascular"
replace strcod = "   Other circulatory" 	if strcod=="othercirculatory"
replace strcod = " External causes" 		if strcod=="external"
replace strcod = "   Suicide" 				if strcod=="suicide"
replace strcod = " Respiratory" 			if strcod=="respiratory"
replace strcod = " Nervous" 			if strcod=="nervous"
replace strcod = " Digestive" 				if strcod=="digestive"
replace strcod = " Infection" 				if strcod=="infection"
replace strcod = " Perinatal" 				if strcod=="perinatal"
replace strcod = " Endocrine" 				if strcod=="endocrine"
replace strcod = " Genitourinary" 			if strcod=="genitourinary"
replace strcod = " Mental" 					if strcod=="mental"
replace strcod = " Musculoskeletal" 		if strcod=="muscoskeletal"
replace strcod = " Blood" 					if strcod=="blood"
replace strcod = " Other" 					if strcod=="other"

smraer _d _e _y , desaer(1) dessmr(1)

*-----------------------------------------------------------------------------------------
* ADDITIONAL VARIABLES FOR IN TABLE 
*-----------------------------------------------------------------------------------------
gen seq =_n
gsort -seq


gen strobsperc = string(_d) + " (" + string((_d/_d[1])*100 ,"%9.1f") + "%" + ")"
gen obsrate = string(_d/(_y ) ,"%9.1f") 


*-----------------------------------------------------------------------------------------
* STATA INTO WORD
*-----------------------------------------------------------------------------------------
cap: putdocx clear
putdocx begin,  font(arial narrow, "10")
putdocx paragraph
putdocx table tbl1 = data("strcod strobsperc strexp smrstr aerstr") , varnames ///
border(start, nil) border(insideV, nil) border(end, nil) width(100%) 				   ///
layout(autofitcontents) border(all, nil) 											   ///
title("Table 2. Observed and expected deaths, standardised mortality ratio and absolute excess risk for specific causes of death") ///
note("*per 10,000 person-years" , font(arial narrow, "8"))


//bold and centre columns
putdocx table tbl1(1,.), bold
putdocx table tbl1(2,.), bold
putdocx table tbl1(.,1), halign(left) bold
putdocx table tbl1(.,2), halign(left)
putdocx table tbl1(.,3), halign(center)
putdocx table tbl1(.,4), halign(center)

putdocx save 	"$temp/table2_smrsaers.docx" , replace
shell 			"$temp/table2_smrsaers.docx" 

exit










