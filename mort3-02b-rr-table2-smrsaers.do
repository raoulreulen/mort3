/*========================================================================================
TABLE 2: SIRs by SPN - creates table with SIRs for all separate SPNs in one 
small table 
*=======================================================================================*/

local cod allcauses neoplasm spn infection blood endocrine mental nervous circulation ///
respiratory digestive muscoskeletal genitourinary perinatal other external suicide


*=========================================================================================
* SIRs check CALCULATE SIRs MANUALLY FOR EACH SITE (RCR 11 Dec 2020)
*=========================================================================================
estimates clear

local i = 0 

foreach x of local cod { 
use "$temp/x-mort3-stset-`x'"  , clear

*stsplit age50 , at(50) after(dob)
*keep if age50>=50

	gen _e  = (_t- _t0)*rate_`x'
	gen _y  = (_t-_t0)/10000
	
	collapse (sum) _d _e _y 
	gen strcod = "`x'"

	if `i'>0 append using `table2'
	tempfile table2
	save `table2'
	local i =1
	}
	
	
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
title("Table 2. Standardised mortality ratios and absolute excess risks for specific causes of death") ///
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











































exit
