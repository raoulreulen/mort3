/*==============================================================================
SUICIDE (SEPARATE PAPER? FOR HAMA) BY DIAG
COUDL BE USED TO GENERATE BY DIAG ALL SMRs 
==============================================================================*/

local cod suicide 

*===============================================================================
* SIRs check CALCULATE SIRs MANUALLY FOR EACH SITE (RCR 11 Dec 2020)
*===============================================================================
estimates clear

local i = 0 

foreach x of local cod { 
use "$temp/x-mort3-stset-`x'"  , clear

	gen _e  = (_t- _t0)*rate_`x'
	gen _y  = (_t-_t0)/10000
	
	//expand to create total
	expand 2 , gen(tag)
	replace diag = 0 if tag==0 //overall 
	label define ldiag 0 "overall" , modify
	
	
	//by diag 
	collapse (sum) _d _e _y , by(diag)

		
	gen strcod = "`x'"

	if `i'>0 append using `table2'
	tempfile table2
	save `table2'
	local i =1
	}
	
	
smraer _d _e _y , desaer(1) dessmr(1)
sort diag 

*-------------------------------------------------------------------------------
* ADDITIONAL VARIABLES FOR IN TABLE 
*-------------------------------------------------------------------------------
gen seq =_n
gsort seq 

gen strobsperc = string(_d) + " (" + string((_d/_d[1])*100 ,"%9.1f") + "%" + ")"
gen obsrate = string(_d/(_y ) ,"%9.1f") 

*-------------------------------------------------------------------------------
* STATA INTO WORD
*-------------------------------------------------------------------------------
cap: putdocx clear
putdocx begin,  font(arial narrow, "10")
putdocx paragraph
putdocx table tbl1 = data("strcod diag strobsperc strexp smrstr aerstr") , varnames ///
border(start, nil) border(insideV, nil) border(end, nil) width(100%) 				   ///
layout(autofitcontents) border(all, nil) 											   ///
title("Table x. Standardised mortality ratios and absolute excess risks for suicide deaths") ///
note("*per 10,000 person-years" , font(arial narrow, "8"))

//bold and centre columns
putdocx table tbl1(1,.), bold
putdocx table tbl1(2,.), bold
putdocx table tbl1(.,1), halign(left) bold
putdocx table tbl1(.,2), halign(left)
putdocx table tbl1(.,3), halign(center)
putdocx table tbl1(.,4), halign(center)

putdocx save 	"$temp/tablex_smrsaerssuicide.docx" , replace
shell 			"$temp/tablex_smrsaerssuicide.docx" 

exit











































exit
