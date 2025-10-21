/*==============================================================================
TABLE 3 - AER% OF TOTAL FOR EACH COD BY FOLLOW-UP OR ATTAINED AGE
(LOOK AT SPN PAPER - TAKE CODE FROM THERE)
==============================================================================*/


version 16

local cod allcauses neoplasm spn infection blood endocrine mental nervous circulation ///
respiratory digestive muscoskeletal genitourinary perinatal other external 
	
 							
local i = 0
local j = 0 

*-------------------------------------------------------------------------------
* AERS by ATT AGE FOR EACH SPT
*-------------------------------------------------------------------------------
foreach x in `cod' {
	use "$temp/x-mort3-stset-`x'"  , clear
	
	/*
	//period analysis
	stsplit cutperiod , after(time=d(1/1/1900)) at(110) 
	replace cutperiod  = cutperiod + 1900
	keep if cutperiod==2010
	*/
	
	//age cut offs
	stsplit agec, at(0 30(10)60 110) after(time=dob)
	
	gen _e 	= rate_`x'*(_t-_t0)
	gen _y 	= (_t-_t0) /10000
		
	collapse (sum) _d _y _e, by(agec)
	gen str cod = "`x'"
	gen j = `++j'
	
	
if `i'>0 append using "$temp/x-mort3-aerprop" , force
save "$temp/x-mort3-aerprop" , replace
local i=1
}



*-------------------------------------------------------------------------------
* READ CREATED FILE
*------------------------------------------------------------------------------
use "$temp/x-mort3-aerprop" , clear

smraer _d _e _y , desaer(0) dessmr(1)

*-------------------------------------------------------------------------------
* RECALCULATE AERs
*------------------------------------------------------------------------------
*gen aer = (_d-_e)/_y



*-------------------------------------------------------------------------------
* CALCULATE PERCENTAGES
*-------------------------------------------------------------------------------	
bysort agec (j): gen pcaer = (aer/aer[1])*100
replace pcaer= 0 if pcaer<0
gen strpcaer = string(round(pcaer,1)) + "%"

/*----------------------------------------------------------------------------------------
* OUTPUT FOR EXPORT TO R (RUN R-SCRIPT: SPN-RR-01A-stackedbarchart.R)
* NEEDS TO BE RUN N WINDOWS MACHINE
* TIFF FORMAT FOR WORD DOCUMENT /  EPS FOR JOURNAL
*---------------------------------------------------------------------------------------*/
*drop if agec==0
replace pcaer = 0 if pcaer<0
replace aer = 0 if aer<0

bysort agec (aer): gen cod_order = (_N-_n) +1
bysort cod (agec): replace cod_order=cod_order[_N]

/*
savesome agec aer aerll aerul cod strpcaer j cod_order using ///
"$rdata/x-barchart-R" if cod!="allcauses" & agec>0 , replace
*/

*-------------------------------------------------------------------------------
* RESHAPE for TABLE
*-------------------------------------------------------------------------------	
keep agec cod _d  aerstr pcaer strpcaer 
order agec cod _d  aerstr pcaer strpcaer 
levelsof  agec , local(age)
reshape wide _d aerstr pcaer strpcaer  , i(cod) j(agec) 

gsort -pcaer60

order cod strpcaer*  pcaer* _d* aerstr*


*-------------------------------------------------------------------------------
* TABLE WORD
*-------------------------------------------------------------------------------	
cap: putdocx clear
putdocx begin,  font(arial narrow, "9") landscape  
putdocx paragraph

putdocx table tbl1 = data("cod _d0 aerstr0 strpcaer0  _d30 aerstr30 strpcaer30 _d40 aerstr40 strpcaer40 _d50 aerstr50 strpcaer50 _d60 aerstr60 strpcaer60") , varnames 				///
border(start, nil) border(insideV, nil) border(end, nil) width(100%) 		 ///
layout(autofitcontents) border(all, nil) 									 ///
title("Table3. Absolute Excess Risk for Specific Causes of Death by Attained Age as a Proportion of the Total Specific Absolute Excess Risk") 

//bold and centre columns
putdocx table tbl1(1,.), bold
putdocx table tbl1(2,.), bold

putdocx table tbl1(.,1), halign(left) bold
putdocx table tbl1(.,2), halign(left)
putdocx table tbl1(.,3), halign(center)

putdocx save 	"$temp/table3_aerperc_cod3paper.docx" , replace
shell 			"$temp/table3_aerperc_cod3paper.docx"
	
exit	
	




































exit
