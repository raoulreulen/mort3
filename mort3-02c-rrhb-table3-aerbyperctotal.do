/*==============================================================================
TABLE 3 - AER% OF TOTAL FOR EACH COD BY FOLLOW-UP OR ATTAINED AGE
(LOOK AT SPN PAPER - TAKE CODE FROM THERE)

MODIFIED BY HAMA 26/08/26
- local cod updated to match circulatory subcategories and neoplasm/recur split - 
  see stset do-file for full rationale (mort3-02a-hbrr-stsetdata-circsubcats)
- added label formatting: nested cardiac/cerebrovascular/othercirculatory under 
  circulatory, suicide under external
- changed AER denominator from 100,000 to 10,000 person-years, with desaer(1) 
  to keep one decimal place of precision at the smaller scale
- reorganised the order results appear in the table (pinned All causes / 
  Recurrence or progression / All causes except recurrence to the top); 
  some remaining ordering and formatting still done manually in Word
 ***IMPORTANT***
 I added p-value for linear trend in AER across attained age (from Raoul's 
  updated version of this do-file) --> CANT GET ANY RESULTS FOR THE P-VALUE --> Looks like it needs a package called rs.ado --> ASK RAOUL 
Original table 3 do-file by Raoul untouched, this is a modified copy 

==============================================================================*/


version 16

local cod allcauses recur allcodexcrecur spn ///
circulation cardiac cerebrovascular othercirculatory ///
external suicide ///
respiratory nervous digestive infection perinatal endocrine genitourinary mental muscoskeletal blood other
	
 							
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

	*---------------------------------------------------------------------------
	* P-VALUE FOR LINEAR TREND IN THE AER ACROSS ATTAINED AGE
	*---------------------------------------------------------------------------
	sort agec
	gen double agemid = (agec + agec[_n+1]) / 2
	replace agemid = agec + 10 if missing(agemid)

	local paer = .
	capture {
		quietly glm _d agemid, family(poisson) link(rs _e) lnoffset(_y)
		estimates store _fa
		quietly glm _d, family(poisson) link(rs _e) lnoffset(_y)
		estimates store _na
		quietly lrtest _fa _na
		local paer = r(p)
	}
	estimates clear

	gen double ptraer = `paer'
	drop agemid

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

//tidy up labels for display, indent nested subcategories under their parent
replace cod = "All causes" 				if cod=="allcauses"
replace cod = "Recurrence or progression" if cod=="recur"
replace cod = "All causes except recurrence" if cod=="allcodexcrecur"
replace cod = " Subsequent primary neoplasm" 	if cod=="spn"
replace cod = " Circulatory (all)" 			if cod=="circulation"
replace cod = "   Cardiac" 				if cod=="cardiac"
replace cod = "   Cerebrovascular" 		if cod=="cerebrovascular"
replace cod = "   Other circulatory" 		if cod=="othercirculatory"
replace cod = " External causes" 			if cod=="external"
replace cod = "   Suicide" 				if cod=="suicide"
replace cod = " Respiratory" 				if cod=="respiratory"
replace cod = " Nervous" 					if cod=="nervous"
replace cod = " Digestive" 					if cod=="digestive"
replace cod = " Infection" 					if cod=="infection"
replace cod = " Perinatal" 					if cod=="perinatal"
replace cod = " Endocrine" 					if cod=="endocrine"
replace cod = " Genitourinary" 				if cod=="genitourinary"
replace cod = " Mental" 						if cod=="mental"
replace cod = " Musculoskeletal" 			if cod=="muscoskeletal"
replace cod = " Blood" 						if cod=="blood"
replace cod = " Other" 						if cod=="other"


smraer _d _e _y , desaer(1) dessmr(1)

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
keep agec cod _d  aerstr pcaer strpcaer ptraer
order agec cod _d  aerstr pcaer strpcaer
levelsof  agec , local(age)
reshape wide _d aerstr pcaer strpcaer  , i(cod ptraer) j(agec)

* Pin these three to the top, in this exact order; everything else sorts by pcaer60
gen _pin = .
replace _pin = 1 if cod=="All causes"
replace _pin = 2 if cod=="Recurrence or progression"
replace _pin = 3 if cod=="All causes except recurrence"

gsort _pin -pcaer60
drop _pin

order cod strpcaer*  pcaer* _d* aerstr*

*-------------------------------------------------------------------------------
*FORMAT THE AER TREND P-VALUE
*    "<0.001", then "<0.01", otherwise 2 decimal places. Left blank where the
*    model could not be fitted (no positive excess risk to model).
*-------------------------------------------------------------------------------
gen str10 ptraerstr = ""
replace ptraerstr = string(ptraer, "%4.2f") if ptraer >= 0.01 & !missing(ptraer)
replace ptraerstr = "<0.01"                 if ptraer <  0.01 & ptraer >= 0.001
replace ptraerstr = "<0.001"                if ptraer <  0.001

*-------------------------------------------------------------------------------
* TABLE WORD
*-------------------------------------------------------------------------------	
cap: putdocx clear
putdocx begin,  font(arial narrow, "8") landscape  
putdocx paragraph

putdocx table tbl1 = data("cod _d0 aerstr0 strpcaer0  _d30 aerstr30 strpcaer30 _d40 aerstr40 strpcaer40 _d50 aerstr50 strpcaer50 _d60 aerstr60 strpcaer60 ptraerstr") , varnames 				///
border(start, nil) border(insideV, nil) border(end, nil) width(100%) 		 ///
layout(autofitcontents) border(all, nil) 									 ///
title("Table3. Absolute Excess Risk for Specific Causes of Death by Attained Age as a Proportion of the Total Specific Absolute Excess Risk") ///
note("*per 10,000 person-years. P for trend from a likelihood-ratio test comparing additive excess-hazard (relative survival) Poisson models with and without attained age, fitted as a continuous term using the midpoint of each band." , font(arial narrow, "7"))

putdocx table tbl1(1,.), bold
putdocx table tbl1(.,1), halign(left) bold
putdocx table tbl1(.,2), halign(left)
putdocx table tbl1(.,3), halign(center)

putdocx save 	"$temp/table3_aerperc_cod3paper.docx" , replace
shell 			"$temp/table3_aerperc_cod3paper.docx"
exit
