/*==============================================================================
TABLE 3 - AER% OF TOTAL FOR EACH COD BY FOLLOW-UP OR ATTAINED AGE
(LOOK AT SPN PAPER - TAKE CODE FROM THERE)
==============================================================================*/


version 16

local cod allcauses neoplasm spn infection blood endocrine mental nervous circulation ///
respiratory digestive muscoskeletal genitourinary perinatal other external

*-------------------------------------------------------------------------------
* ATTAINED-AGE BANDS - change these two settings and the rest of the do-file
* (the table columns, the sort order and the trend test) follows automatically
*-------------------------------------------------------------------------------
// Cut-points passed to stsplit. Bands run from each cut-point to the next,
// so "0 30(10)60 110" gives 0-29, 30-39, 40-49, 50-59, 60+.
local agecuts 0 30(10)60 110

// Assumed half-width of the open-ended TOP band, used only to give it a
// midpoint for the trend test (60 + 10 = 70).
local tophalf 10


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
	stsplit agec, at(`agecuts') after(time=dob)
	
	gen _e 	= rate_`x'*(_t-_t0)
	gen _y 	= (_t-_t0) /100000
		
	collapse (sum) _d _y _e, by(agec)

	*---------------------------------------------------------------------------
	* P-VALUE FOR LINEAR TREND IN THE AER ACROSS ATTAINED AGE
	*   Additive excess-hazard model (rs link, person-years offset) with and
	*   without attained age, compared by likelihood-ratio test.
	*   Trend covariate = midpoint of each band, derived from the band edges
	*   themselves so it follows `agecuts' automatically. -agec- holds the
	*   lower edge, so the midpoint is halfway to the next band; the open-ended
	*   top band uses `tophalf'.
	*   Fitted on the collapsed rows: attained age is the only covariate, so
	*   the per-band totals are sufficient statistics and give the same LR test.
	*   Needs rs.ado (the relative-survival link for -glm-).
	*---------------------------------------------------------------------------
	sort agec
	gen double agemid = (agec + agec[_n+1]) / 2
	replace agemid = agec + `tophalf' if missing(agemid)

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

	// Constant within COD, so it rides through the reshape in i()
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
keep agec cod _d  aerstr pcaer strpcaer ptraer
order agec cod _d  aerstr pcaer strpcaer
levelsof  agec , local(age)
// ptraer is constant within COD, so it goes in i() and stays a single
// column instead of being reshaped into one per age band
reshape wide _d aerstr pcaer strpcaer  , i(cod ptraer) j(agec)

// Sort rows by the percentage in the OLDEST band (whichever that now is)
local nbands : word count `age'
local topband : word `nbands' of `age'
gsort -pcaer`topband'

order cod strpcaer*  pcaer* _d* aerstr*

*-------------------------------------------------------------------------------
* FORMAT THE AER TREND P-VALUE
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

// Build the column list from the age bands actually present, so the table
// follows `agecuts' automatically; the trend p-value is the final column
local vlist
foreach a of local age {
	local vlist `vlist' _d`a' aerstr`a' strpcaer`a'
}
local vlist `vlist' ptraerstr

local ncols : word count `vlist'
local ncols = `ncols' + 1          // +1 for the cod column

local tnote3 "P for trend from a likelihood-ratio test comparing additive excess-hazard (relative survival) Poisson models with and without attained age, fitted as a continuous term using the midpoint of each band. Blank where the model could not be fitted because the AER is not positive."

putdocx table tbl1 = data("cod `vlist'") , varnames 						 ///
border(start, nil) border(insideV, nil) border(end, nil) width(100%) 		 ///
layout(autofitcontents) border(all, nil) 									 ///
title("Table3. Absolute Excess Risk for Specific Causes of Death by Attained Age as a Proportion of the Total Specific Absolute Excess Risk") ///
note("`tnote3'" , font(arial narrow, "7"))

//bold and centre columns
putdocx table tbl1(1,.), bold
putdocx table tbl1(2,.), bold

putdocx table tbl1(.,1), halign(left) bold
putdocx table tbl1(.,2), halign(left)
putdocx table tbl1(.,3), halign(center)

// P for trend column (last one, whatever the number of age bands)
putdocx table tbl1(.,`ncols'), halign(center)

putdocx save 	"$temp/table3_aerperc_cod3paper.docx" , replace
shell 			"$temp/table3_aerperc_cod3paper.docx"
	
exit	
	




































exit
