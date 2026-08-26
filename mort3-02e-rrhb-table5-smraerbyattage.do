/*==============================================================================
TABLE 5: SMRs and AERs by ATTAINED AGE

PURPOSE
	For each cause of death (COD), this do-file:
	  1. Loads the survival-time (stset) dataset created in do-file 02a
	  2. Assigns each record to an attained-age band from its -ageband-
	  3. Calculates expected deaths (_e) and person-years (_y)
	  4. Collapses data by attained-age category (attagecat)
	  5. Tests for a linear trend in the SMR and the AER across attained age
	  6. Stacks all CODs into one dataset
	  7. Computes SMRs and AERs  (using the custom -smraer- command)
	  8. Reshapes the data wide (one row per COD, columns per attained-age band)
	  9. Exports a formatted Word table via putdocx

	Attained-age bands: <20, 20-29, 30-39, 40-49, 50-59, 60+

REQUIRES
	- Globals: $temp  (path to temp/output folder)
	- Datasets: $temp/x-mort3-stset-<cod>.dta  (from do-file 02a)
	- Ado-files: smraer.ado, rs.ado (the relative-survival link used by -glm-
	  for the AER trend test)
	  
MODIFIED BY HAMA 26/08/26
- local cod updated to match circulatory subcategories and neoplasm/recur split
- attained-age bands changed from <20/20-29 to a single combined 0-29 band, 
  to match Table 3
- changed AER denominator from 100,000 to 10,000 person-years (aerdec raised to 1)
- added nested label formatting (circulatory subcategories, suicide under external)
- NOTE: rs.ado not currently installed - AER trend p-values (ptraerstr) blank 
  until resolved with Raoul; SMR trend p-values unaffected

Original table 5 do-file by Raoul untouched, this is a modified copy
==============================================================================*/


*-------------------------------------------------------------------------------
* 0. SETTINGS - change these to alter the scaling / rounding of the table
*-------------------------------------------------------------------------------
// Person-year denominator for the AER (e.g. 10000 or 100000)
local aerper 10000

// Decimal places. NB: per 10,000 person-years the AERs are 10x smaller than
// per 100,000, so at 0 dp a small AER can display as "0" (and a small negative
// lower limit as "-0"). Raise aerdec to 1 if that matters for the rarer CODs.
local aerdec 1
local smrdec 1

// Comma-formatted version of `aerper' for the table note, e.g. 10,000
local aerperc = trim(string(`aerper', "%15.0fc"))


*-------------------------------------------------------------------------------
* 1. DEFINE CAUSES OF DEATH
*    Each name matches a stset dataset and a rate variable (rate_<cod>).
*-------------------------------------------------------------------------------
local cod allcauses recur allcodexcrecur spn ///
	circulation cardiac cerebrovascular othercirculatory ///
	external suicide ///
	respiratory nervous digestive infection perinatal endocrine genitourinary mental muscoskeletal blood other


*-------------------------------------------------------------------------------
* 2. LOOP OVER EACH COD: ASSIGN ATTAINED-AGE BAND, THEN CALCULATE
*    EXPECTED DEATHS & PERSON-YEARS. Then stack all CODs into one file.
*-------------------------------------------------------------------------------
tempfile stacked
local i = 0
local n = 0

foreach x of local cod {
	local ++n

	// Load only the variables needed - these files carry all the covariates
	use _d _t _t0 ageband rate_`x' using "$temp/x-mort3-stset-`x'", clear

	// Attained-age band. 02a has already stsplit at 0 1 5(5)85 110, and
	// 20/30/40/50/60 are all cut-points in that list, so every record already
	// lies entirely within one band. -ageband- (the band's lower edge) is
	// therefore enough - a second -stsplit- here would cut nothing and is by
	// far the slowest step in this do-file, so it is not used.
	gen byte attagecat = 1 + (ageband>=30) + (ageband>=40) ///
	                       + (ageband>=50) + (ageband>=60)

	// Guard: fails loudly if 02a's agebands ever stop nesting inside these
	// bands (i.e. if a record were to straddle an attained-age boundary)
	assert _t <= 10*attagecat + 20 + 0.001 if attagecat < 5

	// _e = expected deaths (individual person-time * expected rate)
	// _y = person-years in units of `aerper' (for AER scaling)
	gen _e = rate_`x' * (_t - _t0)
	gen _y = (_t - _t0) / `aerper'

	// Sum observed deaths (_d), person-years (_y), expected deaths (_e)
	// within each attained-age category
	collapse (sum) _d _y _e, by(attagecat)

	*---------------------------------------------------------------------------
	* P-VALUES FOR LINEAR TREND ACROSS ATTAINED AGE
	*   Fitted on the collapsed rows: attained age is the only covariate, so the
	*   per-band totals are sufficient statistics and give the same likelihood
	*   ratio as the record-level data.
	*   Trend covariate = midpoint of each band (60+ taken as 70).
	*---------------------------------------------------------------------------
	recode attagecat (1=15) (2=35) (3=45) (4=55) (5=70), gen(attmid)

	// SMR trend: multiplicative Poisson model, expected deaths as offset
	local psmr = .
	capture {
		quietly glm _d attmid, family(poisson) lnoffset(_e)
		estimates store _fs
		quietly glm _d, family(poisson) lnoffset(_e)
		estimates store _ns
		quietly lrtest _fs _ns
		local psmr = r(p)
	}

	// AER trend: additive excess-hazard model (rs link), person-years as offset
	local paer = .
	capture {
		quietly glm _d attmid, family(poisson) link(rs _e) lnoffset(_y)
		estimates store _fa
		quietly glm _d, family(poisson) link(rs _e) lnoffset(_y)
		estimates store _na
		quietly lrtest _fa _na
		local paer = r(p)
	}
	estimates clear

	// Constant within COD, so these ride through the reshape in i()
	gen double ptrsmr = `psmr'
	gen double ptraer = `paer'
	drop attmid

	// Tag rows with the COD name and preserve original order
	gen str cod = "`x'"
	gen _order = `n'

	// Stack: append previous iterations then save
	if `i'>0 append using `stacked', force
	save `stacked', replace
	local i = 1
}


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

label define attagecat 1 "0-29" 2 "30-39" 3 "40-49" 4 "50-59" 5 "60+"
label values attagecat attagecat



*-------------------------------------------------------------------------------
* 3. COMPUTE SMRs AND AERs
*    smraer creates: smrstr (SMR + 95% CI), aerstr (AER + 95% CI), strobs
*    (observed deaths). It also makes obsexp ("obs/exp") - not used here, the
*    table reports observed only.
*-------------------------------------------------------------------------------
smraer _d _e _y, desaer(`aerdec') dessmr(`smrdec')


*-------------------------------------------------------------------------------
* 4. RESHAPE WIDE: one row per COD, columns per attained-age category
*    After reshape, variables are named strobs1 smrstr1 aerstr1 ... strobs6 etc.
*-------------------------------------------------------------------------------
keep  attagecat cod strobs aerstr smrstr ptrsmr ptraer _order
order attagecat cod strobs smrstr aerstr

// Grab the age-band levels AND their value labels before reshaping, so the
// labels can be used as column headings in the Word table below
levelsof attagecat, local(ages)
foreach a of local ages {
	local lbl`a' : label (attagecat) `a'
}
local nbands : word count `ages'

// ptrsmr/ptraer are constant within COD, so they go in i() and stay as
// single columns rather than being reshaped into one per age band
reshape wide strobs smrstr aerstr, i(cod _order ptrsmr ptraer) j(attagecat)

// Sort rows to match the order defined in the cod local
sort _order
drop _order

// Format the trend p-values: "<0.001", then "<0.01", otherwise 2 decimals.
// (A missing p-value - model failed to converge - is left blank.)
foreach p in smr aer {
	gen str10 ptr`p'str = ""
	replace ptr`p'str = string(ptr`p', "%4.2f") if ptr`p' >= 0.01 & !missing(ptr`p')
	replace ptr`p'str = "<0.01"                 if ptr`p' <  0.01 & ptr`p' >= 0.001
	replace ptr`p'str = "<0.001"                if ptr`p' <  0.001
}

// Build column list dynamically from the attained-age levels so the table
// adapts automatically if the age bands are changed in future
local vlist
foreach a of local ages {
	local vlist `vlist' strobs`a' smrstr`a' aerstr`a'
}

// ... then the two trend p-value columns at the far right
local vlist `vlist' ptrsmrstr ptraerstr

// Count rows and columns (needed for borders and alignment below)
local nrows = _N
local ncols : word count `vlist'
local ncols = `ncols' + 1          // +1 for the cod column


*-------------------------------------------------------------------------------
* 5. EXPORT TO WORD TABLE (putdocx)
*    The table gets a two-row header: row 1 = attained-age band (spanning its
*    three columns), row 2 = Obs, SMR, AER.
*-------------------------------------------------------------------------------
cap putdocx clear
putdocx begin, font(arial narrow, "7") landscape
putdocx paragraph

// --- Title: "Table 5" in bold, remainder in normal weight ----
putdocx text ("Table 5"), bold
putdocx text (". SMRs and AERs for Specific Causes of Death by Attained Age")

local tnote "Obs = observed deaths. SMR = standardised mortality ratio. AER = absolute excess risk per `aerperc' person-years. P for trend from a likelihood-ratio test comparing Poisson models with and without attained age (band midpoint) fitted as a continuous term; the SMR test uses expected deaths as offset, the AER test an additive excess-hazard model with person-years as offset."

// --- Create the table from data ---
// border(all, nil) removes all borders; we add specific ones below
putdocx table tbl1 = data("cod `vlist'"), varnames ///
	border(start, nil) border(insideV, nil) border(end, nil) ///
	border(insideH, nil) border(all, nil) ///
	width(100%) layout(autofitcontents) ///
	note("`tnote'", font(arial narrow, "6"))

// --- Column alignment: do this now, while the table is still rectangular ---
putdocx table tbl1(.,1), halign(left) bold
forvalues c = 2/`ncols' {
	putdocx table tbl1(.,`c'), halign(center)
}

// --- Overwrite the varnames row (row 1) with readable sub-headings ---
putdocx table tbl1(1,1) = ("Cause of death")
local c = 1
foreach a of local ages {
	local ++c
	putdocx table tbl1(1,`c') = ("Obs")
	local ++c
	putdocx table tbl1(1,`c') = ("SMR (95% CI)")
	local ++c
	putdocx table tbl1(1,`c') = ("AER (95% CI)")
}

// --- Sub-headings for the two trend p-value columns ---
local ++c
putdocx table tbl1(1,`c') = ("SMR")
local ++c
putdocx table tbl1(1,`c') = ("AER")

// --- Insert a spanner row above it holding the attained-age band labels ---
putdocx table tbl1(1,.), addrows(1, before)

// Fill the spanner right-to-left: merging cells removes the ones to the
// right, so working backwards keeps the column indices still to come valid.
// The "P for trend" spanner is rightmost, so it goes first.
local pcol = 2 + 3*`nbands'
putdocx table tbl1(1,`pcol') = ("P for trend"), halign(center)
putdocx table tbl1(1,`pcol'), colspan(2)
putdocx table tbl1(1,`pcol'), border(bottom, single)

forvalues k = `nbands'(-1)1 {
	local a : word `k' of `ages'
	local c = 2 + 3*(`k'-1)
	putdocx table tbl1(1,`c') = ("`lbl`a'' years"), halign(center)
	putdocx table tbl1(1,`c'), colspan(3)
	putdocx table tbl1(1,`c'), border(bottom, single)
}

// --- Both header rows bold ---
putdocx table tbl1(1,.), bold
putdocx table tbl1(2,.), bold

// --- Horizontal line below the sub-heading row ---
putdocx table tbl1(2,.), border(bottom, single)

// --- Horizontal line at the bottom of the table ---
// Data rows now run from 3 to nrows+2 (rows 1-2 are the header)
local lastrow = `nrows' + 2
putdocx table tbl1(`lastrow',.), border(bottom, single)

// --- Save and open the document ---
putdocx save "$temp/table5_smraerbyattage_cod3paper.docx", replace
shell "$temp/table5_smraerbyattage_cod3paper.docx"

exit
