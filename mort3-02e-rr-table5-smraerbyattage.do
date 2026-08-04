/*==============================================================================
TABLE 5: SMRs and AERs by ATTAINED AGE

PURPOSE
	For each cause of death (COD), this do-file:
	  1. Loads the survival-time (stset) dataset created in do-file 02a
	  2. Splits follow-up time at attained ages 20, 30, 40, 50 and 60 (stsplit)
	  3. Calculates expected deaths (_e) and person-years (_y)
	  4. Collapses data by attained-age category (attagecat)
	  5. Stacks all CODs into one dataset
	  6. Computes SMRs and AERs  (using the custom -smraer- command)
	  7. Reshapes the data wide (one row per COD, columns per attained-age band)
	  8. Exports a formatted Word table via putdocx

	Attained-age bands: <20, 20-29, 30-39, 40-49, 50-59, 60+

REQUIRES
	- Globals: $temp  (path to temp/output folder)
	- Datasets: $temp/x-mort3-stset-<cod>.dta  (from do-file 02a)
	- Ado-file: smraer.ado
==============================================================================*/


*-------------------------------------------------------------------------------
* 0. SETTINGS - change these to alter the scaling / rounding of the table
*-------------------------------------------------------------------------------
// Person-year denominator for the AER (e.g. 10000 or 100000)
local aerper 10000

// Decimal places. NB: per 10,000 person-years the AERs are 10x smaller than
// per 100,000, so 1 dp stops rare CODs rounding away to "0". Set to 0 if you
// switch `aerper' back to 100000.
local aerdec 1
local smrdec 1

// Comma-formatted version of `aerper' for the table note, e.g. 10,000
local aerperc = trim(string(`aerper', "%15.0fc"))


*-------------------------------------------------------------------------------
* 1. DEFINE CAUSES OF DEATH
*    Each name matches a stset dataset and a rate variable (rate_<cod>).
*-------------------------------------------------------------------------------
local cod allcauses neoplasm spn infection blood endocrine mental nervous ///
	circulation respiratory digestive muscoskeletal genitourinary perinatal ///
	other external suicide


*-------------------------------------------------------------------------------
* 2. LOOP OVER EACH COD: SPLIT ON ATTAINED AGE, THEN CALCULATE
*    EXPECTED DEATHS & PERSON-YEARS. Then stack all CODs into one file.
*-------------------------------------------------------------------------------
tempfile stacked
local i = 0
local n = 0

foreach x of local cod {
	local ++n

	// Load the stset dataset for this COD
	use "$temp/x-mort3-stset-`x'", clear

	// Analysis time is attained age (origin(dob), scale(365.25) in 02a), so
	// splitting at 20/30/40/50/60 allocates each person's follow-up to the
	// correct attained-age band. These cut-points nest inside the 5-year
	// agebands already created in 02a, so no record is cut in half twice.
	stsplit attage, at(0 20 30 40 50 60 110)

	recode attage (0 = 1) (20 = 2) (30 = 3) (40 = 4) (50 = 5) (60 110 = 6), ///
		gen(attagecat)

	// _e = expected deaths (individual person-time * expected rate)
	// _y = person-years in units of `aerper' (for AER scaling)
	gen _e = rate_`x' * (_t - _t0)
	gen _y = (_t - _t0) / `aerper'

	// Sum observed deaths (_d), person-years (_y), expected deaths (_e)
	// within each attained-age category
	collapse (sum) _d _y _e, by(attagecat)

	// Tag rows with the COD name and preserve original order
	gen str cod = "`x'"
	gen _order = `n'

	// Stack: append previous iterations then save
	if `i'>0 append using `stacked', force
	save `stacked', replace
	local i = 1
}

label define attagecat 1 "<20" 2 "20-29" 3 "30-39" 4 "40-49" 5 "50-59" 6 "60+"
label values attagecat attagecat


*-------------------------------------------------------------------------------
* 3. COMPUTE SMRs AND AERs
*    smraer creates: smrstr (SMR + 95% CI), aerstr (AER + 95% CI), obsexp
*-------------------------------------------------------------------------------
smraer _d _e _y, desaer(`aerdec') dessmr(`smrdec')


*-------------------------------------------------------------------------------
* 4. RESHAPE WIDE: one row per COD, columns per attained-age category
*    After reshape, variables are named obsexp1 smrstr1 aerstr1 ... obsexp6 etc.
*-------------------------------------------------------------------------------
keep  attagecat cod obsexp aerstr smrstr _order
order attagecat cod obsexp smrstr aerstr

// Grab the age-band levels AND their value labels before reshaping, so the
// labels can be used as column headings in the Word table below
levelsof attagecat, local(ages)
foreach a of local ages {
	local lbl`a' : label (attagecat) `a'
}
local nbands : word count `ages'

reshape wide obsexp smrstr aerstr, i(cod _order) j(attagecat)

// Sort rows to match the order defined in the cod local
sort _order
drop _order

// Build column list dynamically from the attained-age levels so the table
// adapts automatically if the age bands are changed in future
local vlist
foreach a of local ages {
	local vlist `vlist' obsexp`a' smrstr`a' aerstr`a'
}

// Count rows and columns (needed for borders and alignment below)
local nrows = _N
local ncols : word count `vlist'
local ncols = `ncols' + 1          // +1 for the cod column


*-------------------------------------------------------------------------------
* 5. EXPORT TO WORD TABLE (putdocx)
*    The table gets a two-row header: row 1 = attained-age band (spanning its
*    three columns), row 2 = Obs/Exp, SMR, AER.
*-------------------------------------------------------------------------------
cap putdocx clear
putdocx begin, font(arial narrow, "7") landscape
putdocx paragraph

// --- Title: "Table 5" in bold, remainder in normal weight ----
putdocx text ("Table 5"), bold
putdocx text (". SMRs and AERs for Specific Causes of Death by Attained Age")

local tnote "Obs/Exp = observed/expected deaths. SMR = standardised mortality ratio. AER = absolute excess risk per `aerperc' person-years."

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
	putdocx table tbl1(1,`c') = ("Obs/Exp")
	local ++c
	putdocx table tbl1(1,`c') = ("SMR (95% CI)")
	local ++c
	putdocx table tbl1(1,`c') = ("AER (95% CI)")
}

// --- Insert a spanner row above it holding the attained-age band labels ---
putdocx table tbl1(1,.), addrows(1, before)

// Fill the spanner right-to-left: merging cells removes the ones to the
// right, so working backwards keeps the column indices still to come valid
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
