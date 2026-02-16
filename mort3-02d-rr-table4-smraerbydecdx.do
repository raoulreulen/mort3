/*==============================================================================
TABLE 4: SMRs and AERs by Decade of Diagnosis (DECDX)

PURPOSE
	For each cause of death (COD), this do-file:
	  1. Loads the survival-time (stset) dataset created in do-file 02a
	  2. Calculates expected deaths (_e) and person-years (_y)
	  3. Collapses data by decade-of-diagnosis category (decdxcat)
	  4. Stacks all CODs into one dataset
	  5. Computes SMRs and AERs  (using the custom -smraer- command)
	  6. Reshapes the data wide (one row per COD, columns per decade)
	  7. Exports a formatted Word table via putdocx

REQUIRES
	- Globals: $temp  (path to temp/output folder)
	- Datasets: $temp/x-mort3-stset-<cod>.dta  (from do-file 02a)
	- Ado-file: smraer.ado
==============================================================================*/


*-------------------------------------------------------------------------------
* 1. DEFINE CAUSES OF DEATH
*    Each name matches a stset dataset and a rate variable (rate_<cod>).
*-------------------------------------------------------------------------------
local cod allcauses neoplasm spn infection blood endocrine mental nervous ///
	circulation respiratory digestive muscoskeletal genitourinary perinatal ///
	other external suicide


*-------------------------------------------------------------------------------
* 2. LOOP OVER EACH COD: CALCULATE EXPECTED DEATHS & PERSON-YEARS
*    Then stack all CODs into a single temporary file.
*-------------------------------------------------------------------------------
tempfile stacked
local i = 0
local n = 0

foreach x of local cod {
	local ++n

	// Load the stset dataset for this COD
	use "$temp/x-mort3-stset-`x'", clear

	// _e = expected deaths (individual person-time * expected rate)
	// _y = person-years in units of 10,000 (for AER scaling)
	gen _e = rate_`x' * (_t - _t0)
	gen _y = (_t - _t0) / 100000

	// Sum observed deaths (_d), person-years (_y), expected deaths (_e)
	// within each decade-of-diagnosis category
	collapse (sum) _d _y _e, by(decdxcat)

	// Tag rows with the COD name and preserve original order
	gen str cod = "`x'"
	gen _order = `n'

	// Stack: append previous iterations then save
	if `i'>0 append using `stacked', force
	save `stacked', replace
	local i = 1
}


*-------------------------------------------------------------------------------
* 3. COMPUTE SMRs AND AERs
*    smraer creates: smrstr (SMR + 95% CI), aerstr (AER + 95% CI), obsexp
*    desaer(0) = 0 decimal places for AER; dessmr(1) = 1 decimal for SMR.
*-------------------------------------------------------------------------------
smraer _d _e _y, desaer(0) dessmr(1)


*-------------------------------------------------------------------------------
* 4. RESHAPE WIDE: one row per COD, columns per decade category
*    After reshape, variables are named obsexp1 smrstr1 aerstr1 ... obsexp5 etc.
*-------------------------------------------------------------------------------
keep  decdxcat cod obsexp aerstr smrstr _order
order decdxcat cod obsexp smrstr aerstr

levelsof decdxcat, local(dec)
reshape wide obsexp smrstr aerstr, i(cod _order) j(decdxcat)

// Sort rows to match the order defined in the cod local
sort _order
drop _order

// Build column list dynamically from decade levels so the table
// adapts automatically if decades are added/removed in future
local vlist
foreach d of local dec {
	local vlist `vlist' obsexp`d' smrstr`d' aerstr`d'
}

// Count rows and columns (needed for borders and alignment below)
local nrows = _N
local ncols : word count `vlist'
local ncols = `ncols' + 1          // +1 for the cod column


*-------------------------------------------------------------------------------
* 5. EXPORT TO WORD TABLE (putdocx)
*-------------------------------------------------------------------------------
cap putdocx clear
putdocx begin, font(arial narrow, "8") landscape
putdocx paragraph

// --- Title: "Table 4" in bold, remainder in normal weight ----
putdocx text ("Table 4"), bold
putdocx text (". SMRs and AERs for Specific Causes of Death by Decade of Childhood Cancer Diagnosis")

// --- Create the table from data ---
// border(all, nil) removes all borders; we add specific ones below
// +2 rows: row 1 = variable labels (header), rows 2..N+1 = data
putdocx table tbl1 = data("cod `vlist'"), varnames ///
	border(start, nil) border(insideV, nil) border(end, nil) ///
	border(insideH, nil) border(all, nil) ///
	width(100%) layout(autofitcontents)

// --- Header row: bold and centred ---
putdocx table tbl1(1,.), bold

// --- First column (COD names): left-aligned and bold ---
putdocx table tbl1(.,1), halign(left) bold

// --- Remaining columns: centre-aligned ---
forvalues c = 2/`ncols' {
	putdocx table tbl1(.,`c'), halign(center)
}

// --- Horizontal line below header row (bottom border of row 1) ---
putdocx table tbl1(1,.), border(bottom, single)

// --- Horizontal line at the bottom of the table (bottom border of last row) ---
// Data rows run from 2 to nrows+1 (row 1 is header)
local lastrow = `nrows' + 1
putdocx table tbl1(`lastrow',.), border(bottom, single)

// --- Save and open the document ---
putdocx save "$temp/table4_smraerbydecdx_cod3paper.docx", replace
shell "$temp/table4_smraerbydecdx_cod3paper.docx"

exit
