/*==============================================================================
FIGURE 3: CUMULATIVE MORTALITY BY DECADE OF DIAGNOSIS (DECDX)

PURPOSE
	For each cause of death (COD), this do-file:
	  1. Runs stcompet separately for each decdxcat level
	  2. Produces a cumulative incidence figure with one line per decade

REQUIRES
	- Globals: $temp
	- Datasets: $temp/x-mort3-prepforstset
	- Ado-files: stcompet
==============================================================================*/

local cod spn

local xtime 70

// Colours to cycle through for each decade (add more if needed)
local colours black blue red green orange purple



*-------------------------------------------------------------------------------
* LOOP OVER EACH COD — one figure per COD
*-------------------------------------------------------------------------------
foreach x of local cod {

	// --- 2a. Load once, recode, stset -----------------------------------
	use "$temp/x-mort3-prepforstset", clear

	// Other COD coded as 2 = competing risk
	recode `x' (0=2) if allcauses == 1

	// stset on this cause of death
	stset dox, fail(`x'==1) id(indexno) origin(dob) entry(doe) scale(365.25)
	assert _st == 1
	
	// Extract value labels from decdxcat for legend
	levelsof decdxcat, local(decades)
	foreach d of local decades {
		local lbl`d' : label (decdxcat) `d'
	}

	// Log-rank test for heterogeneity across decades
	sts test decdxcat
	local chi2 = r(chi2)
	local df   = r(df)
	local pval = chi2tail(`df', `chi2')
	if `pval' < 0.001 {
		local pvaltxt "P<0.001"
	}
	else {
		local pvaltxt "P=`: di %4.2f `pval''"
	}

	// --- 2b. Run stcompet for each decade category ----------------------
	tempfile combined
	local first = 1

	foreach d of local decades {

		// Restore the master data each iteration (avoid re-reading from disk)
		preserve
		keep if decdxcat == `d'

		// Cumulative incidence with competing risks
		stcompet cuminc = ci, compet1(2)
		gen ci_`d' = cuminc * 100 if _d == 1
		keep if _d == 1
		keep ci_`d' _t
		sort _t

		// Stack all decades into one file
		if `first' {
			save `combined', replace
			local first = 0
		}
		else {
			append using `combined', force
			save `combined', replace
		}
		restore
	}

	// --- 2c. Build figure options dynamically from decade levels --------
	use `combined', clear

	// Dynamic y-axis: find max across all CI curves, round up
	egen _cimax = rowmax(ci_*)
	quietly summarize _cimax
	local ymax = ceil(r(max)) + 1
	if `ymax' <= 20      local ytick = 1
	else if `ymax' <= 40 local ytick = 2
	else                 local ytick = 5
	drop _cimax

	local civars
	local conn
	local clps
	local clws
	local clcs
	local legorder
	local k = 0

	foreach d of local decades {
		local ++k
		local civars   `civars' ci_`d'
		local conn     `conn' J
		local clps     `clps' solid
		local clws     `clws' med
		local col : word `k' of `colours'
		local clcs     `clcs' `col'
		local legorder `legorder' `k' "`lbl`d''"
	}


	// --- 2d. Draw figure ------------------------------------------------
	#delimit ;

	twoway

	/* CUMULATIVE INCIDENCE BY DECADE */
	(line `civars' _t if (_t>=10 & _t<=`xtime')
	, sort connect(`conn') clp(`clps') clw(`clws') clc(`clcs'))

		,

		/*LABELS*/
		xtitle("Attained age, years", size(small))
		ytitle("Cumulative mortality, %", size(small))
		ylabel(0(`ytick')`ymax', angle(0)) xlabel(10(5)`xtime')
		title("Deaths due to `x' by decade of childhood cancer diagnosis", size(small))

		/*LEGEND with p-value below last entry*/
		legend(on order(`legorder')
		note("`pvaltxt' (log-rank)", size(small))
		ring(0) position(11) size(small) rowgap(0.1) cols(1))

	;
	#delimit cr

}
