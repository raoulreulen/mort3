/*==============================================================================
 mort3-10b-rr-sex-smraer.do
 PURPOSE
   SMRs and AERs (per 100,000 person-years) by sex for every cause of death,
   with a p for heterogeneity between sexes in the SMR (Poisson model of
   observed deaths with log(expected) offset, Wald test of i.sex) and a p for
   the difference in AER (normal approximation: diff / sqrt(se1^2+se2^2),
   se(AER) = sqrt(O)/PY).
 INPUTS   $temp/x-mort3-stset-<cod>.dta  (do-file 02a); smraer.ado
 OUTPUTS  $output/table_sex_smraer.docx ; $output/table_sex_smraer.csv
 NOTES    recur has no population comparator: SMR/AER blank; the p for
          heterogeneity is for the crude death rate (offset log(person-years)).
 Date 23 Oct 2026   Author RCR (Claude-assisted)
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
version 18
set more off

local cod allcauses allcodexcrecur recur spn infection blood endocrine mental nervous ///
	circulation cardiac cerebrovascular othercirculatory respiratory digestive ///
	muscoskeletal genitourinary perinatal external suicide other

local lab_allcauses       "all causes"
local lab_allcodexcrecur  "all causes excl. recurrence"
local lab_recur           "recurrence/progression"
local lab_spn             "SPN"
local lab_infection       "infection"
local lab_blood           "blood"
local lab_endocrine       "endocrine"
local lab_mental          "mental"
local lab_nervous         "nervous system"
local lab_circulation     "circulatory"
local lab_cardiac         "cardiac"
local lab_cerebrovascular "cerebrovascular"
local lab_othercirculatory "other circulatory"
local lab_respiratory     "respiratory"
local lab_digestive       "digestive"
local lab_muscoskeletal   "musculoskeletal"
local lab_genitourinary   "genitourinary"
local lab_perinatal       "perinatal/congenital"
local lab_external        "external causes"
local lab_suicide         "suicide"
local lab_other           "other"

*-------------------------------------------------------------------------------
* 1. Collapse each cause by sex and stack
*-------------------------------------------------------------------------------
tempfile stacked
local i = 0
local n = 0
foreach x of local cod {
	local ++n
	use _d _t _t0 newrate sex using "$temp/x-mort3-stset-`x'", clear
	gen double _e = (_t-_t0)*newrate
	gen double _y = (_t-_t0)/100000
	collapse (sum) _d _e _y, by(sex)
	gen str30 cod = "`lab_`x''"
	gen str20 codvar = "`x'"
	gen _order = `n'
	if `i'>0 append using `stacked'
	save `stacked', replace
	local i = 1
}
use `stacked', clear
sort _order sex

*-------------------------------------------------------------------------------
* 2. SMR / AER by sex (smraer.ado), plus p-values
*-------------------------------------------------------------------------------
smraer _d _e _y, dessmr(1) desaer(0)
* recur: no comparator -> blank SMR/AER
replace smrstr = "" if codvar=="recur"
replace aerstr = "" if codvar=="recur"
replace obsexp = strobs + "/-" if codvar=="recur"

* se of AER (per 100,000 PY) = sqrt(O)/PY
gen double aerse = sqrt(_d)/_y

gen str10 phet = ""
gen str10 paer = ""
levelsof _order, local(ords)
foreach o of local ords {
	qui levelsof codvar if _order==`o', local(cv)
	local cv : word 1 of `cv'
	preserve
		keep if _order==`o'
		* p for heterogeneity in SMR (crude rate for recur)
		if "`cv'"=="recur" {
			qui poisson _d i.sex, exposure(_y)
		}
		else {
			qui poisson _d i.sex, exposure(_e)
		}
		qui test 2.sex
		local ph = r(p)
		* p for AER difference
		qui sum aer if sex==1
		local a1 = r(mean)
		qui sum aerse if sex==1
		local s1 = r(mean)
		qui sum aer if sex==2
		local a2 = r(mean)
		qui sum aerse if sex==2
		local s2 = r(mean)
		local z = (`a1'-`a2')/sqrt(`s1'^2+`s2'^2)
		local pa = 2*(1-normal(abs(`z')))
		di "`cv': p-het SMR = `ph'  AER male = `a1' female = `a2'  z = `z'  p = `pa'"
	restore
	replace phet = cond(`ph'<0.001,"<0.001",string(`ph',"%9.3f")) if _order==`o'
	if "`cv'"!="recur" replace paer = cond(`pa'<0.001,"<0.001",string(`pa',"%9.3f")) if _order==`o'
}

*-------------------------------------------------------------------------------
* 3. Reshape wide: one row per cause
*-------------------------------------------------------------------------------
keep cod codvar _order sex _d _e _y obsexp smrstr aerstr phet paer
rename (_d _e _y) (obs exp py)
reshape wide obs exp py obsexp smrstr aerstr, i(cod codvar _order phet paer) j(sex)
sort _order
order cod obsexp1 smrstr1 aerstr1 obsexp2 smrstr2 aerstr2 phet paer
label var cod     "Cause of death"
label var obsexp1 "Males O/E"
label var smrstr1 "Males SMR (95% CI)"
label var aerstr1 "Males AER (95% CI)*"
label var obsexp2 "Females O/E"
label var smrstr2 "Females SMR (95% CI)"
label var aerstr2 "Females AER (95% CI)*"
label var phet    "p-het SMR"
label var paer    "p-diff AER"
export delimited cod codvar obs1 exp1 py1 obsexp1 smrstr1 aerstr1 obs2 exp2 py2 obsexp2 smrstr2 aerstr2 phet paer ///
	using "$output/table_sex_smraer.csv", replace
local nrows = _N

*-------------------------------------------------------------------------------
* 4. Word table
*-------------------------------------------------------------------------------
cap putdocx clear
putdocx begin, font("arial narrow", 10) landscape
putdocx paragraph
putdocx text ("Table. "), bold
putdocx text ("Standardised mortality ratios and absolute excess risks for specific causes of death, by sex")
putdocx table tbl1 = data("cod obsexp1 smrstr1 aerstr1 obsexp2 smrstr2 aerstr2 phet paer"), varnames ///
	border(start, nil) border(insideV, nil) border(end, nil) border(insideH, nil) border(all, nil) ///
	width(100%) layout(autofitcontents)
putdocx table tbl1(1,.), bold border(bottom, single)
putdocx table tbl1(.,1), halign(left) bold
forvalues c = 2/9 {
	putdocx table tbl1(.,`c'), halign(center)
}
local last = `nrows'+1
putdocx table tbl1(`last',.), border(bottom, single)
putdocx paragraph
putdocx text ("O/E = observed/expected deaths; SMR with exact Poisson 95% CI; *AER per 100,000 person-years. p-het SMR: Wald test of sex in a Poisson model of observed deaths with log(expected) as offset (for recurrence/progression, which has no population comparator, the test compares crude death rates and SMR/AER are not shown). p-diff AER: normal approximation, (AER males - AER females)/sqrt(se1^2+se2^2) with se = sqrt(O)/PY."), font("arial narrow", 8)
putdocx save "$output/table_sex_smraer.docx", replace
di "saved $output/table_sex_smraer.docx"
exit
