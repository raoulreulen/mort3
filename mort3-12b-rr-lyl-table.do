/*=========================================================================================
mort3-12b-rr-lyl-table.do

PURPOSE
  Turn the life-years-lost (LYL) by cause results produced by R (lillies; Andersen
  decomposition) into a Word table with putdocx.
  Columns: N at risk at the index age; survivors' restricted remaining LE (lillies, 95% bootstrap CI);
  general-population restricted remaining LE (E&W 2022 period rates, sex-weighted);
  excess LYL = population LE - survivors' LE (95% CI); survivors' LYL before tau by cause of death
  (recurrence/progression, SPN, circulatory, respiratory, external, infection, other causes; 95% CI).
  Blocks: age 20 and age 40 (tau = 80) for overall / sex / cancer type;
          decade of diagnosis at age 20 and 30 with tau = 60 and tau = 80 (extrapolated for later decades).

INPUTS   $output/table_lyl_bycause.csv   (from mort3-12a-rr-lyl-bycause.R; first lines are '#' comments)
OUTPUTS  $output/table_lyl_bycause.docx

Date: 23 Oct 2026   Author: RCR (Claude-assisted)
=========================================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"

*-------------------------------------------------------------------------------
* Read CSV (skip the comment header lines: find the first row starting with "type")
*-------------------------------------------------------------------------------
tempname fh
local hdr = 0
local i = 0
file open `fh' using "$output/table_lyl_bycause.csv", read text
file read `fh' line
while r(eof)==0 {
    local ++i
    if substr(`"`macval(line)'"',1,6)==`""type""' & `hdr'==0 local hdr = `i'
    file read `fh' line
}
file close `fh'
di "variable names on row `hdr'"
import delimited "$output/table_lyl_bycause.csv", varnames(`hdr') rowrange(`hdr') clear ///
    stringcols(_all) bindquote(strict)
keep type stratum age tau n_atrisk extrapolated str_*
destring age tau n_atrisk extrapolated, replace

* Row order
gen order = .
local k = 0
foreach s in "all survivors" "male" "female" "leukaemia" "Hodgkin lymphoma" "non-Hodgkin lymphoma" "CNS" ///
             "neuroblastoma" "heritable retinoblastoma" "non-heritable retinoblastoma" "Wilms tumour" "bone" ///
             "soft tissue sarcoma" "other" "<1970" "1970-79" "1980-89" "1990-99" "2000-06" {
    local ++k
    replace order = `k' if stratum=="`s'"
}
assert order<.
format n_atrisk %9.0gc
gen str_n = string(n_atrisk, "%9.0fc")
* mark extrapolated rows (no survivor followed to tau)
replace str_excess = str_excess + "*" if extrapolated==1

*-------------------------------------------------------------------------------
* Assemble table body: blocks of rows with a block-header line
*-------------------------------------------------------------------------------
local blocks `""20 80 main" "40 80 main" "20 60 dec" "30 60 dec" "20 80 dec" "30 80 dec""'
foreach blk of local blocks {
    tokenize `blk'
    local a `1'
    local t `2'
    local set `3'
    preserve
    if "`set'"=="main" {
        keep if age==`a' & tau==`t' & type!="decade of diagnosis"
        local title "Age `a' years, maximum age (tau) `t': overall, by sex and by cancer type"
    }
    else {
        keep if age==`a' & tau==`t' & type=="decade of diagnosis"
        local title "Age `a' years, maximum age (tau) `t': by decade of diagnosis"
    }
    sort order
    keep stratum str_n str_le_surv str_le_pop str_excess str_recur str_spn str_circ str_resp str_external str_infection str_other
    gen rowtype = 0
    gen seq = _n
    set obs `=_N+1'
    replace rowtype = 1 in L
    replace stratum = "`title'" in L
    replace seq = 0 in L
    sort seq
    drop seq
    tempfile b_`a'_`t'_`set'
    save `b_`a'_`t'_`set'', replace
    restore
}
clear
foreach blk of local blocks {
    tokenize `blk'
    append using `b_`1'_`2'_`3''
}
foreach v of varlist str_* {
    replace `v' = "" if rowtype==1
}

label variable stratum       "Group"
label variable str_n         "N at risk"
label variable str_le_surv   "Survivors' remaining LE, years (95% CI)"
label variable str_le_pop    "Population remaining LE, years"
label variable str_excess    "Excess LYL, years (95% CI)"
label variable str_recur     "Recurrence/progression"
label variable str_spn       "SPN"
label variable str_circ      "Circulatory"
label variable str_resp      "Respiratory"
label variable str_external  "External"
label variable str_infection "Infection"
label variable str_other     "Other causes"

*-------------------------------------------------------------------------------
* Word table
*-------------------------------------------------------------------------------
cap: putdocx clear
putdocx begin, font(arial narrow, "10") landscape
putdocx paragraph
putdocx table tbl1 = data("stratum str_n str_le_surv str_le_pop str_excess str_recur str_spn str_circ str_resp str_external str_infection str_other"), ///
    varnames border(start, nil) border(insideV, nil) border(end, nil) width(100%) ///
    layout(autofitcontents) border(all, nil) ///
    title("Table. Life years lost (LYL) before the maximum age tau by 5-year survivors of childhood cancer, compared with the general population, and the survivors' LYL decomposed by cause of death (Andersen decomposition, R package lillies)") ///
    note("Survivors' remaining LE = area under the all-cause survival curve from the index age to tau (lillies; delayed entry at age of 5-year survival). Population remaining LE = same quantity from the England & Wales 2022 period life table, weighted by the sex distribution of survivors alive at the index age. Excess LYL = population LE minus survivors' LE = survivors' total LYL before tau minus population LYL before tau. Cause-specific columns = survivors' life years lost before tau to each cause (sum = tau minus index age minus survivors' LE; this total includes the years that the general population also loses). 95% CIs from 500 bootstrap resamples of survivors (percentile); the population life table is treated as fixed. 'Other causes' includes deaths of unknown cause. *No survivor in this group was followed to tau: the cumulative incidence is carried flat to tau, which understates LYL; later decades of diagnosis (1980 onwards) have not reached tau.", font(arial narrow, "8"))

* putdocx row indexing: row 1 = title, row 2 = column headers, data row r = table row r+2
putdocx table tbl1(1,.), bold
putdocx table tbl1(2,.), bold
local c = 0
foreach v of varlist stratum str_n str_le_surv str_le_pop str_excess str_recur str_spn str_circ str_resp str_external str_infection str_other {
    local ++c
    local lab : variable label `v'
    putdocx table tbl1(2,`c') = ("`lab'")
}
local N = _N
forvalues r = 1/`N' {
    if rowtype[`r']==1 {
        local rr = `r' + 2
        putdocx table tbl1(`rr',1), colspan(12)
        putdocx table tbl1(`rr',.), bold halign(left)
    }
}
putdocx table tbl1(.,1), halign(left)
forvalues c = 2/12 {
    putdocx table tbl1(.,`c'), halign(center)
}
putdocx save "$output/table_lyl_bycause.docx", replace
di "saved $output/table_lyl_bycause.docx"
exit
