/*==============================================================================
 WORD TABLES: LOSS IN EXPECTATION OF LIFE (LEL) FOR 5-YEAR SURVIVORS
 Imports the CSVs written by mort3-11a-rr-lel-stpm2.R and makes putdocx tables.
 All life-expectancy figures are REMAINING years of life from age a0; the
 "age reached" (= a0 + remaining LE) is shown alongside so that nobody reads
 remaining LE as age at death.

 Inputs : $output/table_lel.csv           (overall + by cancer type, a0 20/30/40/50)
          $output/table_lel_bydecade.csv  (by decade of diagnosis, a0 20/40)
          $output/table_lel_1990plus.csv  (era-specific model: 1990-2006 and 1980-89, a0 20/30/40)
 Outputs: $output/table_lel.docx
           Table A: rows = cancer type (overall first); columns = a0 20/30/40/50
           Table B: rows = decade of diagnosis, a0 = 20 and 40
           Table C: era-specific (diag x era model): 1990-2006 block then 1980-89 block,
                    a0 20/30/40, with extrapolation flag (E) and oldest attained age observed
 RCR (Claude-assisted), 23 Oct 2026
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
clear all

*-------------------------------------------------------------------------------
* helper: string cells, then wide by a0 (expects ord set and level string)
*-------------------------------------------------------------------------------
cap program drop lelwide
program define lelwide
    syntax , a0list(numlist) [flag]
    gen expstr  = string(exp_rem_le, "%4.1f") + " (" + string(exp_age_reached, "%4.1f") + ")"
    gen survstr = string(surv_rem_le, "%4.1f") + " (" + string(surv_rem_le_lo, "%4.1f") + " to " + string(surv_rem_le_hi, "%4.1f") + ")" ///
                  + " [" + string(surv_age_reached, "%4.1f") + "]"
    if "`flag'" != "" replace survstr = survstr + " E" if extrap_flag == "E"
    gen lelstr  = string(lel, "%4.1f") + " (" + string(lel_lo, "%4.1f") + " to " + string(lel_hi, "%4.1f") + ")"
    gen pctstr  = string(plost, "%4.1f") + " (" + string(plost_lo, "%4.1f") + " to " + string(plost_hi, "%4.1f") + ")"
    local extra
    if "`flag'" != "" local extra maxobs
    keep ord level a0 expstr survstr lelstr pctstr `extra'
    reshape wide expstr survstr lelstr pctstr, i(level `extra') j(a0)
    sort ord
    local vlist level
    foreach a of numlist `a0list' {
        local vlist `vlist' expstr`a' survstr`a' lelstr`a' pctstr`a'
    }
    order `vlist' `extra'
    keep `vlist' `extra'
end

cap program drop setord
program define setord
    // ord from a quoted list of level names (first = 1)
    gen ord = .
    local i = 0
    foreach l in `0' {
        local i = `i' + 1
        replace ord = `i' if level == "`l'"
    }
    assert ord != .
end

cap program drop lelheader
program define lelheader
    // header rows for table `tbl': first column label, 4 cells per a0, optional trailing column
    syntax , tbl(string) first(string) a0list(numlist) [trail(string)]
    putdocx table `tbl'(1, .), addrows(2, before)
    putdocx table `tbl'(1, 1) = ("`first'")
    local c = 2
    foreach a of numlist `a0list' {
        putdocx table `tbl'(1, `c') = ("Alive at age `a'"), colspan(4)
        local c = `c' + 1
    }
    if "`trail'" != "" putdocx table `tbl'(1, `c') = ("`trail'")
    local k = 0
    foreach a of numlist `a0list' {
        putdocx table `tbl'(2, `=2 + 4*`k'') = ("Expected remaining LE, general population (age reached)")
        putdocx table `tbl'(2, `=3 + 4*`k'') = ("Survivor remaining LE (95% CI) [age reached]")
        putdocx table `tbl'(2, `=4 + 4*`k'') = ("LEL, years (95% CI)")
        putdocx table `tbl'(2, `=5 + 4*`k'') = ("% of expected remaining life lost (95% CI)")
        local k = `k' + 1
    }
    local nc = 1 + 4 * `k' + ("`trail'" != "")
    putdocx table `tbl'(1, .), bold halign(center) border(top)
    putdocx table `tbl'(2, .), bold halign(center) border(bottom)
    putdocx table `tbl'(., 1), halign(left)
    forvalues c = 2/`nc' {
        putdocx table `tbl'(., `c'), halign(center)
    }
end

local diagorder overall leukaemia "Hodgkin lymphoma" "non-Hodgkin lymphoma" CNS neuroblastoma "heritable retinoblastoma" "non-heritable retinoblastoma" Wilms bone "soft tissue sarcoma" other
local lenote "LE = remaining life expectancy, i.e. remaining years of life from the stated age (not age at death); age reached = stated age + remaining LE. LEL = loss in expectation of life (expected minus survivor remaining LE). Expected remaining LE from the England & Wales 2022 period life table (sex-specific). Survivor remaining LE from a flexible parametric relative survival model (attained-age timescale; 5 df restricted cubic spline for the baseline log cumulative excess hazard; cancer type, sex and decade of diagnosis as covariates with proportional excess hazards; time-dependent cancer-type effects were tested and did not converge), integrating the product of expected survival and model relative survival to age 100 (the model is extrapolated beyond the oldest ages observed). 95% CIs from 500 parametric simulations of the model coefficients."

putdocx clear
putdocx begin, font("arial narrow", 10) landscape

*-------------------------------------------------------------------------------
* Table A: by cancer type
*-------------------------------------------------------------------------------
import delimited "$output/table_lel.csv", clear varnames(1) case(preserve) stringcols(2)
setord `diagorder'
lelwide, a0list(20 30 40 50)
replace level = "All survivors" if level == "overall"

putdocx paragraph
putdocx text ("Table A. Remaining life expectancy and loss in expectation of life among 5-year survivors of childhood cancer compared with the general population of England and Wales (2022 rates), conditional on being alive at age 20, 30, 40 and 50 years, by childhood cancer type"), bold
putdocx table tbl1 = data("level expstr20 survstr20 lelstr20 pctstr20 expstr30 survstr30 lelstr30 pctstr30 expstr40 survstr40 lelstr40 pctstr40 expstr50 survstr50 lelstr50 pctstr50"), border(all, nil) border(bottom)
lelheader, tbl(tbl1) first("Childhood cancer type") a0list(20 30 40 50)
putdocx table tbl1(3, .), bold
putdocx paragraph
putdocx text ("`lenote' Estimates are standardised to the sex and decade-of-diagnosis distribution of survivors under follow-up at the given age (overall: also cancer type)."), font("arial narrow", 9)

*-------------------------------------------------------------------------------
* Table B: by decade of diagnosis
*-------------------------------------------------------------------------------
import delimited "$output/table_lel_bydecade.csv", clear varnames(1) case(preserve) stringcols(2)
setord "<1970" 1970-79 1980-89 1990-99 2000-06
lelwide, a0list(20 40)

putdocx pagebreak
putdocx paragraph
putdocx text ("Table B. Remaining life expectancy and loss in expectation of life among 5-year survivors of childhood cancer by decade of diagnosis, conditional on being alive at age 20 and 40 years (standardised to the cancer type and sex distribution of each decade)"), bold
putdocx table tbl2 = data("level expstr20 survstr20 lelstr20 pctstr20 expstr40 survstr40 lelstr40 pctstr40"), border(all, nil) border(bottom)
lelheader, tbl(tbl2) first("Decade of diagnosis") a0list(20 40)
putdocx paragraph
putdocx text ("`lenote' Survivors diagnosed in 2000-06 had not reached age 40 by the end of follow-up (31 Dec 2022); their estimates at age 40 are model extrapolations assuming the decade effect on the excess hazard is proportional over age."), font("arial narrow", 9)

*-------------------------------------------------------------------------------
* Table C: era-specific (cancer type x era interaction model)
*-------------------------------------------------------------------------------
import delimited "$output/table_lel_1990plus.csv", clear varnames(1) case(preserve) stringcols(1 3 7)
local modelspec = model[1]
gen maxobs = string(max_age_obs, "%4.1f")
tempfile era
save `era'
foreach er in "1990-2006" "1980-89" {
    use `era', clear
    keep if era == "`er'"
    setord `diagorder'
    lelwide, a0list(20 30 40) flag
    replace level = "All survivors diagnosed `er'" if level == "overall"
    local nm = subinstr("`er'", "-", "_", .)
    tempfile blk_`nm'
    save `blk_`nm''
}
use `blk_1990_2006', clear
append using `blk_1980_89'

putdocx pagebreak
putdocx paragraph
putdocx text ("Table C. Remaining life expectancy and loss in expectation of life by childhood cancer type for survivors diagnosed 1990-2006 and 1980-89, conditional on being alive at age 20, 30 and 40 years (model with cancer type x era of diagnosis interaction: `modelspec')"), bold
putdocx table tbl3 = data("level expstr20 survstr20 lelstr20 pctstr20 expstr30 survstr30 lelstr30 pctstr30 expstr40 survstr40 lelstr40 pctstr40 maxobs"), border(all, nil) border(bottom)
lelheader, tbl(tbl3) first("Childhood cancer type") a0list(20 30 40) trail("Oldest attained age observed")
* block header rows: before row 3 (1990-2006 block) and before the 1980-89 block (12 data rows later)
putdocx table tbl3(3, .), addrows(1, before)
putdocx table tbl3(3, 1) = ("Diagnosed 1990-2006"), colspan(14) bold halign(left)
putdocx table tbl3(4, .), bold
putdocx table tbl3(16, .), addrows(1, before)
putdocx table tbl3(16, 1) = ("Diagnosed 1980-89"), colspan(14) bold halign(left)
putdocx table tbl3(17, .), bold
putdocx paragraph
putdocx text ("`lenote' Here the model additionally includes a cancer type x era of diagnosis interaction (eras <1980, 1980-89, 1990-2006, with decade main effects retained) so that the excess hazard ratio for each cancer type is era-specific; estimates are standardised to the sex (and, overall, cancer type) and decade distribution of survivors in that era under follow-up at the given age. E = extrapolation flag: the model-based mean age reached exceeds the oldest attained age observed for that era and cancer type (last column), i.e. the estimate relies on extrapolating the era-specific excess hazard beyond the ages at which survivors of that era have been observed; for survivors diagnosed 1990-2006 this applies to every cell."), font("arial narrow", 9)

putdocx save "$output/table_lel.docx", replace
exit
