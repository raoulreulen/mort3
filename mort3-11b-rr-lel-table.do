/*==============================================================================
 WORD TABLES: LOSS IN EXPECTATION OF LIFE (LEL) FOR 5-YEAR SURVIVORS
 Imports the CSVs written by mort3-11a-rr-lel-stpm2.R and makes putdocx tables.

 Inputs : $output/table_lel.csv           (overall + by cancer type, a0 20/30/40/50)
          $output/table_lel_bydecade.csv  (by decade of diagnosis, a0 20/40)
 Outputs: $output/table_lel.docx  - Table A: rows = cancer type (overall first);
                                    columns = a0 20/30/40/50, each with expected LE,
                                    survivor LE, LEL (95% CI), % of expected life lost
                                  - Table B: rows = decade of diagnosis, a0 = 20 and 40
 RCR (Claude-assisted), 23 Oct 2026
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
clear all

cap program drop lelwide
program define lelwide
    // reshape long (group level a0 ...) to one row per level with string cells per a0
    syntax , a0list(numlist)
    gen expstr  = string(exp_le, "%4.1f")
    gen survstr = string(surv_le, "%4.1f") + " (" + string(surv_le_lo, "%4.1f") + " to " + string(surv_le_hi, "%4.1f") + ")"
    gen lelstr  = string(lel, "%4.1f") + " (" + string(lel_lo, "%4.1f") + " to " + string(lel_hi, "%4.1f") + ")"
    gen pctstr  = string(plost, "%4.1f") + " (" + string(plost_lo, "%4.1f") + " to " + string(plost_hi, "%4.1f") + ")"
    keep ord level a0 expstr survstr lelstr pctstr
    reshape wide expstr survstr lelstr pctstr, i(level) j(a0)
    sort ord
    drop ord
    local vlist level
    foreach a of numlist `a0list' {
        local vlist `vlist' expstr`a' survstr`a' lelstr`a' pctstr`a'
    }
    order `vlist'
    keep `vlist'
end

*-------------------------------------------------------------------------------
* Table A: by cancer type
*-------------------------------------------------------------------------------
import delimited "$output/table_lel.csv", clear varnames(1) case(preserve) stringcols(2)
gen ord = .
local i = 0
foreach l in overall leukaemia "Hodgkin lymphoma" "non-Hodgkin lymphoma" CNS neuroblastoma "heritable retinoblastoma" "non-heritable retinoblastoma" Wilms bone "soft tissue sarcoma" other {
    local i = `i' + 1
    replace ord = `i' if level == "`l'"
}
assert ord != .
lelwide, a0list(20 30 40 50)
replace level = "All survivors" if level == "overall"

putdocx clear
putdocx begin, font("arial narrow", 10) landscape
putdocx paragraph
putdocx text ("Table. Residual life expectancy and loss in expectation of life among 5-year survivors of childhood cancer compared with the general population of England and Wales (2022 rates), conditional on being alive at age 20, 30, 40 and 50 years, by childhood cancer type"), bold
putdocx table tbl1 = data("level expstr20 survstr20 lelstr20 pctstr20 expstr30 survstr30 lelstr30 pctstr30 expstr40 survstr40 lelstr40 pctstr40 expstr50 survstr50 lelstr50 pctstr50"), border(all, nil) border(bottom)
putdocx table tbl1(1, .), addrows(2, before)
putdocx table tbl1(1, 1) = ("Childhood cancer type")
putdocx table tbl1(1, 2) = ("Alive at age 20"), colspan(4)
putdocx table tbl1(1, 3) = ("Alive at age 30"), colspan(4)
putdocx table tbl1(1, 4) = ("Alive at age 40"), colspan(4)
putdocx table tbl1(1, 5) = ("Alive at age 50"), colspan(4)
forvalues k = 0/3 {
    putdocx table tbl1(2, `=2 + 4*`k'') = ("Expected LE")
    putdocx table tbl1(2, `=3 + 4*`k'') = ("Survivor LE (95% CI)")
    putdocx table tbl1(2, `=4 + 4*`k'') = ("LEL, years (95% CI)")
    putdocx table tbl1(2, `=5 + 4*`k'') = ("% of expected life lost (95% CI)")
}
putdocx table tbl1(1, .), bold halign(center) border(top)
putdocx table tbl1(2, .), bold halign(center) border(bottom)
putdocx table tbl1(., 1), halign(left)
forvalues c = 2/17 {
    putdocx table tbl1(., `c'), halign(center)
}
putdocx table tbl1(3, .), bold
putdocx paragraph
putdocx text ("LE = residual life expectancy in years; LEL = loss in expectation of life (expected minus survivor LE). Expected LE from the England & Wales 2022 period life table (sex-specific). Survivor LE from a flexible parametric relative survival model (attained-age timescale; 5 df restricted cubic spline for the baseline log cumulative excess hazard; cancer type, sex and decade of diagnosis as covariates), integrating the product of expected survival and model relative survival to age 100 (the model is extrapolated beyond the oldest ages observed). Estimates are standardised to the sex and decade-of-diagnosis distribution of survivors under follow-up at the given age (overall: also cancer type). 95% CIs from 500 parametric simulations of the model coefficients."), font("arial narrow", 9)

*-------------------------------------------------------------------------------
* Table B: by decade of diagnosis
*-------------------------------------------------------------------------------
import delimited "$output/table_lel_bydecade.csv", clear varnames(1) case(preserve) stringcols(2)
gen ord = .
local i = 0
foreach l in "<1970" 1970-79 1980-89 1990-99 2000-06 {
    local i = `i' + 1
    replace ord = `i' if level == "`l'"
}
assert ord != .
lelwide, a0list(20 40)

putdocx pagebreak
putdocx paragraph
putdocx text ("Table. Residual life expectancy and loss in expectation of life among 5-year survivors of childhood cancer by decade of diagnosis, conditional on being alive at age 20 and 40 years (standardised to the cancer type and sex distribution of each decade)"), bold
putdocx table tbl2 = data("level expstr20 survstr20 lelstr20 pctstr20 expstr40 survstr40 lelstr40 pctstr40"), border(all, nil) border(bottom)
putdocx table tbl2(1, .), addrows(2, before)
putdocx table tbl2(1, 1) = ("Decade of diagnosis")
putdocx table tbl2(1, 2) = ("Alive at age 20"), colspan(4)
putdocx table tbl2(1, 3) = ("Alive at age 40"), colspan(4)
forvalues k = 0/1 {
    putdocx table tbl2(2, `=2 + 4*`k'') = ("Expected LE")
    putdocx table tbl2(2, `=3 + 4*`k'') = ("Survivor LE (95% CI)")
    putdocx table tbl2(2, `=4 + 4*`k'') = ("LEL, years (95% CI)")
    putdocx table tbl2(2, `=5 + 4*`k'') = ("% of expected life lost (95% CI)")
}
putdocx table tbl2(1, .), bold halign(center) border(top)
putdocx table tbl2(2, .), bold halign(center) border(bottom)
putdocx table tbl2(., 1), halign(left)
forvalues c = 2/9 {
    putdocx table tbl2(., `c'), halign(center)
}
putdocx paragraph
putdocx text ("Survivors diagnosed in 2000-06 had not reached age 40 by the end of follow-up (31 Dec 2022); their estimates at age 40 are model extrapolations assuming the decade effect on the excess hazard is proportional over age."), font("arial narrow", 9)

putdocx save "$output/table_lel.docx", replace
exit
