/*==============================================================================
 MORTALITY AGE GAP: "a survivor aged a has the all-cause mortality of a general
 population person aged a + gap"

 Inputs : $temp/x-mort3-stset-allcauses.dta  (stset/stsplit all-cause data,
          _t0/_t attained age, _d death, sex, diag)
          $rates/allcauses.dta  (E&W all-cause rates; 2022 used as reference)
 Outputs: $output/fig_mortality_age_gap.png      (gap vs age, by sex, 95% CI)
          $output/table_mortality_age_gap.csv    (age 20,30,40,50,60 x sex)
          $output/table_mortality_age_gap_bydiag.csv (leukaemia, CNS, Hodgkin,
                                                  Wilms; ages 30 and 50, by sex)
          $output/fig_mortality_age_gap_data.csv (all ages 20-65, numbers behind fig)
 Method : observed survivor all-cause rate by single year of attained age
          (floor of record midpoint) and sex, smoothed with a Poisson model with
          restricted cubic spline of age (5 knots, Harrell percentiles) on ages
          20-70. General population 2022 E&W rates by sex: 5-y band rates
          interpolated log-linearly between band midpoints (85+ placed at 90,
          log-linear extrapolation beyond). For survivor age a (exact age a+0.5)
          the population age a* with the same rate is found by linear
          interpolation on the log-rate scale; gap = a* - a. 95% percentile CI
          from 200 bootstrap replicates (Poisson resampling of deaths in the
          collapsed age-sex table; population rates fixed).
 RCR (Claude-assisted), 23 Oct 2026
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
clear all
set seed 20261023
local NBOOT 200
local AGELO 20
local AGEHI 70
local NKNOT 5

*-------------------------------------------------------------------------------
* 1. Population rates 2022, single ages (exact age x = 0.5, 1.5, ..., 100.5)
*-------------------------------------------------------------------------------
use "$rates/allcauses", clear
keep if yeargrp == 2022
keep sex ageband newrate
gen mid = cond(ageband == 85, 90, ageband + 2.5)
gen lr = log(newrate)
sort sex mid
tempfile pop
save `pop'

mata:
mata clear
// log-linear interpolation of band midpoints to exact ages x; extrapolate with last slope
real colvector poplr(real colvector mid, real colvector lr, real colvector x)
{
    real colvector out
    real scalar i, j, n
    n = rows(mid)
    out = J(rows(x), 1, .)
    for (i = 1; i <= rows(x); i++) {
        if (x[i] <= mid[1]) {
            out[i] = lr[1]
        }
        else if (x[i] >= mid[n]) {
            out[i] = lr[n] + (lr[n] - lr[n-1]) / (mid[n] - mid[n-1]) * (x[i] - mid[n])
        }
        else {
            for (j = 1; j < n; j++) {
                if (x[i] >= mid[j] & x[i] < mid[j+1]) {
                    out[i] = lr[j] + (lr[j+1] - lr[j]) / (mid[j+1] - mid[j]) * (x[i] - mid[j])
                }
            }
        }
    }
    return(out)
}
// equivalent population age: age x at which population log rate (monotone from age>=10) equals target
real colvector eqage(real colvector x, real colvector plr, real colvector target)
{
    real colvector out
    real scalar i, j, n, start
    n = rows(x)
    out = J(rows(target), 1, .)
    start = 1
    for (j = 1; j <= n; j++) {
        if (x[j] >= 10) {
            start = j
            j = n
        }
    }
    for (i = 1; i <= rows(target); i++) {
        if (target[i] == .) continue
        if (target[i] < plr[start]) {
            out[i] = x[start]
            continue
        }
        for (j = start; j < n; j++) {
            if (target[i] >= plr[j] & target[i] < plr[j+1]) {
                out[i] = x[j] + (target[i] - plr[j]) / (plr[j+1] - plr[j]) * (x[j+1] - x[j])
            }
        }
        if (out[i] == .) out[i] = x[n] + (target[i] - plr[n]) / (plr[n] - plr[n-1]) * (x[n] - x[n-1])
    }
    return(out)
}
// percentile of a row vector ignoring missing
real scalar pctl(real rowvector v, real scalar p)
{
    real colvector w
    w = sort(select(v', v' :!= .), 1)
    if (rows(w) == 0) return(.)
    return(w[max((1, ceil(p * rows(w))))])
}
// assemble OUT: age, surv rate (lo, hi), pop rate at age, equivalent age (lo, hi), gap (lo, hi)
void gapcalc(string scalar pmidv, string scalar plrv, string scalar selv)
{
    real colvector pmid, plr, xx, plx, xa, lr0, eq0, popr, lo, hi, rlo, rhi, age
    real matrix LRB, EQB
    real scalar b, i
    pmid = st_data(., pmidv, selv)
    plr  = st_data(., plrv, selv)
    xa  = st_matrix("XA")
    lr0 = st_matrix("LR0")
    LRB = st_matrix("LRB")
    xx = (0.5::100.5)
    plx = poplr(pmid, plr, xx)
    eq0 = eqage(xx, plx, lr0)
    popr = exp(poplr(pmid, plr, xa))
    EQB = J(rows(lr0), cols(LRB), .)
    for (b = 1; b <= cols(LRB); b++) EQB[., b] = eqage(xx, plx, LRB[., b])
    lo = J(rows(lr0), 1, .); hi = lo; rlo = lo; rhi = lo
    for (i = 1; i <= rows(lr0); i++) {
        lo[i]  = pctl(EQB[i, .], 0.025)
        hi[i]  = pctl(EQB[i, .], 0.975)
        rlo[i] = exp(pctl(LRB[i, .], 0.025))
        rhi[i] = exp(pctl(LRB[i, .], 0.975))
    }
    age = xa :- 0.5
    st_matrix("OUT", (age, exp(lr0), rlo, rhi, popr, eq0, lo, hi, eq0 - age, lo - age, hi - age))
}
end

*-------------------------------------------------------------------------------
* 2. Survivor deaths and person-years by single year of attained age, sex, diag
*-------------------------------------------------------------------------------
use sex diag _t0 _t _d using "$temp/x-mort3-stset-allcauses", clear
gen py = _t - _t0
gen agemid = floor((_t0 + _t) / 2)
collapse (sum) _d py, by(sex diag agemid)
tempfile coll
save `coll'

*-------------------------------------------------------------------------------
* 3. Program agegap2: smooth survivor rates for one sex (optional diag subset),
*    equivalent population age and gap at requested ages, bootstrap CI;
*    results appended to frame res
*-------------------------------------------------------------------------------
frame create pop
frame pop: use `pop'
frame pop: gen byte s1 = sex == 1
frame pop: gen byte s2 = sex == 2
frame create coll
frame coll: use `coll'
frame create res
frame res {
    set obs 0
    gen str20 group = ""
    gen sex = .
    gen age = .
    gen surv_rate = .
    gen surv_rate_lo = .
    gen surv_rate_hi = .
    gen pop_rate = .
    gen eq_age = .
    gen eq_age_lo = .
    gen eq_age_hi = .
    gen gap = .
    gen gap_lo = .
    gen gap_hi = .
    gen deaths = .
    gen py = .
}

cap program drop agegap2
program define agegap2
    // computes gap at ages `ages' for sex `sex' (diag subset optional), appends to frame res
    syntax , sex(integer) ages(numlist) group(string) [diag(string) nknots(integer 5) nboot(integer 200)]
    frame copy coll work, replace
    frame work {
        keep if sex == `sex'
        if "`diag'" != "" keep if inlist(diag, `diag')
        collapse (sum) _d py, by(agemid)
        keep if agemid >= `=$AGELO' & agemid < `=$AGEHI'
        gen x = agemid + 0.5
        qui summ _d
        local totd = r(sum)
        qui summ py
        local totpy = r(sum)
        // prediction rows for requested ages
        local nag : word count `ages'
        local n0 = _N
        set obs `=_N + `nag''
        local k = 0
        foreach a of numlist `ages' {
            local k = `k' + 1
            replace x = `a' + 0.5 in `=`n0' + `k''
        }
        gen byte pred = _n > `n0'
        // knots from the data rows only (Harrell percentiles), then splines for all rows
        mkspline xs = x if pred == 0, cubic nknots(`nknots')
        matrix K = r(knots)
        local knots
        forvalues j = 1/`=colsof(K)' {
            local knots `knots' `=K[1, `j']'
        }
        drop xs*
        mkspline xs = x, cubic knots(`knots')
        gen d0 = _d
        gen lnpy = log(py)
        // point estimate
        qui poisson _d xs* if pred == 0, offset(lnpy)
        predict double xb0 if pred == 1, xb nooffset
        mkmat x if pred == 1, matrix(XA)
        mkmat xb0 if pred == 1, matrix(LR0)
        // bootstrap
        tempname LRB
        matrix `LRB' = J(`nag', `nboot', .)
        forvalues b = 1/`nboot' {
            qui replace _d = cond(d0 > 0, rpoisson(d0), 0) if pred == 0
            cap qui poisson _d xs* if pred == 0, offset(lnpy) iterate(50)
            if _rc == 0 {
                tempvar xbb
                qui predict double `xbb' if pred == 1, xb nooffset
                mkmat `xbb' if pred == 1, matrix(TMP)
                matrix `LRB'[1, `b'] = TMP
            }
        }
        matrix LRB = `LRB'
    }
    frame pop: mata: gapcalc("mid", "lr", "s`sex'")
    local nag = rowsof(OUT)
    frame res {
        local n0 = _N
        set obs `=_N + `nag''
        forvalues i = 1/`nag' {
            local r = `n0' + `i'
            replace group = "`group'" in `r'
            replace sex = `sex' in `r'
            replace age = OUT[`i', 1] in `r'
            replace surv_rate = OUT[`i', 2] in `r'
            replace surv_rate_lo = OUT[`i', 3] in `r'
            replace surv_rate_hi = OUT[`i', 4] in `r'
            replace pop_rate = OUT[`i', 5] in `r'
            replace eq_age = OUT[`i', 6] in `r'
            replace eq_age_lo = OUT[`i', 7] in `r'
            replace eq_age_hi = OUT[`i', 8] in `r'
            replace gap = OUT[`i', 9] in `r'
            replace gap_lo = OUT[`i', 10] in `r'
            replace gap_hi = OUT[`i', 11] in `r'
            replace deaths = `totd' in `r'
            replace py = `totpy' in `r'
        }
    }
    frame drop work
end

global AGELO `AGELO'
global AGEHI `AGEHI'

*-------------------------------------------------------------------------------
* 4. All survivors, ages 20-65 by sex
*-------------------------------------------------------------------------------
foreach s in 1 2 {
    agegap2, sex(`s') ages(20(1)65) group("all") nknots(`NKNOT') nboot(`NBOOT')
}

*-------------------------------------------------------------------------------
* 5. Four largest groups at ages 30 and 50 (4 knots: sparser data)
*-------------------------------------------------------------------------------
foreach s in 1 2 {
    agegap2, sex(`s') ages(30 50) group("leukaemia") diag(1) nknots(4) nboot(`NBOOT')
    agegap2, sex(`s') ages(30 50) group("CNS")       diag(4) nknots(4) nboot(`NBOOT')
    agegap2, sex(`s') ages(30 50) group("Hodgkin")   diag(2) nknots(4) nboot(`NBOOT')
    agegap2, sex(`s') ages(30 50) group("Wilms")     diag(8) nknots(4) nboot(`NBOOT')
}

*-------------------------------------------------------------------------------
* 6. Outputs
*-------------------------------------------------------------------------------
frame change res
label define lsex 1 "male" 2 "female"
label values sex lsex
gen surv_rate_per1000 = surv_rate * 1000
gen pop_rate_per1000 = pop_rate * 1000
format surv_rate* pop_rate* %9.2f
format eq_age* gap* %9.1f
gen gapstr = string(gap, "%4.1f") + " (" + string(gap_lo, "%4.1f") + " to " + string(gap_hi, "%4.1f") + ")"
gen eqstr  = string(eq_age, "%4.1f") + " (" + string(eq_age_lo, "%4.1f") + " to " + string(eq_age_hi, "%4.1f") + ")"
sort group sex age
export delimited group sex age deaths py surv_rate_per1000 pop_rate_per1000 eq_age eq_age_lo eq_age_hi gap gap_lo gap_hi ///
    using "$output/fig_mortality_age_gap_data.csv" if group == "all", replace
export delimited group sex age surv_rate_per1000 pop_rate_per1000 eqstr gapstr eq_age eq_age_lo eq_age_hi gap gap_lo gap_hi ///
    using "$output/table_mortality_age_gap.csv" if group == "all" & inlist(age, 20, 30, 40, 50, 60), replace
export delimited group sex age deaths py surv_rate_per1000 pop_rate_per1000 eqstr gapstr eq_age eq_age_lo eq_age_hi gap gap_lo gap_hi ///
    using "$output/table_mortality_age_gap_bydiag.csv" if group != "all", replace
list group sex age surv_rate_per1000 pop_rate_per1000 eqstr gapstr if group == "all" & inlist(age, 20, 30, 40, 50, 60), noobs clean
list group sex age deaths surv_rate_per1000 eqstr gapstr if group != "all", noobs clean

* figure: gap vs age by sex with 95% CI band
twoway (rarea gap_lo gap_hi age if group == "all" & sex == 1, color("31 95 160%25") lwidth(none)) ///
       (rarea gap_lo gap_hi age if group == "all" & sex == 2, color("200 86 43%25") lwidth(none)) ///
       (line gap age if group == "all" & sex == 1, lcolor("31 95 160") lwidth(medthick)) ///
       (line gap age if group == "all" & sex == 2, lcolor("200 86 43") lwidth(medthick)), ///
       yline(0, lcolor(gs10) lpattern(dash)) ///
       xtitle("Attained age of survivor (years)") ytitle("Mortality age gap (years)") ///
       xlabel(20(5)65) ylabel(, angle(horizontal) format(%3.0f)) ///
       legend(order(3 "Male" 4 "Female") ring(0) position(1) cols(1) region(lstyle(none))) ///
       title("Years by which a survivor's all-cause mortality exceeds their age", size(medsmall)) ///
       subtitle("Equivalent general-population age (England & Wales 2022) minus attained age, 95% CI", size(small)) ///
       graphregion(color(white)) plotregion(lstyle(none)) scheme(s1color)
graph export "$output/fig_mortality_age_gap.png", replace width(2100) height(1500)

exit
