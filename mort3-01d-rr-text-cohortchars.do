/*==============================================================================
COHORT CHARACTERISTICS: PERSON-YEARS, ATTAINED AGE, FOLLOW-UP TIME
=============================================================================*/

use "$temp/x-mort3-prepforstset", clear
stset dox, fail(allcauses==1) id(indexno) origin(dob) entry(doe) scale(365.25)
assert _st==1

*-------------------------------------------------------------------------------
* 1. TOTAL PERSON-YEARS
*-------------------------------------------------------------------------------
quietly gen _py = _t - _t0
quietly sum _py
scalar s_py = r(sum)

*-------------------------------------------------------------------------------
* 2. ATTAINED AGE AT END OF STUDY (before stsplit - one record per person)
*-------------------------------------------------------------------------------
gen attagex = (dox - dob) / 365.25
quietly sum attagex, detail
scalar s_agex_min    = r(min)
scalar s_agex_max    = r(max)
scalar s_agex_median = r(p50)

*-------------------------------------------------------------------------------
* 3. EMBARKATION (flag2020==2)
*-------------------------------------------------------------------------------
quietly count if flag2020 == 2
scalar s_embark = r(N)

*-------------------------------------------------------------------------------
* 4. % REACHING AGE 40+ AND 50+ (stsplit)
*-------------------------------------------------------------------------------
quietly count
scalar s_n = r(N)

stsplit agegrp, at(0 40 50 200) after(time=dob)

quietly {
	bysort indexno: egen byte reached40 = max(agegrp >= 40)
	bysort indexno: egen byte reached50 = max(agegrp >= 50)
	bysort indexno: keep if _n == _N
	count if reached40
	scalar s_n40 = r(N)
	count if reached50
	scalar s_n50 = r(N)
}

*-------------------------------------------------------------------------------
* 5. MEDIAN FOLLOW-UP TIME (time since entry - re-stset with origin=doe)
*-------------------------------------------------------------------------------
use "$temp/x-mort3-prepforstset", clear
stset dox, fail(allcauses==1) id(indexno) origin(doe) entry(doe) scale(365.25)
assert _st==1

quietly sum _t, detail
scalar s_fup_p25    = r(p25)
scalar s_fup_median = r(p50)
scalar s_fup_p75    = r(p75)

*-------------------------------------------------------------------------------
* BUILD TABLE AND OUTPUT TO WORD
*-------------------------------------------------------------------------------
clear
set obs 8
gen str80 Characteristic = ""
gen str30 Value = ""

replace Characteristic = "Total person-years"                              in 1
replace Characteristic = "Attained age at end of study: minimum"           in 2
replace Characteristic = "Attained age at end of study: median"            in 3
replace Characteristic = "Attained age at end of study: maximum"           in 4
replace Characteristic = "Reached age 40 or older, %"                     in 5
replace Characteristic = "Reached age 50 or older, %"                     in 6
replace Characteristic = "Median follow-up time, years (IQR)"             in 7
replace Characteristic = "Embarkation (flag2020=2), n"                    in 8

replace Value = string(s_py, "%12.0fc")                                              in 1
replace Value = string(s_agex_min,    "%4.1f")                                       in 2
replace Value = string(s_agex_median, "%4.1f")                                       in 3
replace Value = string(s_agex_max,    "%4.1f")                                       in 4
replace Value = string((s_n40/s_n)*100, "%4.1f")                                     in 5
replace Value = string((s_n50/s_n)*100, "%4.1f")                                     in 6
replace Value = string(s_fup_median, "%4.1f") + " (" + ///
    string(s_fup_p25, "%4.1f") + "-" + string(s_fup_p75, "%4.1f") + ")"             in 7
replace Value = string(s_embark, "%12.0fc") + " (" + ///
    string((s_embark/s_n)*100, "%4.1f") + "%)"                                      in 8

cap: putdocx clear
putdocx begin, font(arial narrow, "10")
putdocx paragraph

putdocx table tbl = data("Characteristic Value"), varnames          ///
    border(start, nil) border(insideV, nil) border(end, nil) width(100%)    ///
    layout(autofitcontents) border(all, nil)                                ///
    title("Table x. Cohort follow-up characteristics")

putdocx table tbl(1,.), bold
putdocx table tbl(.,1), halign(left)
putdocx table tbl(.,2), halign(right)

putdocx save "$temp/text-cohortchars.docx", replace

exit
