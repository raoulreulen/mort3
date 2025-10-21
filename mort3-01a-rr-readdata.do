/*========================================================================================
BCCSS MORTALITY DATA - MENTAL HEALTH PROJECT (HAMA)
READ DATA (N=5,922 DEATHS)
DATE: 12JUne2025
========================================================================================*/

*-----------------------------------------------------------------------------------------
* READ DATA AND APPLY OFFSETS
*-----------------------------------------------------------------------------------------
use "$data/BCCSSMort_Apr25.dta" , clear


//merge with recurrence/spn distinction file
merge 1:1 lngpk_bccss o_indexno using "$data/BCCSSMort_Apr25_NEOP.dta" 
assert _merge==3 | _merge==1 //check whether merge was correct
drop _merge

//read COD categories (new file 04sept2025; see email DLW )
merge 1:1 lngpk_bccss o_indexno using "$data/BCCSSMortCODCat.dta" , ///
keepusing(cat suicide circcat neop_recurr)
assert _merge==3 //check whether merge was correct
drop _merge

//* merge with offset file (using common primary key) and keep correctly merged records
merge 1:1 lngpk_bccss using "$temp/OFFSET_DATA_BCCSS.dta"
assert _merge!=1
keep if _merge == 3
drop _merge

// generate native variables
gen indexno = o_indexno - indexno_offset
gen cohort	= o_cohort 	- cohort_offset
gen country = o_country - country_offset
gen icdver 	= o_icdver 	- icdver_offset

gen dod_d = o_dod_day - dod_offset
gen dod_m = o_dod_mth - dod_offset
gen dod_y = o_dod_yr  - dod_offset

gen dod = mdy(dod_m,dod_d,dod_y)	//create dod variable from 3 variables
format dod %td  					//format into Stata date format
drop dod_d dod_m dod_y

// tidy up
drop o_* *offset


*-----------------------------------------------------------------------------------------
* LABEL VARIABLES
*-----------------------------------------------------------------------------------------
label var indexno 	"unique identifier"
label var ucause  	"underlying cause of death"
label var icdver	"ICD version applicable at time death"
label var country 	"country of cancer diagnosis"
label var cohort 	"cohort up to 1991 (1); cohort 1992-2006 (2)"

label define lblyesno 0 "no"  1 "yes"

*label define lblcohort 1 "up to 1991 (1)" 2 "1992-2006 (2)"
*label values cohort lblcohort

label define lblcountry 1 "England (1)" 2 "Scotland (2)" 3 "Wales (3)"  //check with Dave 
label values country lblcountry 


order indexno lngpk_bccss cohort country dod icdver
save "$temp/x-mort3-readandlabel" , replace

exit




