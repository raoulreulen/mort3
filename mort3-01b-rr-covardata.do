*-----------------------------------------------------------------------------------------
* STEP 1: Cohort Charactaristic-BCCSS
*-----------------------------------------------------------------------------------------

use "$data/BCCSS_CoVariates_Jan26" , clear

/* Merge with revised offset file with additional columns 
(using common primary key) and keep correctly merged records */

merge 1:1 lngpk_bccss using "$temp/OFFSET_DATA_BCCSS_MORT.dta"

//kepp only matched values in both dataset, however, here all (34,488) matched
keep if _merge == 3  


// generate native variables
gen indexno = o_indexno - indexno_offset
gen country = o_country - country_offset
gen sex = o_sex - sex_offset

gen dob_d = o_dob_day - dob_offset
gen dob_m = o_dob_mth - dob_offset
gen dob_y = o_dob_yr - dob_offset

gen dob = mdy(dob_m,dob_d,dob_y)
format dob %td
drop dob_d dob_m dob_y

gen fpt_d = o_fpt_day - fpt_offset
gen fpt_m = o_fpt_mth - fpt_offset
gen fpt_y = o_fpt_yr - fpt_offset

gen fpt = mdy(fpt_m,fpt_d,fpt_y)
format fpt %td
drop fpt_d fpt_m fpt_y


gen agefpn = o_fpt_age - fpt_offset

gen rt = o_rt - trt_offset
gen rtgroupvar = .
gen rtcranioabdo = .

replace rtgroupvar = o_rtgroupvar - trt_offset if o_rtgroupvar != .
replace rtcranioabdo = o_rtcranioabdo - trt_offset if o_rtcranioabdo != .

gen ct = o_ct - trt_offset

gen vitstat = o_vitstat - exit_offset

gen exitd = o_exit_day - exit_offset
gen exitm = o_exit_mth - exit_offset
gen exity = o_exit_yr - exit_offset
gen dox = mdy(exitm,exitd,exity)    //make varible exit date 
format dox %td
drop exitd exitm exity


gen mediccc = o_mediccc - diag_offset
gen diaggrp = o_diaggrp - diag_offset
gen genretino = o_genretino - diag_offset

//drop unnecessary variables
drop o_*
drop _merge
drop *offset


save "$temp/x-mort3-covarbccssdata" , replace
