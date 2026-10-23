/*==============================================================================
 STSET DATA FOR EACH DEATH OUTCOME SEPARATELY - RATES EXTENDED TO 2024
 Copy of mort3-02a-rrhb-stsetdata-codrefinements.do with ONE change:
   population rates are now used up to 2024 (all rate files were updated to
   1950-2024 on 26 Aug 2026), instead of being capped at 2019. Follow-up ends
   31 Dec 2022, so this means person-time in 2020-2022 is compared with the
   contemporaneous (COVID-era) general population rates.
 Everything else (cod list, recur rate = 0, allcodexcrecur reusing allcauses)
 is unchanged.  RCR 23 Oct 2026
==============================================================================*/
use "$temp/x-mort3-prepforstset" , clear
tempfile mort3 
save `mort3'

local ratecap 2024

local cod allcauses allcodexcrecur spn recur infection blood endocrine mental nervous ///
circulation cardiac cerebrovascular othercirculatory ///
respiratory digestive muscoskeletal genitourinary perinatal ///
external suicide other

foreach x in `cod' {
use `mort3' , clear
	stset dox, fail(`x'==1) id(indexno) origin(dob) entry(doe) scale(365.25)
	assert _st==1
	stsplit ageband, at(0 1 5(5)85 110) 
	stsplit yeargrp, after(time=d(1/1/1900)) at(45(1)122) 
	replace yeargrp = 1900 + yeargrp
	replace yeargrp=1950 if yeargrp <1950
	replace yeargrp=`ratecap' if yeargrp >`ratecap'
	drop if round(_t0, 0.001)==4.999 & _t==5
	sort sex ageband yeargrp
	if "`x'"=="recur" {
		gen newrate = 0
	}
	else if "`x'"=="allcodexcrecur" {
		merge m:1 sex ageband yeargrp using "$rates/allcauses" , keepusing(newrate)
		assert _merge!=1
		keep if _merge==3
	}
	else {
		merge m:1 sex ageband yeargrp using "$rates/`x'" , keepusing(newrate)
		assert _merge!=1
		keep if _merge==3
	}
	strate , smr(newrate)
	gen rate_`x' = newrate
	save  "$temp/x-mort3-stset-`x'"  , replace
}
exit 
