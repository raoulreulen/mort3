/*==============================================================================
 STSET DATA FOR EACH DEATH OUTCOME SEPARATELY 

==============================================================================*/

use "$temp/x-mort3-prepforstset" , clear
tempfile mort3 
save `mort3'

*-------------------------------------------------------------------------------
* DEFINE LOCAL (PUT IN HERE THE COD THAT YOU WANT TO RUN STSET FOR)
*-------------------------------------------------------------------------------
local cod allcauses spn neoplasm infection blood endocrine mental nervous circulation ///
respiratory digestive muscoskeletal genitourinary perinatal other external suicide

*-------------------------------------------------------------------------------
* LOOP THROUGH EACH COD AND STSET FOR EACH COD CREATING A NEW DATASET
*-------------------------------------------------------------------------------
foreach x in `cod' {
use `mort3' , clear //read the original file
	
	//stset on each cause of death
	stset dox, fail(`x'==1) id(indexno) origin(dob) entry(doe) scale(365.25)
	assert _st==1

	//stsplit
	stsplit ageband, at(0 1 5(5)85 110) 
	stsplit yeargrp, after(time=d(1/1/1900)) at(45(1)122) 
	replace yeargrp = 1900 + yeargrp
	
	//replace years no rates available
	replace yeargrp=1950 if yeargrp <1950
	replace yeargrp=2019 if yeargrp >2019. //need updated mortality rates

	drop if round(_t0, 0.001)==4.999 & _t==5
	
	sort sex ageband yeargrp
	merge m:1 sex ageband yeargrp using "$rates/`x'.dta" , keepusing(newrate)
	assert _merge!=1
	keep if _merge==3
	

	strate , smr(newrate)
	gen rate_`x' = newrate
	
	//safe file for each cause of death	
	save  "$temp/x-mort3-stset-`x'"  , replace
}


exit 
