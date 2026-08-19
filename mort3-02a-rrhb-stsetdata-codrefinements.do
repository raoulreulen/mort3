/*==============================================================================
 STSET DATA FOR EACH DEATH OUTCOME SEPARATELY 

 MODIFIED BY HAMA 
- included circulatory subcategories: added cardiac, cerebrovascular, othercirculatory 
  to local cod list. Rates capped at 2019 for consistency with all other categories.
- dropped neoplasm from local cod list (it double-counted spn + recur combined); 
  replaced with spn and recur as separate categories
- recur (recurrence/progression of primary cancer) has no valid population comparator, 
  so expected rate is set to 0 (SMR = NA, AER reduces to crude rate)
- added allcodexcrecur (all-cause mortality excluding recurrence); reuses the 
  allcauses rates file since the general population comparator is unchanged
Original stset do-file by Raoul untouched, this is a modified copy 
==============================================================================*/

use "$temp/x-mort3-prepforstset" , clear
tempfile mort3 
save `mort3'

*-------------------------------------------------------------------------------
* DEFINE LOCAL (PUT IN HERE THE COD THAT YOU WANT TO RUN STSET FOR)
*-------------------------------------------------------------------------------
local cod allcauses allcodexcrecur spn recur infection blood endocrine mental nervous ///
circulation cardiac cerebrovascular othercirculatory ///
respiratory digestive muscoskeletal genitourinary perinatal ///
external suicide other

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
	replace yeargrp=2019 if yeargrp >2019 //need updated mortality rates
	drop if round(_t0, 0.001)==4.999 & _t==5
	
	sort sex ageband yeargrp
	
	//special handling for categories without a standard rates file
	if "`x'"=="recur" {
		//recurrence/progression of primary cancer has no valid population 
		//comparator (general population cannot "recur" a cancer they never had)
		//expected rate = 0 -> SMR undefined (NA), AER reduces to crude rate
		gen newrate = 0
	}
	else if "`x'"=="allcodexcrecur" {
		//all-cause mortality excluding recurrence - reuse allcauses rates,
		//since the general population comparator is the same regardless of 
		//which subset of survivor deaths we're comparing
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
	
	//safe file for each cause of death	
	save  "$temp/x-mort3-stset-`x'"  , replace
}
exit 


