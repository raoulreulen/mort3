/*=================================================================
RENAME VARIABLES IN THE NEW RATES FILES TO MATCH EXISTING RATE FILES 
*****RUN ONLY ONCE - no not re-run as part of the regular pipeline - change has already been applied permanently******
Prepares the rate files cardiac.dta, cerebrovascular.dta, othercirculatory.dta for use in the STSET pipeline for the circulatory subcategories. 

======================================================================*/

use "$rates/cardiac.dta", clear 
rename native_rate newrate
save "$rates/cardiac.dta", replace 

use "$rates/cerebrovascular.dta", clear 
rename native_rate newrate 
save "$rates/cerebrovascular.dta", replace 

use "$rates/othercirculatory.dta", clear
rename native_rate newrate
save "$rates/othercirculatory.dta", replace


/*========================================================================================
CHECK YEAR RANGE ACROSS ALL RATES FILES
========================================================================================*/
local ratefiles allcauses spn neoplasm infection blood endocrine mental nervous circulation ///
respiratory digestive muscoskeletal genitourinary perinatal other external suicide ///
cardiac cerebrovascular othercirculatory

foreach x of local ratefiles {
	capture confirm file "$rates/`x'.dta"
	if !_rc {
		use "$rates/`x'.dta", clear
		qui: sum yeargrp
		di as text "`x'" _col(25) "min=" r(min) _col(35) "max=" r(max) _col(45) "N=" r(N)
	}
	else {
		di as error "`x'.dta not found in $rates"
	}
}

/*---------------------------------------------------------------------
FINDINGS OF THE ABOVE CHECK - DONE ON 19/08/26 
All 15 original cause of death rates files (allcauses, spn, neoplasm.... suicide) span 
from 1950-2019, consistent with the exisiting yeargrp cap of 2019 in the STSET do-file 

The 3 new circulatory subcategory rate files (cardiac, cerebrovascular and othercirculatory)
extend further to 2024. 

DECISION: For the initial analysis, all categories including the new subcategories will be 
capped at 2019 in the STSET do-file, to maintain consistency. 
But I think it is better to ask Dave for updated rates uo to 2024 for all 15 original categories
and re-run the analysis with them 
-------------------------------------------------------------------------------------*/ 


exit 
