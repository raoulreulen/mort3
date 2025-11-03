/*==============================================================================
MERGE WITH BCCSS COHORT FILE AND DEFINE COD CATEGORIES

==============================================================================*/

*-------------------------------------------------------------------------------
* DEFINE EACH COD CATEGORY
*-------------------------------------------------------------------------------

use "$temp/x-mort3-readandlabel" , clear

rename (cat) (codcat)

* relabel the codcat variable
label define CAT 0 "unknown" 1 "infection" 2 "neoplasm" 3 "blood" 			 ///
	4 "endocrine" 5 "mental" 6 "nervous" 7 "circulation" 8 "respiratory"	 ///
	9 "digestive" 10 "muscoskeletal" 11 "genitourinary" 12 "perinatal" 13 "other" 14 	///
	"external" , modify
	
* Create the dummy variables with a temporary prefix
tabulate codcat, generate(codcat_)

* Rename them based on the lowercase value labels
levelsof codcat, local(levels)

* create dummies for each cod category
local i = 1
foreach val of local levels {
    local label : label CAT `val'
    rename (codcat_`i') (`label')
	local ++i
	}	
	
//any death even if cause unknown	
gen allcauses = 1 	

//definite and possible suicide 
recode suicide ( 1 2 = 1) 

//COD =SPN/recurrence
gen spn 	= cond(codcat==2 & neop_cod==2,1,0)
gen recur 	= cond(codcat==2 & neop_cod==1,1,0)

//allcauses except recucurrence
gen allcodexcrecur = cond(allcauses==1 & recur!=1 , 1, 0)

tempfile codcat
save `codcat'




*-------------------------------------------------------------------------------
* MERGE BASIC COHORT FILE BCCSS WITH COD FILE 
*-------------------------------------------------------------------------------

//read spn cohort file (need new file from Dave)
use "$data/bccssbasiccohortfromspn"  , clear 

// merge with entire cohort from SPN (thif file doesnt have dob for example)
merge 1:1 indexno using `codcat'
drop _merge


rename (exit2020 survdate) (dox doe)


//add one day if doe and dox are same
replace dox = dox + (1/365.25) if doe==dox

//calculate age at exit 
gen agex = (dox - dob)/365

save "$temp/x-mort3-prepforstset" , replace

exit
