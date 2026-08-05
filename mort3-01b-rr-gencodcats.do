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
	9 "digestive" 10 "muscoskeletal" 11 "genitourinary" 12 "perinatal" 13 "other" 14 ///
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
use "$temp/x-mort3-covarbccssdata"  , clear 

// merge with entire cohort from SPN (thif file doesnt have dob for example)
merge 1:1 indexno using `codcat'
drop _merge


//TYEP OF CHILDHOOD CANCER (MEDICCC)
gen diag = .
replace diag = 1 	if inrange(mediccc, 11, 15) & mediccc!=12 //12=AML
replace diag = 2 	if mediccc==12
replace diag = 3 	if mediccc==21
replace diag = 4 	if mediccc==22
replace diag = 5 	if inrange(mediccc , 31, 36)
replace diag = 6 	if mediccc==41
replace diag = 7 	if mediccc==51 & genretino==1
replace diag = 8 	if mediccc==51 & genretino==0
replace diag = 9 	if mediccc==61
replace diag = 10 	if inrange(mediccc, 81, 85)
replace diag = 11 	if inrange(mediccc, 91, 95)
replace diag = 12 	if inlist(mediccc,23,24, 25, 42, 62, 63, 71, 72, 73) | 	///
					   inrange(mediccc, 101, 122)  | medicc==292
					   
//NB: check with Dave what is code mediccc 292?					   
					   
label define ldiag 0 "overall" 1 "Leukaemia (except AML)"  2  "AML" 3 "Hodgkins"	///
 4 "NHL" 5 "CNS"  6 "neuroblastoma" 7 "hretino" 8 "nhretino" 						///
9 "Wilms" 10 "bone" 11 "softtissue" 12 "other"

label values diag ldiag

assert !mi(diag)

//DATE OF ENTRY
gen doe = mdy(month(fpt), day(fpt) , year(fpt) + 5)
replace doe = mdy(3, 1, year(fpt) + 5) if mi(doe) & month(fpt)==2 & day(fpt)==29
format %td doe
label var doe "date of entry cohort (aka 5-year survival)"

//FPT DATE
rename fpt fptdate

gen fptyear = year(fptdate)
egen decdxcat =  cut(fptyear) , at(1900 1970 1980 1990 2000 2020)

//AGE AT DIAGNOSIS
gen agedx = (fptdate-dob)/365.24
egen agedxcat = cut(agedx) , at(0 4 8 12 20)


//add one day if doe and dox are same
replace dox = dox + (1/365.25) if doe==dox

//calculate age at exit 
gen agex = (dox - dob)/365

save "$temp/x-mort3-prepforstset" , replace

exit
