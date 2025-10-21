/*==============================================================================
TABLE 1 - COHORT CHARACTERISTICS


=============================================================================*/

version 16
cap: erase "$temp/tab1_pcsf.xls"
pause


local cod allcauses 

use  "$temp/x-mort3-prepforstset" , clear
stset dox, fail(`cod'==1) id(indexno) origin(dob) entry(doe) scale(365.25)


//attained age 
stsplit attage , at(0 20(10)60 110) after(time=dob)
label var attage "attained age (categorical)"

//followup time
stsplit fu , at(0 10(10)40 )  after(time=doe) trim  //Survdate is constant, doe is not
label var fu "time since diagnosis (categorical)"

//WRONG 
gen fux = (dox-fptdate)/365.25
label var fux "time since diagnosis (continuous)"

//age at each exit (take last one!)
gen agex = (dox-dob)/365.25
label var agex "Attained age (continuous)"

stvary sex diag doe //check these variables do not vary
assert  r(varies) ==  0

tempfile dataset
save `dataset' , replace



*-------------------------------------------------------------------------------
* VARIABLES REPORTED IN TABLE
*-------------------------------------------------------------------------------
#delimit ;
local table1opts sex cat\diag cat\
agedxcat cat\ agedx contn %2.1f\ agedx conts %2.1f \
decdxcat cat\ decdx contn %2.1f \ decdx conts %2.1f \
attage cat\agex contn %2.1f\ agex conts %2.1f \
fu cat\ fux contn %2.1f\ fux conts %2.1f \
rt cat \ rtcran cat ;
#delimit cr

*===============================================================================
* TABLE1 APPROACH
*===============================================================================
use `dataset' , clear
bysort indexno: keep if _n==_N //need to keep last record
sum indexno
assert `r(N)' ==  34488

*keep only above age 50
*keep if attage >= 50 

*keep if decdxcat==1

//COHORT CHARACTERISTICS (TABLE1_MC for continuous outcomes )
table1, format(%2.1f) missing vars(`table1opts') saving("$temp/tab1_pcsf.xls" , sheet(overall) replace)

//EVENTS - BY SURVIVOR
use `dataset' if _d==1 , clear
*bysort pcsfid:  keep if _N>1
bysort indexno: keep if _n==_N //if more than one digestive per individual
table1  ,  vars(`table1opts') format(%2.1f)  miss saving("$temp/tab1_pcsf.xls", sheet(survivors)) 

//EVENTS ALL
use `dataset' if _d==1 , clear
*bysort pcsfid:  keep if _N>1
*bysort pcsfid: keep if _n==_N //if more than one digestive per individual
table1  ,  vars(`table1opts') format(%2.1f)  miss saving("$temp/tab1_pcsf.xls", sheet(events)) 

9888



*-------------------------------------------------------------------------------
* NEW TABLE COMMAND
*-------------------------------------------------------------------------------
version 17
gen overall =1


table (var) (overall), statistic(fvfrequency overall  sex )  statistic(fvpercent sex)  statistic(mean agedx)  statistic(sd agedx)  nototal



	
*-------------------------------------------------------------------------------
* keep only 
*-------------------------------------------------------------------------------	
ds case , not

#delimit ;
table (var) (case) , nototal style(table-1) 
statistic(fvfreq  	`r(varlist)') 	
statistic(fvpercent  `r(varlist)')
;
#delimit cr

*collect dims
*collect label list result, all
collect recode result fvfrequency = column1 fvpercent = column2
*collect label list result, all
collect layout (var) (case#result[column1 column2])
collect style cell var#result[column2],  nformat(%6.1f) sformat("(%s%%)")
collect preview

//change label header column
collect label dim case "status", modify
collect label levels case 2 "No. of controls (%)" 1 "No. of cases (%)"
collect style header result, level(hide)
collect preview

collect style row split , dups(first)
*collect style row stack, nobinder spacer //adds space between different levels created with eac statistics() option
collect style cell border_block, border(right, pattern(nil))
collect preview

*-------------------------------------------------------------------------------
* TO WORD; PUTDOCX
*-------------------------------------------------------------------------------
cap: putdocx clear
putdocx begin,  font(arial narrow, "9")
collect style putdocx, layout(autofitcontents) ///
title("eTable x. Frequency of individual chemotherapy drugs given by case and control status. Table includes drugs for which at least five individuals were exposed")
putdocx collect
putdocx save "$temp/x-etablex-chemotherapyfreq.docx", replace



88
*-------------------------------------------------------------------------------
* Additional cohort characteristics Yuehan Wang (Male Breast cancer)
*-------------------------------------------------------------------------------
/*
table diag, c(mean agex median agex p25 agex p75 agex )
table diag, c(mean fux median fux p25 fux p75 fux )
table diag fpnchemo ,miss
table diag fpnrt ,miss

gen rtonly   = 1 if fpnrt==1 & fpnchemo==2
gen ctonly   = 1 if fpnrt==0 & fpnchemo==2
gen rtctonly = 1 if fpnrt==1 & fpnchemo==1
tab diag if ctonly==1
tab diag if rtonly==1
tab diag if rtctonly==1

*/




use `dataset' , clear
bysort pcsfid: gen x =1  if _n==1
collapse (sum) x _d , by(sex diag)


888
*===============================================================================
* FOOTNOTE TABLE 1 - OTHER CATEGORY
*===============================================================================
use  `dataset' , clear
bysort pcsfid: keep if _n==_N
keep if diag==1
gen x =1 
recode mediccc (. 104 105  =122) 

gcollapse (sum) x  , by(medicc)

*===============================================================================
* ADDITIONAL FOR YUEHAN -> MALE BREAST CANCER
*===============================================================================
use  `dataset' , clear

sum agedx , detail
sum agedx if _d==1 , detail





exit



















exit
