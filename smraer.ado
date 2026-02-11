*! version 1.1 RCR - SMRAER
cap program drop smraer
program define smraer , sortpreserve 

syntax varlist(min=3 max=3 numeric) [, dessmr(integer 1) desaer(integer 1)]

if "`varlist'"!="" {
		tokenize `varlist' 
		}
		
	tempvar obs exp pyrs
	quietly {
	gen `obs'	=`1'
	gen `exp'	=`2'	
	gen `pyrs' 	=`3'


*-------------------------------------------------------------------------------
*--CALCULATE SMRs(Exact confidence intervals)
*-------------------------------------------------------------------------------
gen smr		= (`obs'/`exp')
gen smrll	= ((invgammap( `obs',(0.05)/2))/`exp') 
gen smrul	= ((invgammap((`obs'+ 1), (1.95)/2))/`exp') 
	
*-------------------------------------------------------------------------------	
*--CALCULATE AERs (Exact confidence intervals)
*-------------------------------------------------------------------------------
gen aer = ((`obs'- `exp')/`pyrs')
*gen aer		= cond((`obs' - `exp')>0 , ((`obs'- `exp')/`pyrs') , 0)
*gen aerll	= aer - (1.96*(sqrt(`obs')/y))
*gen aerul	= aer + (1.96*(sqrt(`obs')/y))
gen aerll = ((invgammap( abs(`obs' - `exp')	  , (0.05)/2))/`pyrs') 
gen aerul = ((invgammap(abs((`obs'-  `exp')+ 1)  , (1.95)/2))/`pyrs') 

*-------------------------------------------------------------------------------	
* PRESENT AS STRING VARIABLE
*-------------------------------------------------------------------------------
gen str smrstr = string(smr , "%9.`dessmr'f") + " (" + string(smrll , "%9.`dessmr'f") ///
+ "," + string(smrul , "%9.`dessmr'f") + ")" /*+ " (" + string(`obs') + ")"*/
*drop smr smrll smrul 
	
gen str aerstr = string(aer , "%9.`desaer'f") + " (" + string(aerll , "%9.`desaer'f") ///
 + "," + string(aerul , "%9.`desaer'f") + ")"  /*  + " (" + string(diagd) + ")" */
 
*-------------------------------------------------------------------------------	
* PDROP STUFF AND TIDY UP
*-------------------------------------------------------------------------------

gen strobs 	= string(`obs',"%9.0f")
gen strexp 	= string(`exp', "%9.1f")  
gen obsexp 	= string(`obs',"%9.0f") + "/" + string(`exp', "%9.1f")  
replace obsexp  =  "" if obsexp=="./." 
replace smrstr  =  "" if smrstr==". (.,.)"
replace aerstr  =  "" if aerstr==". (.,.)"
 
drop `exp' 
drop `pyrs'
}


end
