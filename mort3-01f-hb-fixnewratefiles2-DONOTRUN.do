/*=================================================================
RENAME A VARIABLE IN THE NEW RATES FILES (UP TO 2024) TO MATCH EXISTING CODE 
*****RUN ONLY ONCE - no not re-run as part of the regular pipeline - change has already been applied permanently******
=================================================================*/


local cod allcauses spn infection blood endocrine mental nervous circulation respiratory digestive muscoskeletal genitourinary perinatal external suicide other

foreach x of local cod {
    capture confirm file "$rates/`x'.dta"
    if !_rc {
        use "$rates/`x'.dta", clear
        capture confirm variable native_rate
        if !_rc {
            rename native_rate newrate
            save "$rates/`x'.dta", replace
            di as result "`x'.dta: renamed native_rate to newrate"
        }
        else {
            di as text "`x'.dta: native_rate not found - skipping"
        }
    }
    else {
        di as error "`x'.dta not found in $rates"
    }
}

/*=================================================================
Checking that all rate files now extend to 2024
Check done 26/08/26 and they all show this rage: (min=1950  max=2024)
=================================================================*/

local cod allcauses spn infection blood endocrine mental nervous circulation cardiac cerebrovascular othercirculatory respiratory digestive muscoskeletal genitourinary perinatal external suicide other

foreach x of local cod {
    use "$rates/`x'.dta", clear
    qui: sum yeargrp
    di as text "`x'" _col(25) "min=" r(min) _col(35) "max=" r(max)
}