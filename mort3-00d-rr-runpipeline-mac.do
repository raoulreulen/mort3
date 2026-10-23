/*==============================================================================
 MASTER: REBUILD THE ANALYSIS DATASETS ON THE MAC (prep + stset)
==============================================================================*/
do "/Users/raoulreulen/Documents/GitHub/mort3/mort3-00c-rr-setwd-mac.do"
do "$scripts/mort3-01a-rrhb-readdata.do"
do "$scripts/mort3-01b-rr-covardata.do"
do "$scripts/mort3-01b-rrhb-gencodcats.do"
do "$scripts/mort3-02a-rr-stsetdata-rates2024.do"
exit
