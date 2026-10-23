/*=========================================================================================
RUN THIS DO-FILE FIRST WHEN STARTING STATA ON RAOUL'S MACBOOK - SETS GLOBALS
(profile.do already sets $temp and $gitpath; this adds the mort3-specific paths)
========================================================================================*/

global od		"/Users/raoulreulen/Library/CloudStorage/OneDrive-UniversityofBirmingham/mort3"
global data 	"$od/data"
global rates 	"$od/rates"
global output 	"$od/output"
global scripts 	"/Users/raoulreulen/Documents/GitHub/mort3"
global temp 	"/Users/raoulreulen/Documents/TEMPWORK"

adopath ++ "$scripts"   // smraer.ado lives in the GitHub folder
cd "$scripts"
exit
