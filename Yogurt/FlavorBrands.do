clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 20000
set more off

cd "/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt"

//-------------------------------------
// Read flavor and brand codes
//-------------------------------------

import delimited "flavors.tsv", clear
rename v1 fl
rename v2 flavor
save "flavors.dta", replace

import delimited "brands.tsv", clear
rename v1 br
rename v2 brand
save "brands.dta", replace
