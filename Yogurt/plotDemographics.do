clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 20000
set more off

cd "/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt"

import delimited "observables/obsUser.tsv"

rename v1 hhid
rename v2 hhsize
rename v3 hhincome
rename v4 smom

forvalues i = 1/6 {
	local pos = `i'+4
	rename v`pos' restype`i'
}

forvalues i = 0/5 {
	local pos = `i'+11
	rename v`pos' mrace`i'
}

forvalues i = 0/5 {
	local pos = `i'+17
	rename v`pos' mheadwstat`i'
}

rename v23 mheadwstat7

forvalues i = 0/5 {
	local pos = `i'+24
	rename v`pos' frace`i'
}

forvalues i = 0/7 {
	local pos = `i'+30
	rename v`pos' fheadwstat`i'
}


save "observables/obsUser.dta", replace

import delimited "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hsigma.tsv", clear

rename v1 number
rename v2 hhid
rename v3 cweight
rename v4 clf
rename v5 cnf

forvalues i = 1/2 {
	local pos = `i'+5
	rename v`pos' cbrand`i'
}

forvalues i = 4/8 {
	local pos = `i'+4
	rename v`pos' cbrand`i'
}

forvalues i = 10/14 {
	local pos = `i'+3
	rename v`pos' cbrand`i'
}

rename v18 cbrand16

forvalues i = 18/22 {
	local pos = `i'+1
	rename v`pos' cbrand`i'
}


forvalues i = 1/6 {
	local pos = `i'+23
	rename v`pos' cflavor`i'
}

forvalues i = 8/18 {
	local pos = `i'+22
	rename v`pos' cflavor`i'
}

forvalues i = 20/22 {
	local pos = `i'+21
	rename v`pos' cflavor`i'
}

forvalues i = 24/30 {
	local pos = `i'+20
	rename v`pos' cflavor`i'
}

forvalues i = 32/32 {
	local pos = `i'+19
	rename v`pos' cflavor`i'
}

forvalues i = 34/34 {
	local pos = `i'+18
	rename v`pos' cflavor`i'
}


merge 1:1 hhid using "observables/obsUser.dta"

twoway (scatter cbrand1 hhincome) (lfit cbrand1 hhincome), ///
	graphregion(color(white)) scheme(s2color) /*legend(off) scale(1.1)*/ ///
	/*xtitle("adsf") ytitle("wefad")*/
graph export "observables/output/income_cweight.pdf", as(pdf) preview(off) replace
