clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 20000
set more off

cd "/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt"

import delimited "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hbeta.tsv", clear
rename v1 itemn
rename v2 item
rename v* beta*
save "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hbeta.dta", replace

import delimited "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hrho.tsv", clear
rename v1 itemn
rename v2 item
rename v* rho*
save "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hrho.dta", replace

import delimited "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/htheta.tsv", clear
rename v1 usern
rename v2 user
rename v* theta*
save "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/htheta.dta", replace

import delimited "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hsigma.tsv", clear
rename v1 usern
rename v2 user
rename v* sigma*
save "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hsigma.dta", replace

import delimited "observables/obsUser.tsv", clear
rename v1 user
rename v* userObs*
save "observables/obsUser.dta", replace

import delimited "observables/obsItem.tsv", clear
rename v1 item
rename v* itemObs*
save "observables/obsItem.dta", replace

use "observables/obsUser.dta", clear

merge 1:1 user using "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/htheta.dta" 
drop _merge
merge 1:1 user using "observables/output/n2599-m275-k25-uc36-ic50-batch-hier-vb/hsigma.dta"
drop _merge

