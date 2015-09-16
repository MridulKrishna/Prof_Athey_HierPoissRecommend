clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 20000
set more off

cd "/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt"

// Check number of products
use "purchases.dta", clear
gen product = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f")+string(lf)+string(nf)+string(weightamount,"%05.0f")
keep product
duplicates drop

// Check number of users
use "purchases.dta", clear
keep hhid
duplicates drop

// Actual code starts here

//-------------------------------------
// Ratings
//-------------------------------------
use "purchases.dta", clear
gen item = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f")+string(lf)+string(nf)+string(weightamount,"%05.0f")

keep hhid item multipack purchased
gen n = multipack*purchased
drop multipack purchased

bys hhid item: egen rating = sum(n)

local name observables
local catVar item

keep hhid `catVar' rating

save "`name'.dta", replace

set seed 10
gen r = runiform()
sort r
gen num = _n
drop r

sort hhid `catVar'

local dbsize = _N

	keep if num < `dbsize'*(0.7)
	drop num

keep hhid
duplicates drop
