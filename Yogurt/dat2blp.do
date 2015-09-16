clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 20000
set more off

cd "/Users/casta2k/Dropbox/2015 Summer/Athey/Matlab_AtheyCastillo/Yogurt"

use "purchases.dta", clear

gen product = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f")+string(lf)+string(nf)+string(weightamount/100,"%03.0f")
gen weight = weightamount*multipack

//keep if week >=198800

// Get number of weeks and products
tab week, nofreq
local nWeeks = r(r)
disp `nWeeks'
tab product, nofreq
local nProducts = r(r)
disp `nProducts'

bys product market week: egen weights = sum(weight)

local varsProd product market week brand weights flav1 flav2 lf nf
keep `varsProd'
duplicates drop

gen marketWeek = string(market)+string(week)
drop market week

reshape wide weights, j(marketWeek) i(product) string

foreach wvar of varlist weights* {
	replace `wvar' = 2000 if `wvar'==.
}

reshape long weights, i(product) j(marketWeek) string

gen market = substr(marketWeek,1,1)
gen week = substr(marketWeek,2,.)

bys market week: egen sales = sum(weights)
gen share = weights/sales

disp `nWeeks'*`nProducts'
