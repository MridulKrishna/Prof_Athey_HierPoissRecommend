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
replace purchased = 1 if purchased == 0
gen item = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f")+string(lf)+string(nf)+string(weightamount,"%05.0f")

keep hhid item multipack purchased
gen n = multipack*purchased
drop multipack purchased

bys hhid item: egen rating = sum(n)

local name observables
local catVar item

keep hhid `catVar' rating
duplicates drop

save "`name'.dta", replace

set seed 10
gen r = runiform()
sort r
gen num = _n
drop r

sort hhid `catVar'

local dbsize = _N

// Export files with ratings
preserve
	keep if num < `dbsize'*(0.7)
	drop num
	export delimited "`name'/train.tsv", delimiter(tab) novar replace

restore
preserve
	keep if (num >= `dbsize'*(0.7) & num < `dbsize'*(0.9))
	drop num
	export delimited "`name'/test.tsv", delimiter(tab) novar replace
	
	keep hhid
	export delimited "`name'/test_users.tsv", delimiter(tab) novar replace	
restore
preserve
	keep if num >= `dbsize'*(0.9)
	drop num
	export delimited "`name'/validation.tsv", delimiter(tab) novar replace
restore

//-------------------------------------
// Item observables
//-------------------------------------
use "purchases.dta", clear
gen item = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f")+string(lf)+string(nf)+string(weightamount,"%05.0f")

keep item brand weightamount flav1 flav2 lf nf 
duplicates drop

// Check that item uniquely identifies data
duplicates report item

// Generates dummies for brand and flavor
rename brand br

forvalues i = 1/22 {
	gen brand`i' = (br==`i')
}
drop br

forvalues i = 1/34 {
	gen flavor`i' = (flav1==`i' | flav2==`i')
}
drop flav1 flav2

forvalues i = 1/22 {
	quietly sum brand`i'
	display "`i'="r(sum)
	
	if r(sum) == 0 {
		drop brand`i'
	}
}

forvalues i = 1/34 {
	quietly sum flavor`i'
	display "`i'="r(sum)
	
	if r(sum) == 0 {
		drop flavor`i'
	}
}

order item

sort item

export delimited "`name'/obsItem.tsv", delimiter(tab) novar replace

//-------------------------------------
// User observables
//-------------------------------------
use "purchases.dta", clear

keep hhid hhsize hhincome restype mrace mheadwstat frace fheadwstat smom
duplicates drop

// Check that hhid uniquely identifies data
duplicates report hhid

forvalues i = 1/6 {
	gen restype`i' = (restype==`i')
	quietly sum restype`i'
	display "`i'="r(sum)
}
drop restype

forvalues i = 0/5 {
	gen mrace`i' = (mrace==`i')
	quietly sum mrace`i'
	display "`i'="r(sum)
}
drop mrace

forvalues i = 0/7 {
	gen mheadwstat`i' = (mheadwstat==`i')
	quietly sum mheadwstat`i'
	display "`i'="r(sum)
	
	if r(sum) == 0 {
		drop mheadwstat`i'
	}
}
drop mheadwstat

forvalues i = 0/5 {
	gen frace`i' = (frace==`i')
	quietly sum frace`i'
	display "`i'="r(sum)
}
drop frace

forvalues i = 0/7 {
	gen fheadwstat`i' = (fheadwstat==`i')
	quietly sum fheadwstat`i'
	display "`i'="r(sum)
}
drop fheadwstat

sort hhid

export delimited "`name'/obsUser.tsv", delimiter(tab) novar replace

//-------------------------------------
// Check most common products
//-------------------------------------
use "purchases.dta", clear
replace purchased = 1 if purchased == 0
gen item = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f")+string(lf)+string(nf)+string(weightamount,"%05.0f")

keep item multipack purchased
gen n = multipack*purchased
drop multipack purchased

bys item: egen rating = sum(n)

drop n
duplicates drop

gen brst = substr(item,1,2)
gen fl1st = substr(item,3,2)
gen fl2st = substr(item,5,2)
gen lfst = substr(item,7,1)
gen nfst = substr(item,8,1)
gen wst = substr(item,9,5)
destring brst, gen(br)
destring fl1st, gen(fl1)
destring fl2st, gen(fl2)
destring lfst, gen(lf)
destring nfst, gen(nf)
destring wst, gen(weight)
drop brst fl1st fl2st lfst nfst wst

merge n:1 br using "brands.dta"
drop _merge br

rename fl1 fl
merge n:1 fl using "flavors.dta"
drop _merge

drop fl
rename flavor flavor1

rename fl2 fl
merge n:1 fl using "flavors.dta"
drop _merge

drop fl
rename flavor flavor2

drop if item == ""

gsort -rating

order brand flavor1 flavor2 lf nf weight rating

//-------------------------------------
// Check greatest buyers
//-------------------------------------
use "purchases.dta", clear
replace purchased = 1 if purchased == 0

keep hhid multipack purchased
gen n = multipack*purchased
drop multipack purchased

bys hhid: egen rating = sum(n)

drop n
duplicates drop

gsort -rating

order brand flavor1 flavor2 lf nf weight rating

//-------------------------------------
// Check greatest buyer
//-------------------------------------
use "purchases.dta", clear
replace purchased = 1 if purchased == 0
1151894
