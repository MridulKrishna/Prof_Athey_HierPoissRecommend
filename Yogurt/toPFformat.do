clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 20000
set more off

cd "/Users/casta2k/Dropbox/2015 Summer/Athey/Matlab_AtheyCastillo/Yogurt"

// Find number of brands, flavors, brand/flavors, and households
use "purchases.dta", clear
keep hhid
duplicates drop

use "purchases.dta", clear
keep brand
duplicates drop

use "purchases.dta", clear
gen flavor = 100*flav1+flav2
keep flavor
duplicates drop

use "purchases.dta", clear
gen brfl = 10000*brand+100*flav1+flav2
keep flavor
duplicates drop


//------------------------------------------------------------
// By brand
//------------------------------------------------------------
use "purchases.dta", clear

local name brandPurchases
local nameQ brandQuintiles
local catVar brand

keep hhid `catVar'

bys hhid `catVar': gen purchases = _n
bys hhid `catVar': egen rating = max(purchases)

keep hhid `catVar' rating

duplicates drop

xtile quintile = rating, n(5)

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
	drop quintile
	keep if num < `dbsize'*(0.7)
	drop num
	export delimited "`name'/train.tsv", delimiter(tab) novar replace

restore
preserve

	drop quintile
	keep if (num >= `dbsize'*(0.7) & num < `dbsize'*(0.9))
	drop num
	export delimited "`name'/test.tsv", delimiter(tab) novar replace
	
	keep hhid
	export delimited "`name'/test_users.tsv", delimiter(tab) novar replace	

restore
preserve

	drop quintile
	keep if num >= `dbsize'*(0.9)
	drop num
	export delimited "`name'/validation.tsv", delimiter(tab) novar replace
restore

// Export files with quintiles

preserve
	drop rating
	keep if num < `dbsize'*(0.7)
	drop num
	export delimited "`nameQ'/train.tsv", delimiter(tab) novar replace

restore
preserve

	drop rating
	keep if (num >= `dbsize'*(0.7) & num < `dbsize'*(0.9))
	drop num
	export delimited "`nameQ'/test.tsv", delimiter(tab) novar replace
	
	keep hhid
	export delimited "`nameQ'/test_users.tsv", delimiter(tab) novar replace	

restore
preserve

	drop rating
	keep if num >= `dbsize'*(0.9)
	drop num
	export delimited "`nameQ'/validation.tsv", delimiter(tab) novar replace
restore



//------------------------------------------------------------
// By flavor
//------------------------------------------------------------
use "purchases.dta", clear

gen flavor = 100*flav1+flav2

local name flavorPurchases
local nameQ flavorQuintiles
local catVar flavor

keep hhid `catVar'

bys hhid `catVar': gen purchases = _n
bys hhid `catVar': egen rating = max(purchases)

keep hhid `catVar' rating

duplicates drop

xtile quintile = rating, n(5)

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
	drop quintile
	keep if num < `dbsize'*(0.7)
	drop num
	export delimited "`name'/train.tsv", delimiter(tab) novar replace

restore
preserve

	drop quintile
	keep if (num >= `dbsize'*(0.7) & num < `dbsize'*(0.9))
	drop num
	export delimited "`name'/test.tsv", delimiter(tab) novar replace
	
	keep hhid
	export delimited "`name'/test_users.tsv", delimiter(tab) novar replace	

restore
preserve

	drop quintile
	keep if num >= `dbsize'*(0.9)
	drop num
	export delimited "`name'/validation.tsv", delimiter(tab) novar replace
restore

// Export files with quintiles

preserve
	drop rating
	keep if num < `dbsize'*(0.7)
	drop num
	export delimited "`nameQ'/train.tsv", delimiter(tab) novar replace

restore
preserve

	drop rating
	keep if (num >= `dbsize'*(0.7) & num < `dbsize'*(0.9))
	drop num
	export delimited "`nameQ'/test.tsv", delimiter(tab) novar replace
	
	keep hhid
	export delimited "`nameQ'/test_users.tsv", delimiter(tab) novar replace	

restore
preserve

	drop rating
	keep if num >= `dbsize'*(0.9)
	drop num
	export delimited "`nameQ'/validation.tsv", delimiter(tab) novar replace
restore

//------------------------------------------------------------
// By brand/flavor
//------------------------------------------------------------
use "purchases.dta", clear

gen brfl = 10000*brand+100*flav1+flav2

local name brflPurchases
local nameQ brflQuintiles
local catVar brfl

keep hhid `catVar'

bys hhid `catVar': gen purchases = _n
bys hhid `catVar': egen rating = max(purchases)

keep hhid `catVar' rating

duplicates drop

xtile quintile = rating, n(5)

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
	drop quintile
	keep if num < `dbsize'*(0.7)
	drop num
	export delimited "`name'/train.tsv", delimiter(tab) novar replace

restore
preserve

	drop quintile
	keep if (num >= `dbsize'*(0.7) & num < `dbsize'*(0.9))
	drop num
	export delimited "`name'/test.tsv", delimiter(tab) novar replace
	
	keep hhid
	export delimited "`name'/test_users.tsv", delimiter(tab) novar replace	

restore
preserve

	drop quintile
	keep if num >= `dbsize'*(0.9)
	drop num
	export delimited "`name'/validation.tsv", delimiter(tab) novar replace
restore

// Export files with quintiles

preserve
	drop rating
	keep if num < `dbsize'*(0.7)
	drop num
	export delimited "`nameQ'/train.tsv", delimiter(tab) novar replace

restore
preserve

	drop rating
	keep if (num >= `dbsize'*(0.7) & num < `dbsize'*(0.9))
	drop num
	export delimited "`nameQ'/test.tsv", delimiter(tab) novar replace
	
	keep hhid
	export delimited "`nameQ'/test_users.tsv", delimiter(tab) novar replace	

restore
preserve

	drop rating
	keep if num >= `dbsize'*(0.9)
	drop num
	export delimited "`nameQ'/validation.tsv", delimiter(tab) novar replace
restore
