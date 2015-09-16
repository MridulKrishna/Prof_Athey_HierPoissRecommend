clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 20000
set more off

cd "/Users/casta2k/Dropbox/2015 Summer/Athey/Matlab_AtheyCastillo/Yogurt"

use "purchases.dta", clear

bysort hhid week dayofweek: gen purInDay = _n
bysort hhid week: gen purInWeek = _n

replace hhid = round(hhid)

//duplicates report hhid week dayofweek

// Create index for every transaction
gen purchase = string(hhid)+"0"+string(week)+"0"+string(dayofweek)+"0"+string(purInDay)
gen hhweek = string(hhid)+"0"+string(week)
gen hhday = string(hhid)+"0"+string(week)+"0"+string(dayofweek)
order purchase hhid week dayofweek purInDay purInWeek

//duplicates report hhday

save "purchasesIndices.dta", replace


// General statistics for purchases
use "purchasesIndices.dta", clear

tabulate brand
histogram brand, graphregion(color(white)) ///
		scheme(s2color) xtitle("Purchases by brand") fraction
graph export "brands.pdf", as(pdf) replace

// Frequency of outside option
use "purchasesIndices.dta", clear

//gen weekday = string(week)+"0"+string(dayofweek)
keep hhid week brand

duplicates drop hhid week, force

set more off
reshape wide brand, i(hhid) j(week) 
reshape long brand, i(hhid) j(week) 

gen bought = (brand != .)

tabulate bought
histogram bought, graphregion(color(white)) ///
		scheme(s2color) xtitle("Weeks with purchases") fraction
graph export "FreqPurchWeek.pdf", as(pdf) replace

// Statistics for each day
use "purchasesIndices.dta", clear
	
	gen fl1 = flav1
	gen fl2 = flav2
	replace fl1 = 0 if flav1 == .
	replace fl2 = 0 if flav2 == .
	gen brfl = 10000*brand+100*fl1+fl2
	drop flav1 flav2
	
	local vdrop purchase weightamount multipack nf displ feat xprice ustcoupons valcoup ///
		umancoupons valmcoupons purchased uweight purInWeek

	local towide market stid n brand fl1 fl2 brfl lf 
		 
	local hhvars hhsize hhincome restype mrace mheadwstat frace fheadwstat smom

	disp "`vdrop'"
	drop `vdrop'
		
	reshape wide `towide', i(hhday) j(purInDay)

	gen nPur = 0

	forvalues i = 1/16 {
		gen m`i' = (market`i' != .)
		//replace nPur = if(`i'
		
		replace nPur = nPur + m`i'
		drop m`i'
	}
	
	tabulate nPur
	histogram nPur, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of purchases per day") fraction
	graph export "nPurchasesDay.pdf", as(pdf) replace
	
	preserve
		keep hhday brand*
		
		reshape long brand, i(hhday) j(purInDay) 
		bys hhday brand: gen first = _n == 1
		bys hhday: gen numbrands = sum(first)
		drop first
		reshape wide brand numbrands, i(hhday) j(purInDay)
		gen brands = numbrands16
		gen full = (brand16 != . )
		replace brands = brands - 1 +full
		drop numbrands*
		
		tabulate brands
		histogram brands, graphregion(color(white)) ///
				scheme(s2color) xtitle("Number of brands purchased in a day") fraction
		graph export "nBrandsDay.pdf", as(pdf) replace
	restore
	
	preserve
		keep hhday brfl*
		
		reshape long brfl, i(hhday) j(purInDay) 
		bys hhday brfl: gen first = _n == 1
		bys hhday: gen numbrfls = sum(first)
		drop first
		reshape wide brfl numbrfls, i(hhday) j(purInDay)
		gen brfls = numbrfls16
		gen full = (brfl16 != . )
		replace brfls = brfls - 1 +full
		drop numbrfls*
		
		tabulate brfls
		histogram brfls, graphregion(color(white)) ///
				scheme(s2color) xtitle("Number of brands/flavors purchased in a day") ///
				fraction
		graph export "nBrFlsDay.pdf", as(pdf) replace
	restore


// Statistics for each week
use "purchasesIndices.dta", clear

	local vdrop purchase weightamount multipack nf displ feat xprice ustcoupons valcoup ///
		umancoupons valmcoupons purchased uweight purInDay hhday dayofweek

	local towide market stid n brand flav1 flav2 lf 
		 
	local hhvars hhsize hhincome restype mrace mheadwstat frace fheadwstat smom

	disp "`vdrop'"
	drop `vdrop'
		
	reshape wide `towide', i(hhweek) j(purInWeek)

	gen nPur = 0

	forvalues i = 1/22 {
		gen m`i' = (market`i' != .)
		//replace nPur = if(`i'
		
		replace nPur = nPur + m`i'
		drop m`i'
	}

	tabulate nPur
	histogram nPur, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of purchases per week") fraction
	graph export "nPurchasesWeek.pdf", as(pdf) replace

	keep hhweek brand*
	
	reshape long brand, i(hhweek) j(purInWeek) 
	bys hhweek brand: gen first = _n == 1
	bys hhweek: gen numbrands = sum(first)
	drop first
	reshape wide brand numbrands, i(hhweek) j(purInWeek)
	gen brands = numbrands22
	gen full = (brand22 != . )
	replace brands = brands - 1 +full
	drop numbrands*
	
	tabulate brands
	histogram brands, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of brands purchased in a week") fraction
	graph export "nBrandsWeek.pdf", as(pdf) replace
	

// Test brfl
use "purchasesIndices.dta", clear

	gen fl1 = flav1
	gen fl2 = flav2
	replace fl1 = 0 if flav1 == .
	replace fl2 = 0 if flav2 == .
	gen brfl = 10000*brand+100*fl1+fl2
	drop flav1 flav2
		
	keep hhday brfl

	duplicates drop

	bysort hhday: gen product = _n

	reshape wide brfl, i(hhday) j(product)

	forvalues i = 1/15 {
		gen nonnull`i' = 0
		replace nonnull`i' = (brfl`i'!=.)
	}

	egen nprod = rowtotal(nonnull*)

	tabulate nprod
	histogram nprod, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of brands purchased in a day") fraction
	graph export "nBrFlsDay2.pdf", as(pdf) replace

// Brands bought in any period
use "purchasesIndices.dta", clear

keep hhid brand

duplicates drop

bys hhid: gen count = _n
bys hhid: egen brands = max(count)

drop brand count

duplicates drop

tabulate brands
histogram brands, graphregion(color(white)) ///
		scheme(s2color) xtitle("Number of brands in whole period") fraction
graph export "nBrWhole.pdf", as(pdf) replace

// Flavors bought in any period
use "purchasesIndices.dta", clear

	gen flavor = 100*flav1+flav2
	keep hhid flavor

	duplicates drop

	bys hhid: gen count = _n
	bys hhid: egen flavors = max(count)

	drop flavor count

	duplicates drop

	tabulate flavors
	histogram flavors, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of flavors in whole period") fraction
	graph export "nFlsWhole.pdf", as(pdf) replace
	
	save "nFlsWhole.dta", replace

// Brand-flavors bought in any period
use "purchasesIndices.dta", clear

	gen brfl = 10000*brand+100*flav1+flav2
	//gen flavor = 100*flav1+flav2
	keep hhid brfl // flavor

	duplicates drop

	bys hhid: gen count = _n
	bys hhid: egen brfls = max(count)

	drop brfl count

	duplicates drop

	tabulate brfls
	histogram brfls, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of brand/flavors in whole period") ///
			fraction width(1)
	graph export "nBrFlsWhole.pdf", as(pdf) replace
	
	merge 1:1 hhid using "nFlsWhole.dta"
	
	drop _merge
	gen mayor = brfls < flavors

// Purchases in each market
use "purchasesIndices.dta", clear

	keep market week
	
	bys market week: gen n = _n
	bys market week: egen num = max(n)
	
	drop n
	duplicates drop
	
	tabulate num
	histogram num, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of purchases per market") ///
			fraction
	graph export "purMarket.pdf", as(pdf) replace

// Purchases of each brand in each market
use "purchasesIndices.dta", clear

	gen product = string(brand,"%02.0f") //+string(flav1,"%02.0f")+string(flav2,"%02.0f") //+string(lf)+string(nf)+string(weightamount/100,"%03.0f")
	keep market week product
	
	bys market week product: gen n = _n
	bys market week product: egen num = max(n)
	
	drop n
	duplicates drop
	
	tabulate num
	histogram num, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of purchases or each brand per market") ///
			fraction width(1)
	graph export "purBrandMarket.pdf", as(pdf) replace
	
// Purchases of each brand/flavor in each market
use "purchasesIndices.dta", clear

	gen product = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f") //+string(lf)+string(nf)+string(weightamount/100,"%03.0f")
	keep market week product
	
	bys market week product: gen n = _n
	bys market week product: egen num = max(n)
	
	drop n
	duplicates drop
	
	tabulate num
	histogram num, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of purchases or each brand/flavor per market") ///
			fraction width(1)
	graph export "purBrflMarket.pdf", as(pdf) replace
	
// Purchases of each product in each market
use "purchasesIndices.dta", clear

	gen product = string(brand,"%02.0f")+string(flav1,"%02.0f")+string(flav2,"%02.0f")+string(lf)+string(nf)+string(weightamount/100,"%03.0f")
	keep market week product
	
	bys market week product: gen n = _n
	bys market week product: egen num = max(n)
	
	drop n
	duplicates drop
	
	tabulate num
	histogram num, graphregion(color(white)) ///
			scheme(s2color) xtitle("Number of purchases or each product per market") ///
			fraction width(1)
	graph export "purProductMarket.pdf", as(pdf) replace
