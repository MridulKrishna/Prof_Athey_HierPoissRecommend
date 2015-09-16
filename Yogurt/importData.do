clear
clear matrix
clear mata
set mem 500m
set matsize 2000
set maxvar 5000
set more off

cd "/Users/casta2k/Dropbox/2015_Summer/Athey/Matlab_AtheyCastillo/Yogurt"

import delimited "bigpurchfile.dat", delimiter(space, collapse) asdouble

drop v1

rename v2 market
rename v3 hhid
label variable hhid "Household id"
rename v4 stid
label variable stid "Store id?"
rename v5 week
rename v6 n
label variable n "??"
rename v7 dayofweek
rename v8 brand
rename v9 weightamount
rename v10 multipack
rename v11 flav1
rename v12 flav2
rename v13 lf
label variable lf "Low fat?"
rename v14 nf
label variable nf "No fat?"
rename v15 displ
label variable displ "???"
rename v16 feat
label variable feat "Featured in ad"
rename v17 xprice
label variable xprice "??"
rename v18 ustcoupons
label variable ustcoupons "??"
rename v19 valcoup
label variable valcoup "??"
rename v20 umancoupons
label variable umancoupons "??"
rename v21 valmcoupons
label variable valmcoupons "??"
rename v22 purchased
label variable purchased "??"
rename v23 uweight
label variable uweight "Total weight"
rename v24 hhsize
label variable hhsize "Household size"
rename v25 hhincome
label variable hhincome "Household income"
rename v26 restype
label variable restype "Residence type?"
rename v27 mrace
label variable mrace "Male/mother race?"
rename v28 mheadwstat
label variable mheadwstat "Male/mother ... stat?"
rename v29 frace
label variable frace "Female/father race?"
rename v30 fheadwstat
label variable fheadwstat "Female/father ... stat?"
rename v31 smom
label variable smom "Single mother?"

save "purchases.dta", replace
