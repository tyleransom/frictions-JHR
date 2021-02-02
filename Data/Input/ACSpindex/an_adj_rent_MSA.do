version 14.1
clear all
capture log close
set more off
log using "an_adj_rent_MSA.log", replace

use ACS2005raw, clear

run metarea_CBSA_mapping

drop if qacrehou > 0 | qbedroom > 0 | qbuilty2 > 0 | qkitchen > 0 | qrooms   > 0 | qunitsst > 0 | gq>1 // drop any imputed values and any group quarters
drop if metarea==0 // drop areas that are undisclosed

gen logrent = ln(rentgrs)

qui tab bedrooms, gen(bed_dum)
qui tab rooms, gen(room_dum)
qui tab builtyr2, gen(age_dum)
qui tab unitsstr, gen(units_dum)
qui tab plumbing, gen(plumbing_dum)
qui tab kitchen, gen(kitchen_dum)
qui tab acrehous, gen(acre_dum)
qui tab cbsa_code, gen(metro_dum)

local regressors bed_dum* room_dum* age_dum* units_dum* plumbing_dum* kitchen_dum* acre_dum*

count
areg logrent `regressors', absorb(cbsa_code) // metro_dum2-metro_dum265
matrix coeffs = e(b)'
predict fixed_effects, d

preserve
	collapse fixed_effects , by(cbsa_code)
	mkmat fixed_effects cbsa_code, matrix(FE) nomissing
restore

matrix accum A = `regressors' , noconstant deviations means(matmeans)
matrix matmeans = [matmeans,1]

** need to repmat "matmeans" 265 times ... or collapse the regressor means by location to get a 265 x K matrix ........
matrix matmeans = J(265,1,1)*matmeans

mata
zed = st_matrix("FE")
zed1 = zed[.,1]
zed2 = zed[.,2]
st_matrix("FE1",zed1)
st_matrix("FEnames",zed2)
end
matrix prices = matmeans*coeffs + FE1

svmat prices, names(locprice)
egen mean_price = mean(locprice)
sum mean_price
svmat FEnames, names(locnames)

preserve
keep locprice1 locnames1 mean_price
keep if locprice1<.
ren locnames1 cbsa_code
ren locprice1 grents
* replace grents = grents/mean_price*100 // deflate gross rents by mean location AFTER merging in with ACCRA data
drop mean_price
sum grents
l
save ../AccraData/sourcedata/gross_rents, replace
restore

log close
exit
