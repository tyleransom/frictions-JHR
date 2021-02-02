clear all
set mem 5g
capture log close
capture cd afs/econ.duke.edu/data/vjh3/TylerCPS/AccraData/finaldo
set more off
log using "cr_latlong.log", replace

*****************************************************************************
* Merge in county longitude and latitude
* using data from Census and MABLE/GeoCorr
*****************************************************************************
*============================================================================
* Import crosswalk between census place and county FIPS
*============================================================================
tempfile holder1
insheet using ../sourcedata/place2county.csv, comma clear

drop if placefp==99999
* collapse (first) county stabb cntyname placenm pop2000 afact, by(placefp state) // this gets rid of duplicate placefp codes
* now need to get rid of duplicate counties within the same place name

bys placefp state (county): egen maxplpop = max(pop2000)
drop if maxplpop~=pop2000

bys placenm state (county): egen maxnmpop = max(pop2000)
drop if maxnmpop~=pop2000


isid placefp state
isid placenm state
drop maxplpop maxnmpop
save `holder1', replace

*============================================================================
* Import dataset of census place and latitude/longitude
*============================================================================
insheet using ../sourcedata/placelongitude.csv, comma clear
* duplicates drop placefp state, force
duplicates drop placenm state, force
* isid placefp state
isid placenm state
* merge 1:1 placefp state using `holder1'


*============================================================================
* Merge the two together and collapse into one lat/long coord per county
*============================================================================
merge 1:1 placenm state using `holder1'
keep if _merge==3
drop _merge

** now collapse into one lat/long per county ... use highest-population place
bys county (placefp): egen maxpop = max(pop2000)
drop if maxpop~=pop2000

note: "Census place latitude/longitude downloaded from http://www.census.gov/geo/www/tiger/latlng.txt"

replace county = county-state*1000

save ../sourcedata/county_latlong_place, replace
log close
exit
