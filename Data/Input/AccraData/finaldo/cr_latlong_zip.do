clear all
set mem 5g
capture log close
capture cd afs/econ.duke.edu/data/vjh3/TylerCPS/AccraData/finaldo
set more off
log using "cr_latlong_zip.log", replace

tempfile holder1
*****************************************************************************
* Merge in county longitude and latitude
* using data from Census and MABLE/GeoCorr
*****************************************************************************
*============================================================================
* Import dataset of census place and latitude/longitude
*============================================================================
insheet using ../sourcedata/ziplongitude.csv, comma clear
keep zip latitude longitude stabb
isid zip
save `holder1', replace

*============================================================================
* Import crosswalk between census place and county FIPS
*============================================================================
insheet using ../sourcedata/zip2county.csv, comma clear
genera state = floor(county/1000)
rename county FIPS
genera county = FIPS - state*1000
destring zip, replace force
codebook zip
codebook FIPS
* isid FIPS
* isid zip
* drop maxzippop

*============================================================================
* Merge the two together and collapse into one lat/long coord per county
*============================================================================
merge m:1 zip using `holder1'
mdesc zip
mdesc FIPS
codebook FIPS

* l if _merge==1

bys FIPS (zip): egen numvalidzips = count(zip)
l if numvalidzips==0
bys FIPS (zip): egen numvalidlats = count(latitude)
tab cntyname if numvalidlats==0
l if numvalidlats==0

* keep if _merge==3
drop _merge

** hand-code latitude and longitude for 8 counties in Georgia:
   * Calhoun GA 31°26′22″N 84°43′29″W  // 31.43944 -84.72472
      * Clay GA 31°36′51″N 85°2′54″W   // 31.61417 -85.04833
   * Decatur GA 30°54′17″N 84°34′16″W  // 30.90472 -84.57111
     * Early GA 31°22′36″N 84°56′2″W   // 31.37667 -84.93389
    * Miller GA 31°10′23″N 84°43′43″W  // 31.17306 -84.72861
   * Quitman GA 31°53′02″N 85°06′05″W  // 31.88389 -85.10139
  * Randolph GA 31°46′15″N 84°47′37″W  // 31.77083 -84.79361
  * Seminole GA 31°2′27″N 84°52′42″W   // 31.04083 -84.87889
  *    Baker GA 31°19′0″N 84°20′22″W   // 31.31680 -84.33955
  *    Grady GA 30°53′0″N 84°13′0″W    // 30.87780 -84.20890  
  *  Terrell GA 31°46′26″N 84°26′27″W  // 31.77397 -84.44087  
* Broomfield CO 39°55′55″N 105°03′57″W // 39.93182 -105.065919
* Hoonah-Angoon Census Area (AK)  58°6′34″N 135°26′11″W //  58.109435 -135.436349
* Skagway Municipality (AK)       59°27′30″N 135°18′50″W // 59.458333 -135.313889 
* 2195 (AK)???
* 2198 (AK)???
* 2275 (AK)???


replace latitude = 31.43944 if cntyname=="Calhoun GA"
replace latitude = 31.61417 if cntyname=="Clay GA"
replace latitude = 30.90472 if cntyname=="Decatur GA"
replace latitude = 31.37667 if cntyname=="Early GA"
replace latitude = 31.17306 if cntyname=="Miller GA"
replace latitude = 31.88389 if cntyname=="Quitman GA"
replace latitude = 31.77083 if cntyname=="Randolph GA"
replace latitude = 31.04083 if cntyname=="Seminole GA"
replace latitude = 31.31680 if cntyname=="Baker GA"
replace latitude = 30.87780 if cntyname=="Grady GA"
replace latitude = 31.77397 if cntyname=="Terrell GA"
replace latitude = 39.93182 if cntyname=="Broomfield CO"
replace latitude = 58.10944 if cntyname=="Hoonah-Angoon Census Area"
replace latitude = 59.45833 if cntyname=="Skagway Municipality"

replace longitude = -84.72472 if cntyname=="Calhoun GA"
replace longitude = -85.04833 if cntyname=="Clay GA"
replace longitude = -84.57111 if cntyname=="Decatur GA"
replace longitude = -84.93389 if cntyname=="Early GA"
replace longitude = -84.72861 if cntyname=="Miller GA"
replace longitude = -85.10139 if cntyname=="Quitman GA"
replace longitude = -84.79361 if cntyname=="Randolph GA"
replace longitude = -84.87889 if cntyname=="Seminole GA"
replace longitude = -84.33955 if cntyname=="Baker GA"
replace longitude = -84.20890 if cntyname=="Grady GA"
replace longitude = -84.44087 if cntyname=="Terrell GA"
replace longitude = -105.0659 if cntyname=="Broomfield CO"
replace longitude = -135.4364 if cntyname=="Hoonah-Angoon Census Area"
replace longitude = -135.3139 if cntyname=="Skagway Municipality"

l if FIPS==13007
l if FIPS==8014

** now collapse into one lat/long per county ... use highest-population place
bys FIPS (zip): egen maxzippop = max(pop2000)
drop if maxzippop~=pop2000

l if FIPS==13007
drop if mi(FIPS) | mi(zip) | mi(latitude) | mi(longitude)
l if FIPS==13007

* l if numvalidlats==0
* l if maxzippop>=.

codebook FIPS



* ** now collapse into one lat/long per county ... use highest-population place
* bys county (placefp): egen maxpop = max(pop2000)
* drop if maxpop~=pop2000

note: "ZIP latitude/longitude downloaded from http://federalgovernmentzipcodes.us/free-zipcode-database-Primary.csv"

* replace county = county-state*1000

save ../sourcedata/county_latlong_zip, replace
log close
exit
