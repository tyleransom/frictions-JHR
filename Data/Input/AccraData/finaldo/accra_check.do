/** makeaccra.do

DESCRIPTION: This do-file takes ACCRA data from 1990 through 2008 and creates price index data for each county.
 It assumes Cobb-Douglas and allows substitution between 6 large product groups.   **/

clear all
set mem 100m
set more off

capture log close
log using accra_check.log, replace text

tempfile xwalk




*****************************************************************************
* Initial cleaning of ACCRA data
*****************************************************************************
insheet using /afs/econ.duke.edu/data/vjh3/TylerJMP/AccraData/accra19902008.csv, comma names clear
capture drop v16-v19
replace state_code = "46"    if state_code == "s46"
replace cbsa_code  = "10100" if cbsa_code  == "s10100"
replace city_code  = "700"   if city_code  == "s700"
destring cbsa_code, replace
recode cbsa_code (13644 = 47900)  // Bethesda-Gaithersburg-Frederick MD Metro Div.
recode cbsa_code (14484 = 14460)  // Boston-Quincy MA Metro Div.
recode cbsa_code (15764 = 14460)  // Cambridge-Newton-Framingham MA Metro Div.
recode cbsa_code (16974 = 16980)  // Chicago-Naperville-Joliet IL Metro Div.
recode cbsa_code (19124 = 19100)  // Dallas-Plano-Irving TX Metro Div.
recode cbsa_code (19804 = 19820)  // Detroit-Livonia-Dearborn MI Metro Div.
recode cbsa_code (20764 = 35620)  // Edison NJ Metro Div.
recode cbsa_code (22744 = 33100)  // Fort Lauderdale-Pompano Beach-Deerfield Beach FL Metro Div.
recode cbsa_code (23020 = 18880)  // Fort Walton Beach-Crestview-Destin FL Metro
recode cbsa_code (23104 = 19100)  // Fort Worth-Arlington TX Metro Div.
recode cbsa_code (29404 = 16980)  // Lake County-Kenosha County IL-WI Metropolitan Div.
recode cbsa_code (30540 = 45640)  // Lexington-Thomasville NC Micro
recode cbsa_code (31084 = 31100)  // Los Angeles-Long Beach-Glendale CA Metro Div.
recode cbsa_code (33124 = 33100)  // Miami-Miami Beach-Kendall FL Metro Div.
recode cbsa_code (35004 = 35620)  // Nassau-Suffolk NY Metro Div.
recode cbsa_code (35084 = 35620)  // Newark-Union NJ-PA Metro Div.
recode cbsa_code (35644 = 35620)  // New York-White Plains-Wayne NY-NJ Metro Div.
recode cbsa_code (36084 = 41860)  // Oakland-Fremont-Hayward CA Metro Div.
recode cbsa_code (37964 = 37980)  // Philadelphia PA Metro Div.
recode cbsa_code (41884 = 41860)  // San Francisco-San Mateo-Redwood City CA Metro Div.
recode cbsa_code (42044 = 31100)  // Santa Ana-Anaheim-Irvine CA Metro Div.
recode cbsa_code (42260 = 35840)  // Sarasota-Bradenton-Venice FL Metro
recode cbsa_code (42644 = 42660)  // Seattle-Bellevue-Everett WA Metro Div.
recode cbsa_code (45104 = 42660)  // Tacoma WA Metro Div.
recode cbsa_code (46940 = 42680)  // Vero Beach FL Metro
recode cbsa_code (47644 = 19820)  // Warren-Farmington Hills-Troy MI Metro Div.
recode cbsa_code (47850 = 47580)  // Warner Robins GA Metro
recode cbsa_code (47894 = 47900)  // Washington-Arlington-Alexandria DC-VA-MD-WV Metro Div.
recode cbsa_code (48424 = 33100)  // West Palm Beach-Boca Raton-Boynton Beach FL
recode cbsa_code (48864 = 37980)  // Wilmington DE-MD-NJ Metro Div.
recode cbsa_code (49220 = 32270)  // Wisconsin Rapids-Marshfield WI Micro

gen flag = metro_micro_name == "NonMetro US"
tab state_name, mi
tab urban_area_name if flag
l year quarter urban_area_name cbsa_code metro_micro_name if flag, sep(0)

drop if metro_micro_name == "Nonmetropolitan Census Area Canada" 
drop if metro_micro_name == "Saskatoon SK CMA"
drop if state_name       == "Virgin Islands"
drop if state_name       == "British Columbia"


keep metro_micro_name cbsa_code
collapse (first) metro_micro_name, by(cbsa_code)
sort cbsa_code
codebook cbsa_code

preserve
*****************************************************************************
* read in county crosswalk data for future merging 
* [source: Missouri's MABLE/Geocorr12 website]
*****************************************************************************
insheet using /afs/econ.duke.edu/data/vjh3/TylerJMP/AccraData/countyxwalk.csv, comma names clear
ren cbsa cbsa_code
isid county
replace county = county-state*1000

keep cbsaname cbsa_code
collapse (first) cbsaname, by(cbsa_code)

codebook cbsa_code

save `xwalk', replace
restore

merge 1:1 cbsa_code using `xwalk'

l if _merge==1, sep(0)
l if _merge==2, sep(0)

log close
exit
