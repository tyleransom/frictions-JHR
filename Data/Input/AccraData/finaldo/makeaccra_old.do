/** makeaccra.do

DESCRIPTION: This do-file takes ACCRA data from 1990 through 2008 and creates price index data for each county.
 It assumes Cobb-Douglas and allows substitution between 6 large product groups. 
 
DATA SOURCES: MSA GDP comes from the BEA; ACCRA data comes from ACCRA; county unemployment data comes from BLS (bls.gov/lau); 
      county income per worker by industry comes from BEA **/

clear all
set mem 100m
set more off

capture log close
log using makeaccra.log, replace text

tempfile xwalk accra_quarterly accra_weights accra_annual msa_gdp msa_gdp_per_capita temp1 temp2 temp3 counties countiesRural csize_impute region_impute cityp ruralp LLM countyunemp countyunemplong MSAunemplong MSAlatlong countylatlong county_grents city_grents countyLandUse MSAlandUseLong countyLandUseLong

global pathstem /afs/econ.duke.edu/data/vjh3/TylerJMP/AccraData/
global source   sourcedata/
global final    finaldata/



*****************************************************************************
* Initial cleaning of ACCRA data
*****************************************************************************
insheet using ${pathstem}${source}accra19902008.csv, comma names clear
capture noisily drop v16-v19
* tab cbsa_code
drop if metro_micro_name == "Nonmetropolitan Census Area Canada" 
drop if metro_micro_name == "Saskatoon SK CMA"
drop if metro_micro_name == "Vancouver BC"
drop if state_name       == "Virgin Islands"
drop if state_name       == "British Columbia"

*============================================================================
* change quarters from strings to numbers
*============================================================================
replace quarter = "1" if quarter=="1Q"
replace quarter = "2" if quarter=="2Q"
replace quarter = "3" if quarter=="3Q"
replace quarter = "4" if quarter=="4Q"
destring quarter, replace
drop if mi(quarter)
sort year quarter cbsa_code

*============================================================================
* merge in quarterly ACCRA weights
*============================================================================
preserve
	insheet using ${pathstem}${source}accraweights19902008.csv, comma names clear
	drop if mi(quarter)
	sort year quarter
	local weight_vars weight_composite_index weight_grocery_items weight_housing weight_utilities weight_transportation weight_health_care weight_misc_goods_services
	* foreach x in `weight_vars' {
		* replace `x' = substr( `x',1,length( `x')-1)
		* destring `x', replace
		* replace `x' = `x'/100
	* }
	save `accra_weights', replace
restore
merge m:1 year quarter using `accra_weights'
capture noisily drop _merge

*============================================================================
* make sure location codes are string variables
*============================================================================
replace state_code = "46"    if state_code == "s46"
replace cbsa_code  = "10100" if cbsa_code  == "s10100"
replace city_code  = "700"   if city_code  == "s700"
destring cbsa_code, replace

*============================================================================
* fix some errors in the product variables
*============================================================================
replace transportation = "86.7"  if transportation== "86.7."
replace transportation = "119.9" if transportation== "119..9"
destring transportation, replace

*============================================================================
* clean some incorrect city codes
*============================================================================
recode  cbsa_code (13644 = 47900)  // Bethesda-Gaithersburg-Frederick MD Metro Div.
recode  cbsa_code (14484 = 14460)  // Boston-Quincy MA Metro Div.
recode  cbsa_code (15764 = 14460)  // Cambridge-Newton-Framingham MA Metro Div.
recode  cbsa_code (16974 = 16980)  // Chicago-Naperville-Joliet IL Metro Div.
recode  cbsa_code (19124 = 19100)  // Dallas-Plano-Irving TX Metro Div.
recode  cbsa_code (19804 = 19820)  // Detroit-Livonia-Dearborn MI Metro Div.
recode  cbsa_code (20764 = 35620)  // Edison NJ Metro Div.
recode  cbsa_code (22744 = 33100)  // Fort Lauderdale-Pompano Beach-Deerfield Beach FL Metro Div.
recode  cbsa_code (23020 = 18880)  // Fort Walton Beach-Crestview-Destin FL Metro
recode  cbsa_code (23104 = 19100)  // Fort Worth-Arlington TX Metro Div.
recode  cbsa_code (29404 = 16980)  // Lake County-Kenosha County IL-WI Metropolitan Div.
recode  cbsa_code (30540 = 45640)  // Lexington-Thomasville NC Micro
recode  cbsa_code (31084 = 31100)  // Los Angeles-Long Beach-Glendale CA Metro Div.
recode  cbsa_code (33124 = 33100)  // Miami-Miami Beach-Kendall FL Metro Div.
recode  cbsa_code (35004 = 35620)  // Nassau-Suffolk NY Metro Div.
recode  cbsa_code (35084 = 35620)  // Newark-Union NJ-PA Metro Div.
recode  cbsa_code (35644 = 35620)  // New York-White Plains-Wayne NY-NJ Metro Div.
recode  cbsa_code (36084 = 41860)  // Oakland-Fremont-Hayward CA Metro Div.
recode  cbsa_code (37964 = 37980)  // Philadelphia PA Metro Div.
recode  cbsa_code (41884 = 41860)  // San Francisco-San Mateo-Redwood City CA Metro Div.
recode  cbsa_code (42044 = 31100)  // Santa Ana-Anaheim-Irvine CA Metro Div.
recode  cbsa_code (42260 = 35840)  // Sarasota-Bradenton-Venice FL Metro
recode  cbsa_code (42644 = 42660)  // Seattle-Bellevue-Everett WA Metro Div.
recode  cbsa_code (45104 = 42660)  // Tacoma WA Metro Div.
recode  cbsa_code (46940 = 42680)  // Vero Beach FL Metro
recode  cbsa_code (47644 = 19820)  // Warren-Farmington Hills-Troy MI Metro Div.
recode  cbsa_code (47850 = 47580)  // Warner Robins GA Metro
recode  cbsa_code (47894 = 47900)  // Washington-Arlington-Alexandria DC-VA-MD-WV Metro Div.
recode  cbsa_code (48424 = 33100)  // West Palm Beach-Boca Raton-Boynton Beach FL
recode  cbsa_code (48864 = 37980)  // Wilmington DE-MD-NJ Metro Div.
recode  cbsa_code (49220 = 32270)  // Wisconsin Rapids-Marshfield WI Micro
replace cbsa_code = 28140 if urban_area_name == "Clinton MO"
replace cbsa_code = 17660 if urban_area_name == "Coeur d'Alene ID"
replace cbsa_code = 19140 if urban_area_name == "Dalton GA"
replace cbsa_code = 24660 if urban_area_name == "Eden NC"
* replace cbsa_code = 51111 if urban_area_name == "Glenwood Springs CO"
* replace cbsa_code = 51111 if urban_area_name == "Gunnison CO"
replace cbsa_code = 36740 if urban_area_name == "Lake County FL"
* replace cbsa_code = 52222 if urban_area_name == "Lexington-Buena Vista-Rockbridge VA"
* replace cbsa_code = 53333 if urban_area_name == "Lincoln County OR"
* replace cbsa_code = 54444 if urban_area_name == "Marion-McDowell County NC"
* replace cbsa_code = 55555 if urban_area_name == "McCormick County SC"
* replace cbsa_code = 56666 if urban_area_name == "Nevada MO"
* replace cbsa_code = 57777 if urban_area_name == "Pikeville-Pike County KY"
* replace cbsa_code = 58888 if urban_area_name == "Pryor Creek OK"
replace cbsa_code = 40760 if urban_area_name == "Ruidoso NM"
* replace cbsa_code = 59999 if urban_area_name == "Sullivan County NY"
replace cbsa_code = 11700 if urban_area_name == "Waynesville-Haywood County NC"

*============================================================================
* convert year-quarter pairing to quarterly date variable
*============================================================================
genera date = yq(year,quarter)
format date %tq
// drop year quarter // state_code city_code
isid date state_code cbsa_code city_code

*============================================================================
* average ACCRA over cities with multiple CBSA codes:
*============================================================================
collapse (mean) date composite_index grocery_items housing utilities transportation health_care misc_goods_services weight*, by(cbsa_code year quarter)
isid date cbsa_code

*============================================================================
* average ACCRA over quarters within cities to get annualized version:
*============================================================================
preserve
	collapse (mean) composite_index grocery_items housing utilities transportation health_care misc_goods_services weight*, by(cbsa_code year)
	rename  composite_index     accra
	rename  grocery_items       grocery
	rename  transportation      transport
	rename  health_care         health
	rename  misc_goods_services misc
	local   vars accra grocery housing utilities transport health misc weight_composite_index weight_grocery_items weight_housing weight_utilities weight_transportation weight_health_care weight_misc_goods_services 
	reshape wide `vars', i(cbsa_code) j(year) // if you want to undo this: reshape long `vars', i(state_code cbsa_code city_code ) j(date)
	save `accra_annual'
restore

*============================================================================
* rename vars and reshape wide for ease of merging with other datasets
*============================================================================
rename  composite_index     accra
rename  grocery_items       grocery
rename  transportation      transport
rename  health_care         health
rename  misc_goods_services misc
drop    year quarter
reshape wide `vars', i(cbsa_code) j(date) // if you want to undo this: reshape long `vars', i(state_code cbsa_code city_code ) j(date)
save `accra_quarterly', replace





*****************************************************************************
* read in county crosswalk data for future merging 
* [source: Missouri's MABLE/Geocorr12 website]
*****************************************************************************
insheet using ${pathstem}${source}countyxwalk.csv, comma names clear
isid county
replace county = county-state*1000
* replace cbsa   = 51111 if cntyname=="Garfield CO"
* replace cbsa   = 51111 if cntyname=="Gunnison CO"
* replace cbsa   = 52222 if cntyname=="Rockbridge VA"
* replace cbsa   = 53333 if cntyname=="Lincoln OR"
* replace cbsa   = 54444 if cntyname=="McDowell NC"
* replace cbsa   = 55555 if cntyname=="McCormick SC"
* replace cbsa   = 56666 if cntyname=="Vernon MO"
* replace cbsa   = 57777 if cntyname=="Pike KY"
* replace cbsa   = 58888 if cntyname=="Mayes OK"
* replace cbsa   = 59999 if cntyname=="Sullivan NY"
* tab cbsa
*============================================================================
* make a file with all counties (but don't merge this with the ACCRA or GDP data just yet)
*============================================================================
preserve
	keep cbsa state county stabb cntyname cbsaname
	isid county state
	save `counties', replace
restore

preserve // make one to merge the data for just the rural counties
	keep cbsa state county stabb cntyname cbsaname pop00 pop10 landarea
	rena cbsa cbsa_code
	rena pop00 popR00
	rena pop10 popR10
	rena landarea landareaR
	isid county state
	save `countiesRural', replace
restore

*============================================================================
* prepare to merge in county unemployment data for aggregation to city level
*============================================================================
preserve
	insheet using ${pathstem}${source}county_unemp.csv, comma names clear
	drop cntyname urate
	drop if year==2011
	drop if state==72 // throw out Puerto Rico
	destring labor_force, force replace
	destring employed, force replace
	destring unemployed, force replace
	isid county state year
	save `countyunemp', replace
restore

preserve
	merge 1:m county state using `countyunemp'
	drop if cbsa==99999
	* tab1 state county if _merge==2
	keep if _merge==3
	collapse (sum) labor_force employed unemployed, by(cbsa cbsaname year)
	rename cbsa cbsa_code
	sort cbsa_code year
	isid cbsa_code year
	save `MSAunemplong', replace
restore

preserve
	merge 1:m county state using `countyunemp'
	drop if cbsa~=99999
	keep if _merge==3
	sort state county year
	ren labor_force labor_force_rural
	ren employed   employed_rural
	ren unemployed unemployed_rural
	keep labor_force_rural employed_rural unemployed_rural state county year
	isid state county year
	save `countyunemplong', replace
restore


*============================================================================
* prepare to merge in county land use data for aggregation to city level
*============================================================================
preserve
	insheet using ${pathstem}${source}PctUrbanRural_County.csv, comma names clear
	drop countyname statename
	* key variables: area_cou area_ua --- aggregate these to MSA, then generate MSA-level areapct_ua
	drop if state==72 // throw out Puerto Rico
	isid county state
	save `countyLandUse', replace
restore

preserve
	merge 1:m county state using `countyLandUse'
	drop if cbsa==99999
	* tab1 state county if _merge==2
	keep if _merge==3
	collapse (sum) area_cou area_ua, by(cbsa cbsaname)
	rename cbsa cbsa_code
	sort cbsa_code
	isid cbsa_code
	save `MSAlandUseLong', replace
restore

preserve
	merge 1:m county state using `countyLandUse'
	drop if cbsa~=99999
	keep if _merge==3
	sort state county
	ren area_cou area_cou_rural
	ren area_ua  area_ua_rural
	keep area_cou_rural area_ua_rural state county
	isid state county
	save `countyLandUseLong', replace
restore


*============================================================================
* prepare to merge in county latitude/longitude data for aggregation to city level
*============================================================================
preserve
	drop if cbsa==99999
	capture drop _merge
	gen FIPS = state*1000+county
	merge m:1 state county using ${pathstem}${source}county_latlong_zip, keepusing(latitude longitude state county)
	keep if _merge==3
	keep     cbsa cbsaname latitude longitude pop00
	bys cbsa: egen maxcpop = max(pop00)
	keep if pop00==maxcpop
	drop pop00 maxcpop
	rename   cbsa cbsa_code
	sort     cbsa_code
	save `MSAlatlong', replace
restore

preserve
	drop if cbsa~=99999
	capture drop _merge
	gen FIPS = state*1000+county
	merge m:1 state county using ${pathstem}${source}county_latlong_zip, keepusing(latitude longitude state county)
	replace latitude = 58.10944 if cntyname=="Hoonah-Angoon Census Area"
	replace latitude = 59.45833 if cntyname=="Skagway Municipality"
	replace longitude = -135.4364 if cntyname=="Hoonah-Angoon Census Area"
	replace longitude = -135.3139 if cntyname=="Skagway Municipality"
	keep if _merge==3 | cntyname=="Hoonah-Angoon Census Area" | cntyname=="Skagway Municipality"
	ren latitude Rlatitude
	ren longitude Rlongitude
	drop _merge
	save `countylatlong', replace
restore


*============================================================================
* collapse county-level data (pop 2000, pop 2010, land area) into city-level data
*============================================================================
drop if cbsa==99999
keep     cbsa cbsaname county pop00 pop10 landarea
collapse (sum) pop00 pop10 landarea, by(cbsa cbsaname)
generate dens00 = pop00/landarea
generate dens10 = pop10/landarea
rename   cbsa cbsa_code
sort     cbsa_code
save `xwalk', replace





*****************************************************************************
* read in crosswalk of CBSA code to IPUMS city code (MSA/CMSA/PMSA code)
* [source: Missouri's MABLE/Geocorr12 website]
*****************************************************************************
preserve
	insheet using ${pathstem}${source}county_cps_msa_xwalk.csv, comma names clear
	replace county = county-state*1000
	destring pmsa, replace force
	drop if msacmsa==9999 | mi(county)
	
	* Do some hand-coding of metarea to more seamlessly match up with the ACS
	recode msacmsa (2320 = 2310) // fix El Paso-Elkhart-Elmira-Enid mix-up
	recode msacmsa (2330 = 2320) // fix El Paso-Elkhart-Elmira-Enid mix-up
	recode msacmsa (2985 = 2990) // fix "Grands" mixup
	recode msacmsa (2995 = 3010) // fix "Grands" mixup
	recode msacmsa (3285 = 3300) // fix Hattiesburg mixup
	recode msacmsa (3600 = 3595) // fix Jacksonville mixup
	recode msacmsa (7480 = 7470) // fix Santa Barbara mixup
	replace msacmsa = floor(msacmsa/10) if mi(pmsa)
	replace msacmsa = floor(test1/10) if ~mi(pmsa)
	ren msacmsa metarea
	
	* Do some more hand-coding of metarea to more seamlessly match up with the ACS
	recode  metarea (416 456 645 = 112) // group Boston suburbs with Boston
	recode  metarea (130         = 131) // fix Burlington, VT
	recode  metarea (296         = 160) // group Gary, IN with Chicago
	recode  metarea (280         = 192) // group Ft Worth with Dallas
	recode  metarea (114         = 336) // group Brazoria with Houston
	recode  metarea (594         = 448) // group Orange County with LA
	recode  metarea (501         = 848) // group central NJ counties together
	recode  metarea (364 564 538 = 560) // group Newark, Hoboken, and suburban Long Island with New York
	recode  metarea (872 577     = 736) // group Oakland and Vallejo with San Francisco
	
	* merge
	merge m:1 metarea using ${pathstem}${source}gross_rents
	keep if _merge==3
	keep  metarea state county msaname afact grents
	
	* concatenate so that each county has one price index (since there are some counties that contribute to multiple PMSAs and hence have multiple price indices)
	collapse (mean) grents (first) metarea msaname [aweight = afact], by(state county)
	save `city_grents', replace
restore

preserve
	* merge city-level gross rents into county-level file
	use `city_grents', clear // variables: state, county, grents, metarea, msaname
	merge 1:1 state county using `counties', keepusing (state county cntyname cbsaname)
	drop _merge
	ren grents grents_urban
	
	* merge state-level gross rents into county-level file and assign state rents to locations where city rents are unobserved
	merge m:1 state using ${pathstem}${source}gross_rents_state
	replace grents = grents_urban if ~mi(grents_urban)
	l state county grents_urban metarea msaname cntyname cbsaname grents if _merge~=3
	keep  state county cntyname cbsaname metarea msaname grents
	save `county_grents', replace
restore




*****************************************************************************
* now merge with ACCRA and GDP datasets
*****************************************************************************
merge 1:1 cbsa_code using `accra_annual' // merge==2 is NonMetro US
capture noisily drop _merge
mdesc accra*

preserve
	insheet using ${pathstem}${source}msagdp.csv, comma names clear
	* tab cbsa
	ren cbsa cbsa_code
	sort cbsa_code
	forvalues X = 2001/2010 {
		replace gdp`X' = gdp`X'*1000000 // to convert gdp into dollars (instead of millions of dollars)
	}
	save `msa_gdp', replace
restore

// NOTE: SOME CITIES (E.G. PALM COAST, FL) APPEAR IN MSAGDP.CSV THAT DON'T APPEAR IN MSAGDP_PER_CAPITA.CSV (AND VICE VERSA)!!!
preserve
	insheet using ${pathstem}${source}msagdp_per_capita.csv, comma names clear
	* tab cbsa
	ren cbsa cbsa_code
	sort cbsa_code
	forvalues X = 2001/2010 {
		ren gdp`X' gdpPerCapita`X' // gdpPerCapita is in dollars, so this is OK
	}
	save `msa_gdp_per_capita', replace
restore

merge 1:1 cbsa_code using `msa_gdp'
drop if _merge==2 // this is the US average MSA GDP
capture noisily drop _merge
mdesc accra*

merge 1:1 cbsa_code using `msa_gdp_per_capita'
drop if _merge==2 // this is the US average MSA GDP
capture noisily drop _merge
mdesc accra*

merge 1:1 cbsa_code using `MSAlatlong'
keep if _merge==3

order cbsa_code cbsaname pop* land* gdp* lat* long*
drop _merge
save ${pathstem}${final}city_chars.dta, replace






*****************************************************************************
* now merge this with the county-level data (filled in with city-level values)
*****************************************************************************
use `counties', clear
isid county state
ren cbsa cbsa_code

merge m:1 cbsa_code using ${pathstem}${final}city_chars.dta
drop if _merge==2 // these are listings in the ACCRA data that don't exist in the other crosswalks
capture noisily drop _merge
isid county state

mdesc accra*

order cbsa_code cbsaname state stabb county cntyname pop* land* den* gdp*

merge 1:1 county state using `countiesRural'
replace pop10    = popR10         if cbsa_code==99999
replace pop00    = popR00         if cbsa_code==99999
replace landarea = landareaR      if cbsa_code==99999
replace dens00   = pop00/landarea if cbsa_code==99999
replace dens10   = pop10/landarea if cbsa_code==99999
drop popR10 popR00 landareaR _merge

merge m:1 county state using `countylatlong'
replace latitude  = Rlatitude  if cbsa_code==99999
replace longitude = Rlongitude if cbsa_code==99999
drop _merge

merge m:1 county state using `county_grents'
keep if _merge==3
drop _merge




*****************************************************************************
* impute missing ACCRAs 
*****************************************************************************
* Flag rural counties for imputation purposes
gen rural1 = mi(gdp2001) // missing MSA GDPs
gen shorty = substr(cbsaname,-29,5) // flag last 29 characters from cbsaname
gen rural2 = shorty == "Micro" | mi(shorty)

gen rural = rural1

foreach x in pop00 pop10 dens00 dens10 {
	ren `x' `x'c
}

*============================================================================
* add City Size category and Census Region to data (for imputation purposes)
*============================================================================
generat csize = .
replace csize = 5 if ~rural & inrange(pop00,1.5e6,2e8) // 1.5mil+
replace csize = 4 if ~rural & inrange(pop00,5e5,1.5e6) // 500k - 1.5mil
replace csize = 3 if ~rural & inrange(pop00,2e5,5e5)   // 200k - 500k
replace csize = 2 if ~rural & pop00<2e5                // <200k
replace csize = 1 if rural
assert  ~mi(csize)

lab def vlcsize 1 "Rural" 2 "<200k" 3 "200k - 500k" 4 "500k - 1.5mil" 5 "1.5mil+"
lab def vlregio 1 "Northeast" 2 "South" 3 "Midwest" 4 "West" 5 "Pacific"

generat region = .
replace region = 1 if inlist(state,23,33,50,25,44,9,36,42,34) // ME,NH,VT,MA,RI,CT,NY,PA,NJ
replace region = 2 if inlist(state,10,24,11,51,54,21,47,37,45,13,1,28,12,5,22,40,48) // DE,MD,DC,VA,WV,KY,TN,NC,SC,GA,AL,MS,FL,AR,LA,OK,TX
replace region = 3 if inlist(state,39,18,26,55,17,29,20,31,19,27,46,38) // OH,IN,MI,WI,IL,MO,KS,NE,IA,MN,SD,ND
replace region = 4 if inlist(state,35,8,56,30,16,49,4,32,6,41,53) // NM,CO,WY,MT,ID,UT,AZ,NV,CA,OR,WA
replace region = 5 if inlist(state,2,15) // AK,HI
assert  ~mi(region)

lab val csize  vlcsize
lab val region vlregio

*fix so Honolulu can be imputed with Anchorage:
replace csize = 3 if csize>3 & region==5

reshape long gdpPerCapita gdp accra grocery housing utilities transport health misc weight_composite_index weight_grocery_items weight_housing weight_utilities weight_transportation weight_health_care weight_misc_goods_services, i(state county) j(year)

ren weight_grocery_items       sgrt
ren weight_housing             shot
ren weight_utilities           sutt
ren weight_transportation      shet
ren weight_health_care         strt
ren weight_misc_goods_services smit

preserve
	collapse Maccra=accra Mgrocery=grocery Mhousing=housing Mutilities=utilities Mtransport=transport Mhealth=health Mmisc=misc Msgrt=sgrt Mshot=shot Msutt=sutt Mshet=shet Mstrt=strt Msmit=smit, by(state csize year)
	save `csize_impute', replace
restore

preserve
	collapse Raccra=accra Rgrocery=grocery Rhousing=housing Rutilities=utilities Rtransport=transport Rhealth=health Rmisc=misc Rsgrt=sgrt Rshot=shot Rsutt=sutt Rshet=shet Rstrt=strt Rsmit=smit, by(region csize year)
	save `region_impute', replace
restore

merge m:1 state  csize year using `csize_impute' , nogen
merge m:1 region csize year using `region_impute', nogen

gen qAccraState  = 0
gen qAccraRegion = 0
foreach var in accra grocery housing utilities transport health misc sgrt shot sutt shet strt smit {
	qui replace qAccraState = 1 if mi(`var') & ~mi(M`var')
	qui replace qAccraRegion= 1 if mi(M`var')& ~mi(R`var')
	qui replace `var'  = M`var' if mi(`var') & ~mi(M`var')
	qui replace `var'  = R`var' if mi(`var') & ~mi(R`var')
}

drop R* M* // drop imputing means


* generat AREA=0  // for use with Baum-Snow and Pavan (2012) data
* replace AREA=1  if cbsa_code==11260 // loc=="Anchorage"
* replace AREA=2  if cbsa_code==12060 // loc=="Atlanta"
* replace AREA=3  if cbsa_code==14460 // loc=="Boston"
* replace AREA=4  if cbsa_code==16980 // loc=="Chicago"
* replace AREA=5  if cbsa_code==17140 // loc=="Cincinnati"
* replace AREA=6  if cbsa_code==17380 // loc=="Cleveland"
* replace AREA=7  if cbsa_code==19100 // loc=="Dallas"
* replace AREA=8  if cbsa_code==19740 // loc=="Denver"
* replace AREA=9  if cbsa_code==19820 // loc=="Detroit"
* replace AREA=10 if cbsa_code==26180 // loc=="Honolulu"
* replace AREA=11 if cbsa_code==26420 // loc=="Houston"
* replace AREA=12 if cbsa_code==28140 // loc=="Kansas City"
* replace AREA=13 if cbsa_code==31100 // loc=="Los Angeles"
* replace AREA=14 if cbsa_code==33100 // loc=="Miami"
* replace AREA=18 if cbsa_code==33340 // loc=="Milwaukee"
* replace AREA=19 if cbsa_code==33460 // loc=="Minneapolis"
* replace AREA=20 if cbsa_code==35620 // loc=="New York"
* replace AREA=23 if cbsa_code==37980 // loc=="Philly"
* replace AREA=24 if cbsa_code==38060 // loc=="Phoenix"
* replace AREA=25 if cbsa_code==38300 // loc=="Pittsburgh"
* replace AREA=26 if cbsa_code==38900 // loc=="Portland"
* replace AREA=27 if cbsa_code==41740 // loc=="San Diego"
* replace AREA=28 if cbsa_code==41860 // loc=="San Francisco"
* replace AREA=29 if cbsa_code==42660 // loc=="Seattle"
* replace AREA=34 if cbsa_code==41180 // loc=="St Louis"
* replace AREA=35 if cbsa_code==45300 // loc=="Tampa"
* replace AREA=37 if cbsa_code==47900 // loc=="Washington"

* replace AREA=21 if ~rural & AREA==0 & region==1 // loc=="NorthEast A"
* replace AREA=22 if  rural & AREA==0 & region==1 // loc=="NorthEast B/C"

* replace AREA=31 if ~rural & AREA==0 & region==2 // loc=="South A"
* replace AREA=32 if ~rural & AREA==0 & region==2 // loc=="South B/C"
* replace AREA=33 if  rural & AREA==0 & region==2 // loc=="South D"

* replace AREA=15 if ~rural & AREA==0 & region==3 // loc=="Midwest A"
* replace AREA=16 if ~rural & AREA==0 & region==3 // loc=="Midwest B/C"
* replace AREA=17 if  rural & AREA==0 & region==3 // loc=="Midwest D"

* replace AREA=38 if ~rural & AREA==0 & region==4 // loc=="West A"
* replace AREA=39 if  rural & AREA==0 & region==4 // loc=="West B/C"
* * replace AREA=36 if rural & AREA==0 & region==1 // loc=="US"
* * replace AREA=30 if rural & AREA==0 & region==1 // loc=="Size D"

mdesc accra grocery housing grents utilities transport health misc sgrt shot sutt shet strt smit if year<2009 // I don't have any ACCRA data after 2008

*****************************************************************************
* Form price index based on shares from ACCRA (shares are same across all 
* locations and taken from CEX, though I'm not sure how exactly)
*****************************************************************************
*** ACCRA Data
gen CPI_Ta=(grocery)^sgrt*(housing)^shot*(utilities)^sutt*(health)^shet*(transport)^strt*(misc)^smit
*** Accra Data but with housing component from the ACS following Winters (2009)
gen CPI_Tb=(grocery)^sgrt*(grents)^(shot+sutt)*(health)^shet*(transport)^strt*(misc)^smit // include share of utilities in housing because rents are gross
* gen CPI_Tavr=(GR)^sgrt*(var2*bvalrent)^shot*(UT)^sutt*(HE)^shet*(TR)^strt*(MI)^smit
* *** Only Spatial differentiation is from the housing component from the census
* gen CPI_Tvr=(GRn)^sgrt*(var2*bvalrent)^shot*(UTn)^sutt*(HEn)^shet*(TRn)^strt*(MIn)^smit
* *** Only spatial diff from housing component from the census & grocery from Weinstein/Handbury
* gen CPI_Tgvr=(GR2)^sgrt*(var2*bvalrent)^shot*(UTn)^sutt*(HEn)^shet*(TRn)^strt*(MIn)^smit
* *** Housing from census, Nonhousing component based on accra regression
* egen M1 = mean(bvalrent^shr_nh1) if year==2001 & accrasamp==1
* egen MM1 = max(M1)
* gen CPI_T1=(var2*bvalrent)^shot*((var7/MM1)*bvalrent^shr_nh1)^(1-shot)
* *** Housing from census, Nonhousing component based on Albouy's number 
* egen M2 = mean(bvalrent^shr_nh2) if year==2001 & accrasamp==1
* egen MM2 = max(M2) 
* gen CPI_T2=(var2*bvalrent)^shot*((var7/MM2)*bvalrent^shr_nh2)^(1-shot)
* keep statefips cntyfips year CPI_T*

*============================================================================
* Deflate to the mean location and 2000 dollars
*============================================================================
run temporalCPI.do
tab year, sum(CPI_Ta) mean // before any deflation
tab year, sum(CPI_Tb) mean // before any deflation
l cbsaname year CPI* gdp gdpPerCapita tempcpi if state==49 & county==49
replace CPI_Ta = CPI_Ta*tempcpi/100 // temporally deflate (2000 is unchanged since it's the base year)
replace CPI_Tb = CPI_Tb*tempcpi/100 // temporally deflate (2000 is unchanged since it's the base year)
tab year, sum(CPI_Ta) mean
tab year, sum(CPI_Tb) mean

* collapse to MSAs (instead of county components)
preserve
	keep if year==2000 & cbsa_code<50000
	collapse pop00 CPI_T*, by(cbsa_code)
	save `cityp'
restore

* keep all rural counties separate
preserve
	keep if year==2000 & cbsa_code>50000 & ~mi(cbsa_code)
	keep pop00 CPI_T* cbsa_code
	save `ruralp'
restore

* generate the mean location in 2000
preserve
	clear
	use `cityp'
	append using `ruralp'
	foreach X of varlist CPI_T* {
		sum `X' [aw=pop00]
		global avg`X' = r(mean)
	}
restore

* spatially deflate by mean location in 2000
foreach X of varlist CPI_T* {
	generat avgloc`X' = ${avg`X'}
	replace `X'       = `X'/avgloc`X' // *100 ?
}

tab year, sum(CPI_Ta) mean
tab year, sum(CPI_Tb) mean
l cbsaname year CPI* gdp gdpPerCapita tempcpi if state==49 & county==49
tab year, sum(gdp) mean
tab year, sum(gdpPerCapita) mean

* MSA GDP and GDP per capita is in 2005 dollars; convert to 2000 dollars
gen  defA = tempcpi/100 if year==2005
egen def  = mean(defA)
mdesc def
replace gdp = gdp / def
replace gdpPerCapita = gdpPerCapita / def

tab year, sum(gdp) mean
l cbsaname year CPI* gdp* tempcpi if state==49 & county==49


*****************************************************************************
* Merge in city unemployment statistics
*****************************************************************************
capture noisily drop _merge
merge m:1 cbsa_code year using `MSAunemplong'

capture noisily drop _merge
capture noisily drop _merge
capture noisily drop _merge
merge 1:1 county state year using `countyunemplong', gen(merge1)
replace labor_force = labor_force_rural if cbsa_code==99999
replace employed    = employed_rural    if cbsa_code==99999
replace unemployed  = unemployed_rural  if cbsa_code==99999
drop labor_force_rural employed_rural unemployed_rural

mdesc labor_force employed unemployed // missing Kalawao County, HI and some rural Alaskan counties (consolidated together in BLS data but not in FIPS data)

* l state county cntyname cbsa_code year labor_force if mi(labor_force)

gen urate = unemployed/labor_force
gen erate = employed/labor_force
gen lfrate00 = labor_force/pop00
gen lfrate10 = labor_force/pop10



*****************************************************************************
* Merge in city unemployment statistics
*****************************************************************************
capture noisily drop _merge
merge m:1 cbsa_code using `MSAlandUseLong'

capture noisily drop _merge
capture noisily drop _merge
capture noisily drop _merge
merge m:1 county state using `countyLandUseLong', gen(merge11)
replace area_cou = area_cou_rural if cbsa_code==99999
replace area_ua  = area_ua_rural  if cbsa_code==99999
drop area_cou_rural area_ua_rural

mdesc area_cou area_ua // missing ...???

* l state county cntyname cbsa_code year labor_force if mi(labor_force)

gen ua_pct10 = area_ua/area_cou



*****************************************************************************
* Merge in employment shares, income per worker, etc.
*****************************************************************************
local BEAvars growthEmp growthManuf growthService growthInd1 growthInd2 growthInd3 growthInd4 growthInd5 growthInd6 growthInd7 growthInd8 growthInd9 ind1IncPerCapita ind2IncPerCapita ind3IncPerCapita ind4IncPerCapita ind5IncPerCapita ind6IncPerCapita ind7IncPerCapita ind8IncPerCapita ind9IncPerCapita ind1Emp ind2Emp ind3Emp ind4Emp ind5Emp ind6Emp ind7Emp ind8Emp ind9Emp ind1Inc ind2Inc ind3Inc ind4Inc ind5Inc ind6Inc ind7Inc ind8Inc ind9Inc empTot ind1Share ind2Share ind3Share ind4Share ind5Share ind6Share ind7Share ind8Share ind9Share farmShare incPerCapita manufIncPerCapita serviceIncPerCapita m_growthEmp m_growthManuf m_growthService m_incPerCapita m_manufIncPerCapita m_serviceIncPerCapita m_growthInd1 m_growthInd2 m_growthInd3 m_growthInd4 m_growthInd5 m_growthInd6 m_growthInd7 m_growthInd8 m_growthInd9 m_ind1IncPerCapita m_ind2IncPerCapita m_ind3IncPerCapita m_ind4IncPerCapita m_ind5IncPerCapita m_ind6IncPerCapita m_ind7IncPerCapita m_ind8IncPerCapita m_ind9IncPerCapita m_ind1Share m_ind2Share m_ind3Share m_ind4Share m_ind5Share m_ind6Share m_ind7Share m_ind8Share m_ind9Share m_ind1Emp m_ind2Emp m_ind3Emp m_ind4Emp m_ind5Emp m_ind6Emp m_ind7Emp m_ind8Emp m_ind9Emp m_empTot 
preserve
	use ${pathstem}${source}county_data_Tyler.dta, clear
	drop if year<1990
	gen state = fips_st
	gen county = fips_co
	drop if state==0|county==0
	drop FIPS
	isid county state year
	save `LLM', replace
restore

* // note: Maui + Kalawao counties were aggregated in the BEA data to form FIPS 15901
capture noisily drop _merge
merge m:1 cbsa_code year using ${pathstem}${source}BEA_MSA.dta, keepusing( `BEAvars')
drop if _merge==2
foreach var in `BEAvars' {
	ren `var' `var'c
}

capture noisily drop _merge
merge 1:1 county state year using `LLM', keepusing( `BEAvars')
drop if _merge==2
foreach var in `BEAvars' {
	replace `var' = `var'c if cbsa_code<50000
	drop `var'c
}
mdesc `BEAvars'

replace cbsaname = cntyname if cbsaname=="99999" & regexm(substr(cntyname,-2,.),"[A-Z][A-Z]$")
replace cbsaname = cntyname+stabb if cbsaname=="99999"

*============================================================================
* clean up and save
*============================================================================
foreach x in pop00 pop10 dens00 dens10 {
	ren `x'c `x'
}
keep state county year cbsa_code cbsaname metarea msaname stabb cntyname pop?? landarea dens?? latitude longitude gdp* rural csize region CPI_Ta CPI_Tb labor_force employed unemployed urate erate lfrate00 ua_pct10 area_ua area_cou `BEAvars'
l state county cbsa_code cbsaname year pop00 employed unemployed labor_force urate erate lfrate?? if cbsa_code==39340 & year<1991
bys cbsa_code year: gen mainObs = _n==1
save ${pathstem}${final}county_city_chars.dta, replace

mdesc CPI_T* if year<2009
mdesc CPI_T* if year>2008

* Graph price index and other locational characteristics by size category

generat size_cat=1 if pop00<194000
replace size_cat=2 if inrange(pop00,194000,1800000)
replace size_cat=3 if pop00>1800000 & ~mi(pop00)

tab size_cat, sum(pop00)
sum pop00 if size_cat==1
sum pop00 if size_cat==2
sum pop00 if size_cat==3

tab size_cat if year==2004, sum(CPI_Ta)
tab size_cat if year==2004, sum(CPI_Tb)

kdensity CPI_Ta if year==2004 & size_cat==1
graph export smallCOLI.eps, replace

kdensity CPI_Ta if year==2004 & size_cat==2
graph export medCOLI.eps, replace

kdensity CPI_Ta if year==2004 & size_cat==3
graph export bigCOLI.eps, replace

twoway (kdensity CPI_Ta if year==2004 & size_cat==1) (kdensity CPI_Ta if year==2004 & size_cat==2) (kdensity CPI_Ta if year==2004 & size_cat==3), legend(label(1 "< 194,000") label(2 "194,000 ~ 1.8 mil") label(3 "> 1.8 mil") cols(3)) xtitle("Cost of living index") ytitle("Density of PDF") name(COLIdist) graphregion(fcolor(gs16) lcolor(gs16)) // note("Data sources:  ACCRA, BLS-LAUS, BEA).")
graph export COLIdist.eps, replace

twoway (kdensity CPI_Tb if year==2004 & size_cat==1) (kdensity CPI_Tb if year==2004 & size_cat==2) (kdensity CPI_Tb if year==2004 & size_cat==3), legend(label(1 "< 194,000") label(2 "194,000 ~ 1.8 mil") label(3 "> 1.8 mil") cols(3)) xtitle("Cost of living index") ytitle("Density of PDF") name(COLIdistb) graphregion(fcolor(gs16) lcolor(gs16)) // note("Data sources:  ACCRA, BLS-LAUS, BEA).")
graph export COLIdistb.eps, replace

gen gdp000 = gdp/1000
gen gdpPerCapita000 = gdpPerCapita/1000

kdensity gdp000 if year==2004 & size_cat==1
graph export smallgdp.eps, replace

twoway (kdensity gdpPerCapita000 if year==2004 & size_cat==1, bwidth(2.2)) (kdensity gdpPerCapita000 if year==2004 & size_cat==2) (kdensity gdpPerCapita000 if year==2004 & size_cat==3), xtitle("Real GDP per capita") ytitle("Density of PDF") name(gdpdist) graphregion(fcolor(gs16) lcolor(gs16)) // legend(label(1 "< 194,000") label(2 "194,000 ~ 1.8 mil") label(3 "> 1.8 mil") cols(1))
* graph save gdpdist, replace
graph export gdpdist.eps, replace

twoway (kdensity incPerCapita if year==2004 & size_cat==1, bwidth(2)) (kdensity incPerCapita if year==2004 & size_cat==2) (kdensity incPerCapita if year==2004 & size_cat==3), xtitle("Real income per capita") ytitle("Density of PDF") name(incPerCapitadist) graphregion(fcolor(gs16) lcolor(gs16)) // legend(label(1 "< 194,000") label(2 "194,000 ~ 1.8 mil") label(3 "> 1.8 mil") cols(1))
* graph save incPerCapitadist, replace
graph export incPerCapitadist.eps, replace

twoway (kdensity urate if year==2004 & size_cat==1, bwidth(.04)) (kdensity urate if year==2004 & size_cat==2) (kdensity urate if year==2004 & size_cat==3), legend(label(1 "< 194,000") label(2 "194,000 ~ 1.8 mil") label(3 "> 1.8 mil") cols(1)) xtitle("Unemployment rate") ytitle("Density of PDF") name(uratedist) graphregion(fcolor(gs16) lcolor(gs16)) // legend(label(1 "< 194,000") label(2 "194,000 ~ 1.8 mil") label(3 "> 1.8 mil") cols(1))
* graph save uratedist, replace
graph export uratedist.eps, replace


* graph combine "COLIdist" "uratedist" "gdpdist" "incPerCapitadist", title("Distributions of city characteristics by city size")
grc1leg COLIdist uratedist gdpdist incPerCapitadist, title("Distributions of city characteristics by city size") graphregion(fcolor(gs16) lcolor(gs16)) note("Data sources:  ACCRA, BLS-LAUS, BEA")
graph export alldist.eps, replace

* Scatter plot of gdp and population to identify high productivity/low population areas as well as low productivity/high population areas
gen lnpop = ln(pop00)
gen lngdp = ln(gdp000)
graph twoway (scatter lngdp lnpop if year==2004, mlabel(stabb))
graph export gdp_pop_outliers.eps, replace 

* List motivating examples of different portions of the GDP distribution
xtile pop_pctile = pop00        if year==2004, nq(100)
xtile gdp_pctile = gdpPerCapita if year==2004, nq(100)
xtile ipc_pctile = incPerCapita if year==2004, nq(100)
egen  pop_rank = rank(pop00       ) if year==2004 & mainObs, field
egen  gdp_rank = rank(gdpPerCapita) if year==2004 & mainObs, field
egen  ipc_rank = rank(incPerCapita) if year==2004 & mainObs, field
l cbsaname pop00 pop_pctile pop_rank gdp gdp_pctile gdp_rank incPerCapita ipc_pctile ipc_rank if year==2004 & mainObs & ( regexm(cbsaname,"Rivers") | regexm(cbsaname, "Casp") | regexm(cbsaname, "New Y") | regexm(cbsaname ,"Palm C"))

l cbsaname pop00 pop_pctile pop_rank gdpPerCapita gdp_pctile gdp_rank incPerCapita ipc_pctile ipc_rank if year==2004 & mainObs & ( regexm(cbsaname,"Rivers") | regexm(cbsaname, "Casp") | regexm(cbsaname, "New Y") | regexm(cbsaname ,"Palm C"))


preserve
	collapse pop00 gdp incPerCapita CPI_Ta CPI_Tb urate size_cat, by(cbsa_code year)
	gen incPerCapita_Ta = incPerCapita/CPI_Ta
	gen incPerCapita_Tb = incPerCapita/CPI_Tb
	gen gdpt            = gdp/1000
	gen gdp_Ta          = gdp/(1000*CPI_Ta)
	gen gdp_Tb          = gdp/(1000*CPI_Tb)
	tab size_cat if year==2004, sum(CPI_Ta)
	tab size_cat if year==2004, sum(CPI_Tb)
	tab size_cat if year==2004, sum(gdpt)
	tab size_cat if year==2004, sum(gdp_Ta)
	tab size_cat if year==2004, sum(gdp_Tb)
	tab size_cat if year==2004, sum(incPerCapita)
	tab size_cat if year==2004, sum(incPerCapita_Ta)
	tab size_cat if year==2004, sum(incPerCapita_Tb)
	tab size_cat if year==2004, sum(urate)
restore

replace gdp    = 0 if mi(gdp)    // so that gdp doesn't get read in as a string when I import it elsewhere
replace gdp    = 12345.67891 in 1
replace CPI_Ta = 0 if mi(CPI_Ta) // so that CPI doesn't get read in as a string when I import it elsewhere
replace CPI_Tb = 0 if mi(CPI_Tb) // so that CPI doesn't get read in as a string when I import it elsewhere

* l cbsaname cntyname latitude longitude if regexm(cbsaname,"New York") | regexm(cbsaname,"Salt Lake") | regexm(cbsaname,"Provo") | regexm(cbsaname,"Durham") | regexm(cntyname,"Beaver UT") | regexm(cntyname,"Iron UT") | regexm(cntyname,"Grand UT") | regexm(cntyname,"Piute UT")

outsheet using ${pathstem}${final}county_city_chars.csv, comma nol replace
!gzip -f ${pathstem}${final}county_city_chars.csv

* test vincenty code:
* generate NYlat  = 40.64
* generate NYlong = -73.94
* vincenty latitude longitude NYlat NYlong, vin(vindist) hav(havdist) loc(locdist) 

* l cntyname latitude longitude vindist havdist locdist if year==2010 & (regexm(cbsaname,"New York") | regexm(cbsaname,"Salt Lake") | regexm(cbsaname,"Provo") | regexm(cbsaname,"Durham") | regexm(cntyname,"Beaver UT") | regexm(cntyname,"Iron UT") | regexm(cntyname,"Grand UT") | regexm(cntyname,"Piute UT") )

log close
exit
