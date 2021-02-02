/*** cpi.do  
This do-file builds the spatial price deflator for each year
in each county ***/


capture log close
log using cpi.log, text replace

clear all
set more off

global pathstem /afs/econ.duke.edu/data/vjh3/TylerJMP/AccraData/
global source   sourcedata/
global final    finaldata/

****** 0. Set up County Populations in 2000 *********

use ${pathstem}${final}county_city_chars
keep state county pop00
sort state county
save temppop00.dta, replace


******* 1. Assign Monthly or Biannual CPI data to Each Location ***********

*** Merge monthly and biannual data
use ${pathstem}${source}cpi_monthly.dta // from Baum-Snow and Pavan (2012)'s online data appendix
replace loc="St Louis" if loc=="St. Louis"
rename tranportation transportation
sort loc year
merge 1:1 loc year using ${pathstem}${source}cpi_biannual // from Baum-Snow and Pavan (2012)'s online data appendix
tab _merge

*** Assign Biannual CPI numbers to Location/Years Without Monthly CPI information
bys loc (year): replace all= all1 if all[_n-1]==all1[_n-1] & all==. & all1~=. & all[_n-1]~=.
bys loc (year): replace food= food1 if food[_n-1]==food1[_n-1] & food==. & food1~=. & food[_n-1]~=.
bys loc (year): replace  shelter=  shelter1 if  shelter[_n-1]== shelter1[_n-1] &  shelter==. &  shelter1~=. &  shelter[_n-1]~=.
bys loc (year): replace  utilities=  utilities1 if  utilities[_n-1]== utilities1[_n-1] &  utilities==. &  utilities1~=. &  utilities[_n-1]~=.
bys loc (year): replace  medical=  medical1 if  medical[_n-1]== medical1[_n-1] &  medical==. &  medical1~=. &  medical[_n-1]~=.
bys loc (year): replace  transportation=  transportation1 if  transportation[_n-1]== transportation1[_n-1] &  transportation==. & transportation1~=. &  transportation[_n-1]~=.
sum y  all all1
ta loc if all==. & all1~=.

*** Tampa and Phoenix Only have Biannual Data
bys loc (year): replace all= all1 if loc=="Tampa" | loc=="Phoenix"
bys loc (year): replace food= food1 if loc=="Tampa" | loc=="Phoenix" 
bys loc (year): replace  shelter=  shelter1 if   loc=="Tampa" | loc=="Phoenix"
bys loc (year): replace  utilities=  utilities1 if   loc=="Tampa" | loc=="Phoenix"
bys loc (year): replace  medical=  medical1 if   loc=="Tampa" | loc=="Phoenix"
bys loc (year): replace  transportation=  transportation1 if   loc=="Tampa" | loc=="Phoenix"

*** For other locations, only medical needs to be fixed
ta loc if food==. & food1~=.
ta loc if shelter==. & shelter1~=.
ta loc if utilities==. & utilities1~=.
ta loc if transportation==. & transportation1~=.
ta loc if medical==. & medical1~=.
replace medical=medical1 if medical==. & medical1~=.

drop all1 food1 shelter1 utilities1 medical1 transportation1
save cpi.dta,replace

clear


******** 2. Ready Expenditure Shares Data ************

use ${pathstem}${source}shares.dta // from Baum-Snow and Pavan (2012)'s online data appendix
*** Drop number of observations
drop n*
drop if year==1980

*** Use Shares Data from 82=82-83 for 83 and 84, etc.
replace year = year+1

*** These are the shares for the unselected sample
drop sgra suta shoa stra shea smia

expand 2
sort year
replace year=year[_n-1]+1 if year==year[_n-1]
expand 6 if year==2002
replace year=year[_n-1]+1 if year==2002
gsort -year
expand 6 if year==1983
replace year=year[_n-1]-1 if year==1983

sort year
save shares_cpi.dta,replace
clear


**** 3. Merge Datasets and Normalize *******

use cpi.dta
drop _m
sort year
merge m:1 year using shares_cpi.dta
tab _m
drop _m
sort year
merge m:1 year using ${pathstem}${source}cpi_hist_shares.dta  // from Baum-Snow and Pavan (2012)'s online data appendix
drop if year==2008
tab _m
drop _m

*** Calculate Local Nonhousing Index
gen nonhous = (100/(100-CPI_shelter))*(all-CPI_shelter*shelter/100)
** Use this sample for BLS regressions below
gen sampbls = 0
#delimit ;
replace sampbls = 1 if
loc=="Anchorage"
|loc=="Atlanta"
|loc=="Boston"
|loc=="Chicago"
|loc=="Cincinnati"
|loc=="Cleveland"
|loc=="Dallas"
|loc=="Denver"
|loc=="Detroit"
|loc=="Honolulu"
|loc=="Houston"
|loc=="Kansas City"
|loc=="Los Angeles"
|loc=="Miami"
|loc=="Milwaukee"
|loc=="Minneapolis"
|loc=="New York"
|loc=="Philly"
|loc=="Phoenix"
|loc=="Pittsburgh"
|loc=="Portland"
|loc=="San Diego"
|loc=="San Francisco"
|loc=="Seattle"
|loc=="St Louis"
|loc=="Tampa"
|loc=="Washington";
#delimit cr

**** Calculate Miscellaneous Index
gen mish = (100/CPI_mish)*(all-(CPI_food*food/100+CPI_shelter*shelter/100+CPI_utilities*utilities/100+CPI_medical*medical/100+ CPI_transportation*transportation/100))
drop CPI*

keep if year>=1978 & year<2008

**** Renormalize all components to be relative to the first year for which we have all obs in each location
sort loc year
drop if medical==. | transportation==. | shelter==. | food==. | utilities==.|mish==.
by loc: gen rrG=food[1]
by loc: gen rrS=shelter[1]
by loc: gen rrU=utilities[1]
by loc: gen rrH=medical[1]
by loc: gen rrT=transportation[1]
by loc: gen rrM=mish[1]
by loc: gen rrN=nonhous[1]
by loc: gen rrA=all[1]
by loc: replace food=food/rrG
by loc: replace shelter=shelter/rrS
by loc: replace utilities=utilities/rrU
by loc: replace medical=medical/rrH
by loc: replace transportation=transportation/rrT
by loc: replace mish=mish/rrM
by loc: replace nonhous=nonhous/rrN
by loc: replace all = all/rrA
drop rr*

sort loc year


**** 4. Fill in Missing Years for Locations/Year Combinations Without Data

by loc:egen start=min(year)
tab start
tab loc if start>1978

drop if loc=="Size B/C" | loc=="Size A" 
* Expand starting in First year of Data for these locations 
expand 21 if loc=="Midwest B/C" & year==1998
expand 21 if loc=="NorthEast B/C" & year==1998
expand 21 if loc=="West B/C" & year==1998
expand 21 if loc=="South B/C" & year==1998
expand 21 if loc=="Washington" & year==1998
expand 25 if loc=="Phoenix" & year==2002
expand 10 if loc=="Tampa" & year==1987
gsort loc -year
by loc: replace sampbls = 0 if year>=year[_n-1]
by loc: replace year=year[_n-1]-1 if year>=year[_n-1]

rename food var1
rename shelter var2
rename utilities var3
rename transportation var4
rename medical var5
rename mish var6
rename nonhous var7
rename all var8

*** For years without data, replace with value of 10
local i=1
while `i'<=8 {
replace var`i'=10 if var`i'==1 & (var`i'[_n-1]==1 | var`i'[_n-1]==10)
local i=`i'+1
}
gen temp=0

*Assign Missing years of Phoenix to West A....
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="West A"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
replace var`i'=temp if loc=="Phoenix"
local i=`i'+1
}

* Midwest
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="Midwest U"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
gsort loc -year
by loc:replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if var`i'==10 & loc=="Midwest B/C"
local i=`i'+1
}

* NorthEast
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="NorthEast U"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
gsort loc -year
by loc:replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if var`i'==10 & loc=="NorthEast B/C"
local i=`i'+1
}

* West
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="West U"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
gsort loc -year
by loc:replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if var`i'==10 & loc=="West B/C"
local i=`i'+1
}

* South
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="South U"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
gsort loc -year
by loc:replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if var`i'==10 & loc=="South B/C"
local i=`i'+1
}

* Washington and Tampa
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="South A"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
gsort loc -year
by loc:replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if var`i'==10 & (loc=="Washington" | loc=="Tampa")
local i=`i'+1
}

*** Fix up other cases when needing to fill forward
by loc:egen end=max(year)
tab end
tab loc if end<2005
by loc:gen n=_n
by loc:egen nn=max(n)
tab nn loc if nn<28

* Cleveland and St louis are missing 2003. Just take the average of 2002 and 2003.
expand 2 if year==2002 & (loc=="Cleveland" | loc=="St Louis") 
sort loc year
replace year=2003 if year==2002 & year[_n-1]==2002

local i=1
while `i'<=8 {
replace var`i'=(var`i'[_n-1]+var`i'[_n+1])/2 if year==2003 & (loc=="Cleveland" | loc=="St Louis")
local i=`i'+1
}
replace sampbls = 0 if year==2003 & (loc=="Cleveland" | loc=="St Louis")

drop if year>2001 & (loc=="Anchorage" | loc=="Honolulu")

expand 4 if (loc=="Minneapolis") & year==2004
expand 8 if year==2000 & loc=="Honolulu"
expand 7 if year==2001 & loc=="Anchorage" 

sort loc year
by loc:replace year=year[_n-1]+1 if year>=2000

*Honolulu
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="West B/C"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
sort loc year
replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if year>=2001 & loc=="Honolulu"
local i=`i'+1
}
replace sampbls = 0 if year>=2001 & loc=="Honolulu"

*Minneapolis
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="Midwest A"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
gsort loc -year
replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if year==2005 & loc=="Minneapolis"
local i=`i'+1
}
replace sampbls = 0 if year==2005 & loc=="Minneapolis"

*Anchorage
local i=1
while `i'<=8 {
drop temp
gen temp=var`i' if loc=="West B/C"
sort year temp
by year:replace temp=temp[_n-1] if temp==.
sort loc year
replace var`i'=var`i'[_n-1]*temp/temp[_n-1] if year>=2002 & loc=="Anchorage"
local i=`i'+1
}
replace sampbls = 0 if year>=2002 & loc=="Anchorage"

drop temp end start nn n


******* 5. Assign Area Codes for Merging to Accra Data **********

gen AREA=0
replace AREA=1 if loc=="Anchorage"
replace AREA=2 if loc=="Atlanta"
replace AREA=3 if loc=="Boston"
replace AREA=4 if loc=="Chicago"
replace AREA=5 if loc=="Cincinnati"
replace AREA=6 if loc=="Cleveland"
replace AREA=7 if loc=="Dallas"
replace AREA=8 if loc=="Denver"
replace AREA=9 if loc=="Detroit"
replace AREA=10 if loc=="Honolulu"
replace AREA=11 if loc=="Houston"
replace AREA=12 if loc=="Kansas City"
replace AREA=13 if loc=="Los Angeles"
replace AREA=14 if loc=="Miami"
replace AREA=15 if loc=="Midwest A"
replace AREA=16 if loc=="Midwest B/C"
replace AREA=17 if loc=="Midwest D"
replace AREA=18 if loc=="Milwaukee"
replace AREA=19 if loc=="Minneapolis"
replace AREA=20 if loc=="New York"
replace AREA=21 if loc=="NorthEast A"
replace AREA=22 if loc=="NorthEast B/C"
replace AREA=23 if loc=="Philly"
replace AREA=24 if loc=="Phoenix"
replace AREA=25 if loc=="Pittsburgh"
replace AREA=26 if loc=="Portland"
replace AREA=27 if loc=="San Diego"
replace AREA=28 if loc=="San Francisco"
replace AREA=29 if loc=="Seattle"
replace AREA=30 if loc=="Size D"
replace AREA=31 if loc=="South A"
replace AREA=32 if loc=="South B/C"
replace AREA=33 if loc=="South D"
replace AREA=34 if loc=="St Louis"
replace AREA=35 if loc=="Tampa"
replace AREA=36 if loc=="US"
replace AREA=37 if loc=="Washington"
replace AREA=38 if loc=="West A"
replace AREA=39 if loc=="West B/C"

*** Normalize each component such that 2001=1       
drop if AREA==0
sort AREA year
local i=1
while `i'<=8 {
by AREA: gen start=var`i'[24]
replace var`i'=var`i'/start
drop start
local i=`i'+1
}

** Populate each AREA with US numbers
gen US = (AREA==36)
sort year US
local i = 1
while `i'<=8 {
by year: gen varx`i' = var`i'[_N]
local i = `i'+1
}
drop if AREA==36

sort AREA year
save deflator_1978,replace

*** This is the adjustment factor for the census housing deflator
by AREA: gen Hadjust = var2[24]/var2[23]
by AREA: keep if _n==1
keep AREA Hadjust
sort AREA
save adjust.dta, replace

exit
***************** 6. Build Final Index ******************

*** This data set defines the accra sample which should have a mean of 100 in the index
clear
use ${pathstem}${final}county_city_chars.dta
keep cbsa_code state county accra2001 grocery2001 housing2001 utilities2001 transport2001 health2001 misc2001
foreach var in accra grocery housing utilities transport health misc {
	rename `var'2001 `var'
}
sum housing
keep cbsa_code
sort cbsa_code
bys cbsa_code: keep if _n==1
gen mrg = 1
sort cbsa_code mrg
save tempaccra.dta, replace

clear
set mem 50m
use ${pathstem}${final}accra-2001,clear

*** Put Census Housing Variables in 2001$ and normalize the same as the ACCRA data
sort AREA
merge m:1 AREA using adjust
drop _merge
replace brent = brent*Hadjust
replace bval = bval*Hadjust
replace bvalrent = bvalrent*Hadjust
drop Hadjust

*** Rescale to be relative to the average location surveyed by ACCRA
sort cbsa_code
by cbsa_code: gen mrg = _n
sort cbsa_code mrg
merge 1:1 cbsa_code mrg using tempaccra.dta
** This is the stamford-norwalk MSA that got consolidated into New Haven
l cbsa_code if _merge==2
drop if _merge==2
gen accrasamp = (_merge==3)

/*** Make the average ACCRA MSA in our sample the base location -- this
implies very small adjustments for the ACCRA components ***/
foreach X of varlist bval brent bvalrent accra grocery housing utilities transport health misc {
  egen m`X' = mean(`X') if _merge==3
  egen M`X' = mean(m`X')
  replace `X' = (`X'/M`X')
  drop m`X' M`X'
}
drop _merge

gen n=_n
expand 30

/*we generate time that shows the number of years elapsed since the beginning of year 0. Using this 
we also generate year*/

gen year=1978
sort n
by n: replace year=year[_n-1]+1 if year[_n-1]~=.

sort AREA year

merge AREA year using deflator_1978.dta
l loc year if _m==2
drop if _m==2
drop _m
sort AREA year statefips cntyfips
gen sampbls2 = 0
by AREA year: replace sampbls2 = sampbls if _n==1
drop sampbls
rename sampbls2 sampbls

*** Adjust ACCRA Components
gen GR=grocery*var1
gen HO=housing*var2 
gen UT=utilities*var3
gen TR=transport*var4
gen HE=health*var5
gen MI=misc*var6
gen NHO=((accra-.29*housing)/.71)*var7
gen TOT=accra*var8

*** Components that are the same across locations
gen GRn=varx1
gen UTn=varx3
gen TRn=varx4
gen HEn=varx5
gen MIn=varx6
gen NHOn=varx7

*** Take non-housing component based on regression from ACCRA (Albouy)
gen lNHO = log(NHO)
gen lHO = log(HO)
reg lNHO lHO if accrasamp==1 & year==2001
gen shr_nh1 = _b[lHO]
gen xshr_nh1 = 1-_b[lHO]
*** This is what Albouy finds
gen shr_nh2 = .26
gen xshr_nh2 = .74

*** Generate Weinstein/Handbury Index for groceries
** Assume that the average ACCRA place from this regression is average for the index 
gen lpop80 = log(popm80)
reg grocery lpop80 if accrasamp==1 & year==2001
gen mlpop = (1-_b[_cons])/_b[lpop80] if year==2001
gen dlpop = lpop80-mlpop
gen gr2 = 1-0.0283*dlpop if year==2001
sum gr2
** Deal with 0 population counties
replace gr2 = r(max) if year==2001 & gr2==.
sort statefips cntyfips gr2
by statefips cntyfips: replace gr2 = gr2[1]
gen GR2 = gr2*var1

*** ACCRA Data
gen CPI_Ta=(GR)^sgrt*(HO)^shot*(UT)^sutt*(HE)^shet*(TR)^strt*(MI)^smit
*** Accra Data but with housing component from the census
gen CPI_Tavr=(GR)^sgrt*(var2*bvalrent)^shot*(UT)^sutt*(HE)^shet*(TR)^strt*(MI)^smit
*** Only Spatial differentiation is from the housing component from the census
gen CPI_Tvr=(GRn)^sgrt*(var2*bvalrent)^shot*(UTn)^sutt*(HEn)^shet*(TRn)^strt*(MIn)^smit
*** Only spatial diff from housing component from the census & grocery from Weinstein/Handbury
gen CPI_Tgvr=(GR2)^sgrt*(var2*bvalrent)^shot*(UTn)^sutt*(HEn)^shet*(TRn)^strt*(MIn)^smit
*** Housing from census, Nonhousing component based on accra regression
egen M1 = mean(bvalrent^shr_nh1) if year==2001 & accrasamp==1
egen MM1 = max(M1)
gen CPI_T1=(var2*bvalrent)^shot*((var7/MM1)*bvalrent^shr_nh1)^(1-shot)
*** Housing from census, Nonhousing component based on Albouy's number 
egen M2 = mean(bvalrent^shr_nh2) if year==2001 & accrasamp==1
egen MM2 = max(M2) 
gen CPI_T2=(var2*bvalrent)^shot*((var7/MM2)*bvalrent^shr_nh2)^(1-shot)
keep statefips cntyfips year CPI_T*

*** Normalize Indices to be for the average place in 1999
sort statefips cntyfips
merge statefips cntyfips using temppop00.dta
tab _merge
** _merge=1 is Yellowstone NP in MT with no people
keep if _merge ==3
drop _merge

foreach X of varlist CPI_T* {
sum `X' [aw=pop00] if year==1999
replace `X' = 100*`X'/r(mean)
}
drop pop00

sort statefips cntyfips year

save ../finaldata/price-index,replace

log c

erase deflator_1978.dta
erase cpi.dta
erase shares_cpi.dta
erase temppop00.dta
erase tempaccra.dta
erase adjust.dta
