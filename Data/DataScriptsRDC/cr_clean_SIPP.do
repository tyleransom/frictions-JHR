clear all
version 12.0
capture log close
set more off
log using cr_clean_SIPP.log, replace


! gunzip -f [REDACTED]sipp2004NHWmale.csv.gz
insheet using "[REDACTED]sipp2004NHWmale.csv", comma case
! gzip   -f [REDACTED]sipp2004NHWmale.csv

*----------------------------------------------------
* generate ID variable and declare panel data
*----------------------------------------------------
***replace lgtkey = "[REDACTED]"            if lgtkey == "[REDACTED]"
replace puid   = "[REDACTED]" if puid   == "[REDACTED]"
count if gdp==12345.67891
count if gdpPerCapita==12345.67891
recode gdp gdpPerCapita (12345.67891 = 0) // to fix the gdp recode from Stata outside the  RDC
egen ID = group(puid)

capture noisily drop lngdp gdp000 gdpPerCapita000 gdp_pctile gdp_rank
// replace gdp = gdp*1000000
// Check that mean of gdp is on the order of 1e11:
qui sum gdp
if ~inrange(`r(mean)'/1e11,.1,9) {
disp "error in gdp read-in ... gdp levels are off!"
exit
}

xtset ID timmonth

local cityvars CPI_Ta CPI_Tb FIPScounty FIPScounty_new FIPSstate FIPSstate_new cbsa_code cbsaname cntyname county csize dens00 dens10 gdp gdpPerCapita urate landarea pop00 pop10 region rural stabb state m_* growthInd? ind?IncPerCapita ind?Emp ind?Inc ind?Share farmShare incPerCapita longitude latitude

local cityvars_use CPI_Ta CPI_Tb cbsa_code county csize dens00 dens10 gdp gdpPerCapita urate landarea pop00 pop10 region rural state longitude latitude

local dervars brthyr_final brthmn_final quarter* wqc* flag_in_der flag_in_ser totearn* total_der_* working_th* pos_der_fica_* pos_der_nonfica_* 
destring `dervars', ignore("Z") replace

local qvars qjobflag qrace qhgc qged qenrwhich qhrswork qind1 qocc1 qind2 qocc2 qenrlsch qrsn_stop_j1 qrsn_stop_j2 qrsn_no_work qrsn_absent qon_layoff qmonthearnj1 qmonthearnj2 qmultlocj1 qmultlocj2 qfirmsizej1 qfirmsizej2 qstillatj1 qstillatj2 qsdatej1 qsdatej2 qedatej1 qedatej2 qclassj1 qclassj2 qyrenterlf qyrprevlf qyrexitlf qyrstartpvjob qmosoutlf qexitlfcareg* qhourly1 hourly2 qhourly2 look_work qlook_work qyr_move_this_house qyr_move_pvs_house qyr_move_this_state qcollege_major qbirth_state_country qworkedpt1plusweeks qreasonworkpt

local vars LGT* EORIGIN AORIGIN EBORNUS ABORNUS sex wavemap race year calmonth intstatus hgc_cat ged enrsch enrwhich hrswork rhrsweek ind1 occ1 ind2 occ2 job1ID job2ID  numjobs county state  wavemap `cityvars_use' `dervars' hhtype timmonth hhsize age marst enrlsch empstatwk1 empstatwk2 empstatwk3 empstatwk4 empstatwk5 weekswjob weeksabjob rsn_stop_j1 rsn_stop_j2 rsn_no_work rsn_absent on_layoff weight monthearnj1 monthearnj2 multlocj1 multlocj2 firmsizej1 firmsizej2 stillatj1 stillatj2 sdatej1 sdatej2 edatej1 edatej2 classj1 classj2 yrenterlf yrprevlf yrexitlf yrstartpvjob mosoutlf exitlfcaregiver mover hourly1  exper0 tenure0 ID yr_move_this_house yr_move_pvs_house yr_move_this_state college_major birth_state_country prop1 prop2 workedpt1plusweeks reasonworkpt contingent `qvars'
d `vars', f
recode `vars' (-1 = .v)
recode `qvars' (2 3 4 = 1)
recode yr_move_this_house yr_move_pvs_house yr_move_this_state (-5 -3 = .i) (9999 = .i)

*----------------------------
* generate population categories
*----------------------------
// [REDACTED]
gen mispopA = mi(pop00) | county==0
bys ID (timmonth): egen mispopB = mean(mispopA)
gen misPopFlag = mispopB>0

gen pop_cat = .
replace pop_cat = 1 if inrange(pop00,1      ,193500 )
replace pop_cat = 2 if inrange(pop00,193500 ,1797000)
replace pop_cat = 3 if inrange(pop00,1797000,2e+07  )
bys ID (timmonth ): gen pop_cat_lag = pop_cat[_n-1]
bys ID (timmonth ): gen longitude_lag = longitude[_n-1]
bys ID (timmonth ): gen latitude_lag = latitude[_n-1]
vincenty latitude longitude latitude_lag longitude_lag, hav(move_dist)
bys ID (timmonth ): replace move_dist = 0 if _n==1
qui tab pop_cat, gen(pop)
qui tab pop_cat_lag, gen(pop_lag)

replace pop_cat = 3 if [REDACTED]

lab def vlbirthplace 1 "Alabama" 2 "Alaska" 4 "Arizona" 5 "Arkansas" 6 "California" 8 "Colorado" 9 "Connecticut" 10 "Delaware" 11 "DC" 12 "Florida" 13 "Georgia" 15 "Hawaii" 16 "Idaho" 17 "Illinois" 18 "Indiana" 19 "Iowa" 20 "Kansas" 21 "Kentucky" 22 "Louisiana" 23 "Maine" 24 "Maryland" 25 "Massachusetts" 26 "Michigan" 27 "Minnesota" 28 "Mississippi" 29 "Missouri" 30 "Montana" 31 "Nebraska" 32 "Nevada" 33 "New Hampshire" 34 "New Jersey" 35 "New Mexico" 36 "New York" 37 "North Carolina" 38 "North Dakota" 39 "Ohio" 40 "Oklahoma" 41 "Oregon" 42 "Pennsylvania" 44 "Rhode Island" 45 "South Carolina" 46 "South Dakota" 47 "Tennessee" 48 "Texas" 49 "Utah" 50 "Vermont" 51 "Virginia" 53 "Washington" 54 "West Virginia" 55 "Wisconsin" 56 "Wyoming" 72 "Puerto Rico" 78 "US Virgin Islands/American Samoa/Guam" 106 "Denmark" 109 "France" 110 "Germany" 117 "Hungary" 119 "Ireland/Eire" 120 "Italy" 126 "Holland" 126 "Netherlands" 127 "Norway" 128 "Poland" 130 "Azores" 137 "Switzerland" 139 "England" 140 "Scotland" 148 "Europe" 156 "Slovakia/Slovak Republic" 183 "Latvia" 192 "Russia" 200 "Afghanistan" 205 "Burma" 206 "Cambodia" 207 "China" 209 "Hong Kong" 210 "India" 211 "Indonesia" 212 "Iran" 214 "Israel" 215 "Japan" 217 "Korea/South Korea" 224 "Malaysia" 229 "Pakistan" 231 "Philippines" 237 "Syria" 238 "Taiwan" 239 "Thailand" 240 "Turkey" 242 "Vietnam" 245 "Asia" 252 "Middle East" 253 "Palestine" 300 "Bermuda" 301 "Canada" 310 "Belize" 312 "El Salvador" 313 "Guatemala" 315 "Mexico" 316 "Nicaragua" 317 "Panama" 337 "Cuba" 338 "Dominica" 339 "Dominican Republic" 340 "Grenada" 342 "Haiti" 343 "Jamaica" 351 "Trinidad and Tobago" 353 "Caribbean" 376 "Bolivia" 377 "Brazil" 379 "Colombia" 380 "Ecuador" 383 "Guyana" 389 "South America" 415 "Egypt" 417 "Ethiopia" 421 "Ghana" 427 "Kenya" 436 "Morocco" 440 "Nigeria" 449 "South Africa" 462 "Other Africa" 501 "Australia" 555 "Elsewhere"

lab val birth_state_country vlbirthplace

lab def vlcollmaj 1 "Agriculture/Forestry"  2 "Art/Architecture" 3 "Business/Management" 4 "Communications" 5 "Computer and Information Sciences" 6 "Education" 7 "Engineering" 8 "English/Literature" 9 "Foreign Languages" 10 "Health Sciences" 11 "Liberal Arts/Humanities" 12 "Math/Statistics" 13 "Natural Sciences" 14 "Philosophy/Religion/Theology" 15 "Pre-Professional" 16 "Psychology" 17 "Social Sciences/History" 18 "Other"

lab val college_major vlcollmaj

lab var enrsch "Was ... enrolled in school in this month?"
recode enrsch (2 .v = 0)
lab var enrwhich "At what level or grade was...enrolled?"
lab def vlenrwhich 1 "Elementary grades 1-8" 2 "High school grades 9-12" 3 "College year 1 (freshman)" 4 "College year 2 (sophomore)" 5 "College year 3 (junior)" 6 "College year 4 (senior)" 7 "First year graduate or professional school" 8 "Second year or higher in graduate or professional school" 9 "Vocational, technical, or business school beyond high school level" 10 "Enrolled in college, but not working towards a degree"
lab val enrwhich vlenrwhich

save [REDACTED]intermediate2004, replace
! gzip   -f [REDACTED]intermediate2004.dta



/* 
! gunzip -f [REDACTED]intermediate2004.dta.gz
use [REDACTED]intermediate2004, clear
! gzip   -f [REDACTED]intermediate2004.dta
*/

*******************************************************
*********** GET RID OF PEOPLE IN SCHOOL ***************
*******************************************************
*------------------------------------------------------------------
* fix so that only keep people once they complete highest schooling
*------------------------------------------------------------------
// generate variable for last time enrolled in school
gen last_enr_schA = timmonth*enrsch
bys ID (timmonth): egen last_enr_sch  = max(last_enr_schA)
bys ID (timmonth): egen ever_enr_schA = mean(last_enr_schA)
// generate schoolFlag variable to mark that a person was in school
gen schoolFlagOld = timmonth<=last_enr_sch
gen schoolFlag    = ever_enr_schA>0


*******************************************************************************
*********** FLAG PEOPLE WITH GAPS AND/OR MISSING CALENDAR DATES ***************
*******************************************************************************
*--------------------------------------------------------------------------------------------
* generate calendar month from date of reference month (interview date)
*--------------------------------------------------------------------------------------------
gen yrmo = ym(year,calmonth ) if inlist(timmonth,4,8,12,16,20,24,28,32,36,40,44,48)
bys ID (timmonth): ipolate yrmo timmonth , gen(caldate) epolate
bys ID (timmonth): replace caldate = caldate[_n+1]-1 if timmonth[_n+1]-timmonth[_n]==1 & mi(caldate[_n])
bys ID (timmonth): replace caldate = caldate[_n+1]-1 if timmonth[_n+1]-timmonth[_n]==1 & mi(caldate[_n])
bys ID (timmonth): replace caldate = caldate[_n+1]-1 if timmonth[_n+1]-timmonth[_n]==1 & mi(caldate[_n])
bys ID (timmonth): replace caldate = caldate[_n+1]-1 if timmonth[_n+1]-timmonth[_n]==1 & mi(caldate[_n])
bys ID (timmonth): replace caldate = caldate[_n+1]-1 if timmonth[_n+1]-timmonth[_n]==1 & mi(caldate[_n])
bys ID (timmonth): replace caldate = caldate[_n-1]+1 if timmonth[_n]-timmonth[_n-1]==1 & mi(caldate[_n])
format caldate %tm
count if mi(caldate)
tab ID if mi(caldate)
gen noCalDateFlag = inlist(ID,[REDACTED]) // [REDACTED]
count if mi(caldate)
capture noisily assert ~mi(caldate)
replace year = year(dofm(caldate)) if year==0
replace calmonth = month(dofm(caldate))


*******************************************************
*********** FIX BAD AGES ******************************
*******************************************************
*------------------------------------------------------------------
* fix so that only keep people once they complete highest schooling
*------------------------------------------------------------------
gen birthYr = brthyr_final
gen birthMo = brthmn_final
gen DOByrmo = ym(birthYr,birthMo)
format DOByrmo %tm
gen age_true_mo = (caldate-DOByrmo)/12
mdesc age_true_mo age
sum flag_in_der
replace age = age_true_mo


*********************************************************************************
*********** GET RID OF PEOPLE ABSENT IN WAVE 1 (WHEN EMP INFO COLLECTED) ********
*********************************************************************************
bys ID (timmonth): gen  firstObsTemp = _n==1
bys ID (timmonth): gen  in_w1 = timmonth[1]==1
sum in_w1 if firstObsTemp


****************************************************************************
*********** GET RID OF PEOPLE OUT OF AGE RANGE *****************************
****************************************************************************
*------------------------------------------------------------------
* only keep people in prime working ages
*------------------------------------------------------------------
bys ID (timmonth): gen  agebeg = age[1] // generate age at beginning of survey
gen young = agebeg<18
gen old   = agebeg>55
bys ID (timmonth): gen ageflag = ((age[_n+1]-age[_n]>1) | (age[_n+1]-age[_n]<-1)) & timmonth[_n+1]==timmonth[_n]+1 & _n>1 & ID[_n+1]==ID[_n]
bys ID (timmonth ):egen ageFlagA=mean(ageflag )
gen ageFlag = ageFlagA>0 //[REDACTED]
sum ageFlag // [REDACTED]


*--------------------------------------------------------------------------------------------------------------
* generate gapped history indicator and drop observations after the gap (approx 118k obs, or 10% of the sample)
*--------------------------------------------------------------------------------------------------------------
bys ID (timmonth): gen gap = caldate[_n]-caldate[_n-1]~=1 if _n>1
gen t_gap = timmonth*gap
recode t_gap (0 = .)
bys ID (timmonth): egen first_gap = min(t_gap)
gen post_gap = (timmonth>=first_gap)
drop gap t_gap first_gap
gen gapFlag = post_gap


*--------------------------
* fix wage variable(s)
*--------------------------
generat WageHr      = (monthearnj1+monthearnj2)/( weekswjob            *hrswork) if weeksabjob>=.
replace WageHr      = (monthearnj1+monthearnj2)/((weekswjob-weeksabjob)*hrswork) if weeksabjob<.
replace WageHr      = WageHr/1.068524971 if year==2003 // adjust for inflation
replace WageHr      = WageHr/1.096980256 if year==2004 // adjust for inflation
replace WageHr      = WageHr/1.134146341 if year==2005 // adjust for inflation
replace WageHr      = WageHr/1.170731707 if year==2006 // adjust for inflation
replace WageHr      = WageHr/1.204076655 if year==2007 // adjust for inflation
replace WageHr      = WageHr/1.250307782 if year==2008 // adjust for inflation
generat lnWageHr    = log(WageHr)
generat qlnWageHr   = qmonthearnj1 | qmonthearnj2

generat WageHrJ     = (monthearnj1+monthearnj2)/( weekswjob            *hrswork) if weeksabjob>=.
replace WageHrJ     = (monthearnj1+monthearnj2)/((weekswjob-weeksabjob)*hrswork) if weeksabjob<.
replace WageHrJ     = WageHrJ/CPI_Ta // already temporally deflated
generat lnWageHrJ   = log(WageHrJ)

generat WageHrJb    = (monthearnj1+monthearnj2)/( weekswjob            *hrswork) if weeksabjob>=.
replace WageHrJb    = (monthearnj1+monthearnj2)/((weekswjob-weeksabjob)*hrswork) if weeksabjob<.
replace WageHrJb    = WageHrJ/CPI_Tb // already temporally deflated
generat lnWageHrJb  = log(WageHrJb)

generat earn        = monthearnj1+monthearnj2
replace earn        = earn/1.068524971 if year==2003 // adjust for inflation
replace earn        = earn/1.096980256 if year==2004 // adjust for inflation
replace earn        = earn/1.134146341 if year==2005 // adjust for inflation
replace earn        = earn/1.170731707 if year==2006 // adjust for inflation
replace earn        = earn/1.204076655 if year==2007 // adjust for inflation
replace earn        = earn/1.250307782 if year==2008 // adjust for inflation
generat lnearn      = ln(earn) // generate log earnings just to see if my hourly measure is causing problems

generat derearn0    = .
replace derearn0    = total_der_fica_2003+total_der_nonfica_2003 if year==2003
replace derearn0    = total_der_fica_2004+total_der_nonfica_2004 if year==2004
replace derearn0    = total_der_fica_2005+total_der_nonfica_2005 if year==2005
replace derearn0    = total_der_fica_2006+total_der_nonfica_2006 if year==2006
replace derearn0    = total_der_fica_2007+total_der_nonfica_2007 if year==2007
replace derearn0    = total_der_fica_2008+total_der_nonfica_2008 if year==2008
generat lnderearn0  = ln(derearn0) // generate log earnings just to see if my hourly measure is causing problems

generat derearn     = .
replace derearn     = (total_der_fica_2003+total_der_nonfica_2003)/1.068524971 if year==2003 // adjust for inflation
replace derearn     = (total_der_fica_2004+total_der_nonfica_2004)/1.096980256 if year==2004 // adjust for inflation
replace derearn     = (total_der_fica_2005+total_der_nonfica_2005)/1.134146341 if year==2005 // adjust for inflation
replace derearn     = (total_der_fica_2006+total_der_nonfica_2006)/1.170731707 if year==2006 // adjust for inflation
replace derearn     = (total_der_fica_2007+total_der_nonfica_2007)/1.204076655 if year==2007 // adjust for inflation
replace derearn     = (total_der_fica_2008+total_der_nonfica_2008)/1.250307782 if year==2008 // adjust for inflation
generat lnderearn   = ln(derearn) // generate log earnings just to see if my hourly measure is causing problems

generat earnJ       = monthearnj1+monthearnj2
replace earnJ       = earnJ/CPI_Ta // already temporally deflated
generat lnearnJ     = ln(earnJ)

generat derearnJ    = derearn0
replace derearnJ    = derearnJ/CPI_Ta // already temporally deflated
generat lnderearnJ  = ln(derearnJ)

generat earnJb      = monthearnj1+monthearnj2
replace earnJb      = earnJ/CPI_Tb // already temporally deflated
generat lnearnJb    = ln(earnJb)

generat derearnJb   = derearn0
replace derearnJb   = derearnJb/CPI_Tb // already temporally deflated
generat lnderearnJb = ln(derearnJb)


*--------------------------
* fix education variable
*--------------------------
generat hgc = .
replace hgc = 0    if hgc_cat==31
replace hgc = 2.5  if hgc_cat==32
replace hgc = 5.5  if hgc_cat==33
replace hgc = 7.5  if hgc_cat==34
replace hgc = 9    if hgc_cat==35
replace hgc = 10   if hgc_cat==36
replace hgc = 11   if hgc_cat==37
replace hgc = 11.5 if hgc_cat==38
replace hgc = 12   if hgc_cat==39
replace hgc = 13   if hgc_cat==40
replace hgc = 13   if hgc_cat==41
replace hgc = 14   if hgc_cat==43
replace hgc = 16   if hgc_cat==44
replace hgc = 18   if hgc_cat==45
replace hgc = 19   if hgc_cat==46
replace hgc = 20   if hgc_cat==47

* flag increasing hgc out of school
bys ID (timmonth): gen ihgcFlag = hgc[_n+1]~=hgc[_n] & enrsch[_n+1]==0 & enrsch[_n]==0 & _n>1 & ID[_n+1]==ID[_n]
bys ID (timmonth): egen everihgcFlag = mean(ihgcFlag)
replace ihgcFlag = everihgcFlag>0

* impute highest hgc if hgc increased out of schoolFlag
gen hgcp = hgc
bys ID (timmonth): egen hgcimp = max(hgc) if everihgcFlag>0
replace hgcp = hgcimp if everihgcFlag>0 & timmonth>last_enr_sch
replace hgc = hgcp
drop hgcp hgcimp ihgcFlag

gen HSgrad  = inrange(hgc,12,15.99)
gen colgrad = inrange(hgc,16,.)

generat educlevel = .
replace educlevel = 1 if ~HSgrad & ~colgrad //missing hgc is OK because it's people under age 15
replace educlevel = 2 if HSgrad
replace educlevel = 3 if colgrad

lab def vleduc 1 "Didn't graduate HS" 2 "HS graduate" 3 "BA graduate"
lab val educlevel vleduc

tab educlevel, mi


*----------------------------------------------------
* generate labor force participation/outcomes
*----------------------------------------------------
gen selfEmp =  inlist(prop1,1,2) | inlist(prop2,1,2) | numjobs==0 | qjobflag==1 | contingent==1
gen wantPT  =  reasonworkpt==2
gen inlf    = (inrange(rhrsweek,1,6) | (rhrsweek==0 & look_work==1)) & ~selfEmp & ~wantPT
gen empFT   =  inrange(rhrsweek,1,1) & inlf & inlist(reasonworkpt,3,4,10,.v) 
gen emp     =  inrange(rhrsweek,1,2) & inlf & inlist(reasonworkpt,3,4,10,.v) 
gen unemp   =  inlf & ~empFT

gen empFTold = inrange(rhrsweek,1,1) & ~inlist(prop1,1,2) & ~inlist(prop2,1,2) & inrange(numjobs,1,.) & qjobflag~=1 & contingent~=1
gen empold   = inrange(rhrsweek,1,2) & ~inlist(prop1,1,2) & ~inlist(prop2,1,2) & inrange(numjobs,1,.) & qjobflag~=1 & contingent~=1



*----------------------------------------------------------------------------
* get months full-time and part-time employed (to get "monthly" W-2 earnings)
*----------------------------------------------------------------------------
bys ID year (timmonth): egen num_mos_FT      = total(empFTold==1)
bys ID year (timmonth): egen num_mos_PT      = total(empold==1 & empFTold==0)
bys ID year (timmonth): egen num_mos_work    = total(empold==1)
bys ID year (timmonth): egen num_mos_in_data = count(empold==1)
gen tot_month_earn_SIPP = monthearnj1+monthearnj2

gen FTrate     = num_mos_FT/num_mos_work
gen FTrateData = num_mos_FT/num_mos_in_data

gen nmflag1 = num_mos_FT==num_mos_in_data & num_mos_in_data>=5
gen nmflag2 = num_mos_PT==num_mos_in_data & num_mos_in_data>=5
gen nmflag3 = num_mos_FT>0 & num_mos_PT>0 & num_mos_in_data>=5
gen nmflag4 = num_mos_FT==num_mos_work    & num_mos_in_data>=5
gen nmflag5 = num_mos_PT==num_mos_work    & num_mos_in_data>=5

gen derearn_monthly = 0
* Case 1: all months of the year working FT
replace derearn_monthly = derearn/12         if nmflag1 & empFTold
* Case 2: all months of the year working PT
replace derearn_monthly = derearn/12         if nmflag2 & empold & ~empFTold
* Case 3: some months in FT, some months in PT
replace derearn_monthly =      derearn/(12*(2.25*(FTrate)+(1-FTrate)))         if nmflag3 & empold & ~empFTold
replace derearn_monthly = 2.25*derearn/(12*(2.25*(FTrate)+(1-FTrate)))         if nmflag3          &  empFTold
* Case 4: all months of the year working either FT or PT
replace derearn_monthly = derearn/num_mos_work if (nmflag4 | nmflag5) & ~nmflag1 & ~nmflag2 & empold
generat lnderearn_monthly  = ln(derearn_monthly)

gen derearn0_monthly = 0
* Case 1: all months of the year working FT
replace derearn0_monthly = derearn0/12         if nmflag1 & empFTold
* Case 2: all months of the year working PT
replace derearn0_monthly = derearn0/12         if nmflag2 & empold & ~empFTold
* Case 3: some months in FT, some months in PT
replace derearn0_monthly =      derearn0/(12*(2.25*(FTrate)+(1-FTrate)))         if nmflag3 & empold & ~empFTold
replace derearn0_monthly = 2.25*derearn0/(12*(2.25*(FTrate)+(1-FTrate)))         if nmflag3          &  empFTold
* Case 4: all months of the year working either FT or PT
replace derearn0_monthly = derearn0/num_mos_work if (nmflag4 | nmflag5) & ~nmflag1 & ~nmflag2 & empold

generat derearn_monthlyJ    = derearn0_monthly
replace derearn_monthlyJ    = derearn_monthlyJ/CPI_Ta // already temporally deflated
generat lnderearn_monthlyJ  = ln(derearn_monthlyJ)

generat derearn_monthlyJb   = derearn0_monthly
replace derearn_monthlyJb   = derearn_monthlyJb/CPI_Tb // already temporally deflated
generat lnderearn_monthlyJb = ln(derearn_monthlyJb)

generat lnearnfinal   = lnearn
replace lnearnfinal   = lnderearn_monthly   if (mi(lnearn) | qmonthearnj1~=0 | qmonthearnj2~=0) & empold

generat lnearnfinalJ  = lnearnJ
replace lnearnfinalJ  = lnderearn_monthlyJ  if (mi(lnearn) | qmonthearnj1~=0 | qmonthearnj2~=0) & empold

generat lnearnfinalJb = lnearnJb
replace lnearnfinalJb = lnderearn_monthlyJb if (mi(lnearn) | qmonthearnj1~=0 | qmonthearnj2~=0) & empold


*---------------------------------------------------------------------------------------------------------------
* fix experience variable (no longer need to pay attention to gaps in the time series because I've dropped them)
*---------------------------------------------------------------------------------------------------------------
/*
forvalues X=1954/2003 {
	gen work2pQ`X' = inrange(wqc_yrtot_`X',3,4)
}
* ignore really old people
egen exper0_alt1 = rowtotal(work2pQ19?? work2pQ2000 work2pQ2001 work2pQ2002 work2pQ2003)
egen exper0_altC = rowtotal(wqc_yrtot_pre1953 wqc_yrtot_19?? wqc_yrtot_2000 wqc_yrtot_2001 wqc_yrtot_2002 wqc_yrtot_2003)
gen exper0_alt = exper0_altC/4
drop exper0_altC
 
lab var exper0_alt1 "years of working 2+ quarters (IRS)"
*/

forvalues X=1978/2003 {
	egen pos_der_`X' = rowmax(pos_der_fica_`X' pos_der_nonfica_`X')
}

egen exper0_1953_1977C = rowtotal(wqc_yrtot_pre1953 wqc_yrtot_1953 wqc_yrtot_1954 wqc_yrtot_1955 wqc_yrtot_1956 wqc_yrtot_1957 wqc_yrtot_1958 wqc_yrtot_1959 wqc_yrtot_196? wqc_yrtot_1971 wqc_yrtot_1972 wqc_yrtot_1973 wqc_yrtot_1974 wqc_yrtot_1975 wqc_yrtot_1976 wqc_yrtot_1977)
egen exper0_1978_      = rowtotal(pos_der_????)
gen exper0pre1953      = wqc_yrtot_pre1953/4
gen exper0_1953_1977   = exper0_1953_1977C/4

gen exper0_alt = exper0pre1953+exper0_1953_1977+exper0_1978_

lab var exper0_alt  "years of work experience (IRS)"

corr exper0_alt exper0 if wavemap==1 & SREFMON==4 & exper0_alt>0 & flag_in_der==1

bys ID (timmonth): generat lag_empFT = L.empFT
bys ID (timmonth): replace lag_empFT = 0 if _n==1
gen exper_temp = (1/12)*lag_empFT
bys ID (timmonth): gen exper_s = sum(exper_temp)
gen exper     = exper0      + exper_s
gen experAlt  = exper0_alt  + exper_s
drop exper_temp exper_s

bys ID (timmonth): generat lag_emp = L.emp
bys ID (timmonth): replace lag_emp = 0 if _n==1
gen exper_temp = (1/12)*lag_emp
bys ID (timmonth): gen exper_s = sum(exper_temp)
gen experFTorPT = exper0 + exper_s
drop exper_temp exper_s

bys ID (timmonth): gen byte firstper = _n==1
bys ID (timmonth): gen byte  lastper = _n==_N

*----------------------------------------------------------------------------------------
* generate location switching indicator
* note: switching indicator won't necessarily line up with actual switching date because
* switching is only observed at the wave level (and no "date moved" variable exists)
*----------------------------------------------------------------------------------------
bys ID (timmonth): generat byte switch_state  =  state[_n-1] ~=state[_n]  & ID[_n]==ID[_n-1] & _n>1
bys ID (timmonth): generat byte switch_county =  switch_state | (county[_n-1]~=county[_n] & state[_n-1]==state[_n] & ID[_n]==ID[_n-1] & _n>1)
bys ID (timmonth): generat byte switch_MSA    = (cbsa_code[_n-1]~=cbsa_code[_n] & ID[_n]==ID[_n-1] & _n>1)
bys ID (timmonth): generat byte switch_loc    =  switch_MSA
bys ID (timmonth): replace      switch_loc    =  1 if switch_county[_n] & cbsa_code[_n-1]==99999 & cbsa_code[_n]==99999 & ID[_n]==ID[_n-1] & _n>1 // add rural-rural movers as movers if they switched counties

generat migration_type = .
replace migration_type = 1 if switch_loc & ~switch_state
replace migration_type = 2 if switch_loc & switch_state
lab def vlmig 1 "Moved within state, out of MSA" 2 "Moved out of state and out of MSA"
lab val migration_type vlmig

* generate a range of dates when they could have switched location (can only be three possibilities since they are interviewed every 4 months)
gen byte switch_loc_range = .
bys ID (timmonth): replace switch_loc_range = 1 if switch_loc[_n+1]==1
bys ID (timmonth): replace switch_loc_range = 1 if switch_loc[_n+2]==1
bys ID (timmonth): replace switch_loc_range = 1 if switch_loc[_n+3]==1
recode switch_loc_range (. = 0)

gen byte switch_state_range = .
bys ID (timmonth): replace switch_state_range = 1 if switch_state[_n+1]==1
bys ID (timmonth): replace switch_state_range = 1 if switch_state[_n+2]==1
bys ID (timmonth): replace switch_state_range = 1 if switch_state[_n+3]==1
recode switch_state_range (. = 0)

gen byte switch_county_range = .
bys ID (timmonth): replace switch_county_range = 1 if switch_county[_n+1]==1
bys ID (timmonth): replace switch_county_range = 1 if switch_county[_n+2]==1
bys ID (timmonth): replace switch_county_range = 1 if switch_county[_n+3]==1
recode switch_county_range (. = 0)


*------------------------------------------------------
* Generate location variables
*------------------------------------------------------
generat choice29 =  .
replace choice29 =  1 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice29 =  2 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice29 =  3 if pop_cat==3 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice29 =  4 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice29 =  5 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice29 =  6 if pop_cat==3 &  inlist(stabb,"NJ","NY","PA") | regexm(cbsaname,"Philadelphia-Camden")
replace choice29 =  7 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice29 =  8 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington")
replace choice29 =  9 if pop_cat==3 & (inlist(stabb,"IN","IL","MI","OH","WI") | regexm(cbsaname,"Cincinnati")) & ~regexm(cbsaname,"Minneapolis") & ~regexm(cbsaname,"St. Louis")
replace choice29 = 10 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice29 = 11 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale")
replace choice29 = 12 if pop_cat==3 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Minneapolis")| regexm(cbsaname,"St. Louis")
replace choice29 = 13 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice29 = 14 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport")
replace choice29 = 15 if pop_cat==3 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV") & ~regexm(cbsaname,"Philadelphia-Camden")
replace choice29 = 16 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice29 = 17 if pop_cat==2 &  inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville")
replace choice29 = 18 if pop_cat==3 &  inlist(stabb,"AL","KY","MS","TN") & ~regexm(cbsaname,"Cincinnati")
replace choice29 = 19 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice29 = 20 if pop_cat==2 &  inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")
replace choice29 = 21 if pop_cat==3 &  inlist(stabb,"AR","LA","OK","TX")
replace choice29 = 22 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice29 = 23 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice29 = 24 if pop_cat==3 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice29 = 25 if pop_cat==1 &  inlist(stabb,"CA","OR","WA")
replace choice29 = 26 if pop_cat==2 &  inlist(stabb,"CA","OR","WA")
replace choice29 = 27 if pop_cat==3 &  inlist(stabb,"CA","OR","WA")
replace choice29 = 28 if pop_cat<=3 &  inlist(stabb,"AK")
replace choice29 = 29 if pop_cat<=3 &  inlist(stabb,"HI")
lab def vlchoice29 1 "NE small" 2 "NE medium" 3 "NE big" 4 "Mid Atl small" 5 "Mid Atl medium" 6 "Mid Atl big" 7 "E N Central small" 8 "E N Central medium" 9 "E N Central big" 10 "W N Central small" 11 "W N Central medium" 12 "W N Central big" 13 "S Atl small" 14 "S Atl medium" 15 "S Atl big" 16 "E S Central small" 17 "E S Central medium" 18 "E S Central big" 19 "W S Central small" 20 "W S Central medium" 21 "W S Central big" 22 "Mountain small" 23 "Mountain medium" 24 "Mountain big" 25 "Pacific small" 26 "Pacific medium" 27 "Pacific big" 28 "Alaska" 29 "Hawaii"
lab val choice29 vlchoice29


generat choice28 =  .
replace choice28 =  1 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice28 =  2 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice28 =  3 if pop_cat==3 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice28 =  4 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice28 =  5 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice28 =  6 if pop_cat==3 &  inlist(stabb,"NJ","NY","PA") | regexm(cbsaname,"Philadelphia-Camden")
replace choice28 =  7 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice28 =  8 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington")
replace choice28 =  9 if pop_cat==3 & (inlist(stabb,"IN","IL","MI","OH","WI") | regexm(cbsaname,"Cincinnati")) & ~regexm(cbsaname,"Minneapolis") & ~regexm(cbsaname,"St. Louis")
replace choice28 = 10 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice28 = 11 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale")
replace choice28 = 12 if pop_cat==3 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Minneapolis")| regexm(cbsaname,"St. Louis")
replace choice28 = 13 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice28 = 14 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport")
replace choice28 = 15 if pop_cat==3 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV") & ~regexm(cbsaname,"Philadelphia-Camden")
replace choice28 = 16 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice28 = 17 if pop_cat==2 &  inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville")
replace choice28 = 18 if pop_cat==3 &  inlist(stabb,"AL","KY","MS","TN") & ~regexm(cbsaname,"Cincinnati")
replace choice28 = 19 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice28 = 20 if pop_cat==2 &  inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")
replace choice28 = 21 if pop_cat==3 &  inlist(stabb,"AR","LA","OK","TX")
replace choice28 = 22 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice28 = 23 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice28 = 24 if pop_cat==3 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice28 = 25 if pop_cat==1 &  inlist(stabb,"CA","OR","WA")
replace choice28 = 26 if pop_cat==2 &  inlist(stabb,"CA","OR","WA")
replace choice28 = 27 if pop_cat==3 &  inlist(stabb,"CA","OR","WA")
replace choice28 = 28 if pop_cat<=3 &  inlist(stabb,"AK","HI")
lab def vlchoice28 1 "NE small" 2 "NE medium" 3 "NE big" 4 "Mid Atl small" 5 "Mid Atl medium" 6 "Mid Atl big" 7 "E N Central small" 8 "E N Central medium" 9 "E N Central big" 10 "W N Central small" 11 "W N Central medium" 12 "W N Central big" 13 "S Atl small" 14 "S Atl medium" 15 "S Atl big" 16 "E S Central small" 17 "E S Central medium" 18 "E S Central big" 19 "W S Central small" 20 "W S Central medium" 21 "W S Central big" 22 "Mountain small" 23 "Mountain medium" 24 "Mountain big" 25 "Pacific small" 26 "Pacific medium" 27 "Pacific big" 28 "Alaska/Hawaii"
lab val choice28 vlchoice28

generat choice27 =  .
replace choice27 =  1 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice27 =  2 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice27 =  3 if pop_cat==3 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice27 =  4 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice27 =  5 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice27 =  6 if pop_cat==3 &  inlist(stabb,"NJ","NY","PA") | regexm(cbsaname,"Philadelphia-Camden")
replace choice27 =  7 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice27 =  8 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington")
replace choice27 =  9 if pop_cat==3 & (inlist(stabb,"IN","IL","MI","OH","WI") | regexm(cbsaname,"Cincinnati")) & ~regexm(cbsaname,"Minneapolis") & ~regexm(cbsaname,"St. Louis")
replace choice27 = 10 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice27 = 11 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale")
replace choice27 = 12 if pop_cat==3 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Minneapolis")| regexm(cbsaname,"St. Louis")
replace choice27 = 13 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice27 = 14 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport")
replace choice27 = 15 if pop_cat==3 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV") & ~regexm(cbsaname,"Philadelphia-Camden")
replace choice27 = 16 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice27 = 17 if pop_cat==2 &  inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville")
replace choice27 = 18 if pop_cat==3 &  inlist(stabb,"AL","KY","MS","TN") & ~regexm(cbsaname,"Cincinnati")
replace choice27 = 19 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice27 = 20 if pop_cat==2 &  inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")
replace choice27 = 21 if pop_cat==3 &  inlist(stabb,"AR","LA","OK","TX")
replace choice27 = 22 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice27 = 23 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice27 = 24 if pop_cat==3 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice27 = 25 if pop_cat==1 &  inlist(stabb,"CA","OR","WA","HI","AK")
replace choice27 = 26 if pop_cat==2 &  inlist(stabb,"CA","OR","WA","HI","AK")
replace choice27 = 27 if pop_cat==3 &  inlist(stabb,"CA","OR","WA","HI","AK")
lab def vlchoice27 1 "NE small" 2 "NE medium" 3 "NE big" 4 "Mid Atl small" 5 "Mid Atl medium" 6 "Mid Atl big" 7 "E N Central small" 8 "E N Central medium" 9 "E N Central big" 10 "W N Central small" 11 "W N Central medium" 12 "W N Central big" 13 "S Atl small" 14 "S Atl medium" 15 "S Atl big" 16 "E S Central small" 17 "E S Central medium" 18 "E S Central big" 19 "W S Central small" 20 "W S Central medium" 21 "W S Central big" 22 "Mountain small" 23 "Mountain medium" 24 "Mountain big" 25 "Pacific small" 26 "Pacific medium" 27 "Pacific big"
lab val choice27 vlchoice27


generat choice48 =  .
replace choice48 =  1 if regexm(cbsaname,"Atlanta")
replace choice48 =  2 if regexm(cbsaname,"Austin-Round Rock")
replace choice48 =  3 if regexm(cbsaname,"Baltimore")
replace choice48 =  4 if regexm(cbsaname,"Boston")
replace choice48 =  5 if regexm(cbsaname,"Chicago")
replace choice48 =  6 if regexm(cbsaname,"Cincinnati")
replace choice48 =  7 if regexm(cbsaname,"Cleveland")
replace choice48 =  8 if regexm(cbsaname,"Dallas")
replace choice48 =  9 if regexm(cbsaname,"Denver")
replace choice48 = 10 if regexm(cbsaname,"Detroit")
replace choice48 = 11 if regexm(cbsaname,"Houston")
replace choice48 = 12 if regexm(cbsaname,"Indianapolis")
replace choice48 = 13 if regexm(cbsaname,"Kansas City")
replace choice48 = 14 if regexm(cbsaname,"Los Angeles")
replace choice48 = 15 if regexm(cbsaname,"Miami-Fort Lauderdale")
replace choice48 = 16 if regexm(cbsaname,"Milwaukee")
replace choice48 = 17 if regexm(cbsaname,"Minneapolis")
replace choice48 = 18 if regexm(cbsaname,"New York-Northern New Jersey")
replace choice48 = 19 if regexm(cbsaname,"Philadelphia")
replace choice48 = 20 if regexm(cbsaname,"Phoenix")
replace choice48 = 21 if regexm(cbsaname,"Pittsburgh, PA")
replace choice48 = 22 if regexm(cbsaname,"Portland-Vancouver")
replace choice48 = 23 if regexm(cbsaname,"Providence-New Bedford")
replace choice48 = 24 if regexm(cbsaname,"Riverside-San Bernardino")
replace choice48 = 25 if regexm(cbsaname,"San Diego")
replace choice48 = 26 if regexm(cbsaname,"San Francisco")
replace choice48 = 27 if regexm(cbsaname,"Seattle-Tacoma")
replace choice48 = 28 if regexm(cbsaname,"St. Louis")
replace choice48 = 29 if regexm(cbsaname,"Tampa-St. Petersburg")
replace choice48 = 30 if regexm(cbsaname,"Washington-Arlington")
replace choice48 = 31 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice48 = 32 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT") & ~regexm(cbsaname,"Providence-New Bedford")
replace choice48 = 33 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice48 = 34 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice48 = 35 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice48 = 36 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington") & ~regexm(cbsaname,"Indianapolis") & ~regexm(cbsaname,"Milwaukee")
replace choice48 = 37 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice48 = 38 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale") & ~regexm(cbsaname,"Kansas City")
replace choice48 = 39 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice48 = 40 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport")
replace choice48 = 41 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice48 = 42 if pop_cat>=2 &  inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville")
replace choice48 = 43 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice48 = 44 if pop_cat==2 & (inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")) & ~regexm(cbsaname,"Austin-Round Rock")
replace choice48 = 45 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice48 = 46 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice48 = 47 if pop_cat==1 &  inlist(stabb,"CA","OR","WA","AK","HI")
replace choice48 = 48 if pop_cat==2 &  inlist(stabb,"CA","OR","WA","AK","HI") & ~regexm(cbsaname,"Portland-Vancouver")
lab def vlchoice48 1 "Atlanta" 2 "Austin" 3 "Baltimore" 4 "Boston" 5 "Chicago" 6 "Cincinnati" 7 "Cleveland" 8 "Dallas" 9 "Denver" 10 "Detroit" 11 "Houston" 12 "Indianapolis" 13 "Kansas City" 14 "Los Angeles" 15 "Miami" 16 "Milwaukee" 17 "Minneapolis" 18 "New York" 19 "Philadelphia" 20 "Phoenix" 21 "Pittsburgh" 22 "Portland" 23 "Providence" 24 "Riverside" 25 "San Diego" 26 "San Francisco" 27 "Seattle" 28 "St. Louis" 29 "Tampa" 30 "Washington DC" 31 "NE small" 32 "NE medium" 33 "Mid Atl small" 34 "Mid Atl medium" 35 "E N Central small" 36 "E N Central medium" 37 "W N Central small" 38 "W N Central medium" 39 "S Atl small" 40 "S Atl medium" 41 "E S Central small" 42 "E S Central medium" 43 "W S Central small" 44 "W S Central medium" 45 "Mountain small" 46 "Mountain medium" 47 "Pacific small" 48 "Pacific medium"
lab val choice48 vlchoice48


generat choice49 =  .
replace choice49 =  1 if regexm(cbsaname,"Atlanta")
replace choice49 =  2 if regexm(cbsaname,"Austin-Round Rock")
replace choice49 =  3 if regexm(cbsaname,"Baltimore")
replace choice49 =  4 if regexm(cbsaname,"Boston")
replace choice49 =  5 if regexm(cbsaname,"Chicago")
replace choice49 =  6 if regexm(cbsaname,"Cincinnati")
replace choice49 =  7 if regexm(cbsaname,"Cleveland")
replace choice49 =  8 if regexm(cbsaname,"Dallas")
replace choice49 =  9 if regexm(cbsaname,"Denver")
replace choice49 = 10 if regexm(cbsaname,"Detroit")
replace choice49 = 11 if regexm(cbsaname,"Houston")
replace choice49 = 12 if regexm(cbsaname,"Indianapolis")
replace choice49 = 13 if regexm(cbsaname,"Kansas City")
replace choice49 = 14 if regexm(cbsaname,"Los Angeles")
replace choice49 = 15 if regexm(cbsaname,"Miami-Fort Lauderdale")
replace choice49 = 16 if regexm(cbsaname,"Milwaukee")
replace choice49 = 17 if regexm(cbsaname,"Minneapolis")
replace choice49 = 18 if regexm(cbsaname,"New York-Northern New Jersey")
replace choice49 = 19 if regexm(cbsaname,"Philadelphia")
replace choice49 = 20 if regexm(cbsaname,"Phoenix")
replace choice49 = 21 if regexm(cbsaname,"Pittsburgh, PA")
replace choice49 = 22 if regexm(cbsaname,"Portland-Vancouver")
replace choice49 = 23 if regexm(cbsaname,"Providence-New Bedford")
replace choice49 = 24 if regexm(cbsaname,"Riverside-San Bernardino")
replace choice49 = 25 if regexm(cbsaname,"San Diego")
replace choice49 = 26 if regexm(cbsaname,"San Francisco")
replace choice49 = 27 if regexm(cbsaname,"Seattle-Tacoma")
replace choice49 = 28 if regexm(cbsaname,"St. Louis")
replace choice49 = 29 if regexm(cbsaname,"Tampa-St. Petersburg")
replace choice49 = 30 if regexm(cbsaname,"Washington-Arlington")
replace choice49 = 31 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice49 = 32 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT") & ~regexm(cbsaname,"Providence-New Bedford")
replace choice49 = 33 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice49 = 34 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice49 = 35 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice49 = 36 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington") & ~regexm(cbsaname,"Indianapolis") & ~regexm(cbsaname,"Milwaukee")
replace choice49 = 37 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice49 = 38 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale") & ~regexm(cbsaname,"Kansas City")
replace choice49 = 39 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice49 = 40 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport")
replace choice49 = 41 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice49 = 42 if pop_cat>=2 &  inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville")
replace choice49 = 43 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice49 = 44 if pop_cat==2 & (inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")) & ~regexm(cbsaname,"Austin-Round Rock")
replace choice49 = 45 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice49 = 46 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice49 = 47 if pop_cat==1 &  inlist(stabb,"CA","OR","WA")
replace choice49 = 48 if pop_cat==2 &  inlist(stabb,"CA","OR","WA") & ~regexm(cbsaname,"Portland-Vancouver")
replace choice49 = 49 if inlist(stabb,"AK","HI")
lab def vlchoice49 1 "Atlanta" 2 "Austin" 3 "Baltimore" 4 "Boston" 5 "Chicago" 6 "Cincinnati" 7 "Cleveland" 8 "Dallas" 9 "Denver" 10 "Detroit" 11 "Houston" 12 "Indianapolis" 13 "Kansas City" 14 "Los Angeles" 15 "Miami" 16 "Milwaukee" 17 "Minneapolis" 18 "New York" 19 "Philadelphia" 20 "Phoenix" 21 "Pittsburgh" 22 "Portland" 23 "Providence" 24 "Riverside" 25 "San Diego" 26 "San Francisco" 27 "Seattle" 28 "St. Louis" 29 "Tampa" 30 "Washington DC" 31 "NE small" 32 "NE medium" 33 "Mid Atl small" 34 "Mid Atl medium" 35 "E N Central small" 36 "E N Central medium" 37 "W N Central small" 38 "W N Central medium" 39 "S Atl small" 40 "S Atl medium" 41 "E S Central small" 42 "E S Central medium" 43 "W S Central small" 44 "W S Central medium" 45 "Mountain small" 46 "Mountain medium" 47 "Pacific small" 48 "Pacific medium" 49 "Alaska/Hawaii"
lab val choice49 vlchoice49


generat choice50 =  .
replace choice50 =  1 if regexm(cbsaname,"Atlanta")
replace choice50 =  2 if regexm(cbsaname,"Austin-Round Rock")
replace choice50 =  3 if regexm(cbsaname,"Baltimore")
replace choice50 =  4 if regexm(cbsaname,"Boston")
replace choice50 =  5 if regexm(cbsaname,"Chicago")
replace choice50 =  6 if regexm(cbsaname,"Cincinnati")
replace choice50 =  7 if regexm(cbsaname,"Cleveland")
replace choice50 =  8 if regexm(cbsaname,"Dallas")
replace choice50 =  9 if regexm(cbsaname,"Denver")
replace choice50 = 10 if regexm(cbsaname,"Detroit")
replace choice50 = 11 if regexm(cbsaname,"Houston")
replace choice50 = 12 if regexm(cbsaname,"Indianapolis")
replace choice50 = 13 if regexm(cbsaname,"Kansas City")
replace choice50 = 14 if regexm(cbsaname,"Los Angeles")
replace choice50 = 15 if regexm(cbsaname,"Miami-Fort Lauderdale")
replace choice50 = 16 if regexm(cbsaname,"Milwaukee")
replace choice50 = 17 if regexm(cbsaname,"Minneapolis")
replace choice50 = 18 if regexm(cbsaname,"New York-Northern New Jersey")
replace choice50 = 19 if regexm(cbsaname,"Philadelphia")
replace choice50 = 20 if regexm(cbsaname,"Phoenix")
replace choice50 = 21 if regexm(cbsaname,"Pittsburgh, PA")
replace choice50 = 22 if regexm(cbsaname,"Portland-Vancouver")
replace choice50 = 23 if regexm(cbsaname,"Providence-New Bedford")
replace choice50 = 24 if regexm(cbsaname,"Riverside-San Bernardino")
replace choice50 = 25 if regexm(cbsaname,"San Diego")
replace choice50 = 26 if regexm(cbsaname,"San Francisco")
replace choice50 = 27 if regexm(cbsaname,"Seattle-Tacoma")
replace choice50 = 28 if regexm(cbsaname,"St. Louis")
replace choice50 = 29 if regexm(cbsaname,"Tampa-St. Petersburg")
replace choice50 = 30 if regexm(cbsaname,"Washington-Arlington")
replace choice50 = 31 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice50 = 32 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT") & ~regexm(cbsaname,"Providence-New Bedford")
replace choice50 = 33 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice50 = 34 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice50 = 35 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice50 = 36 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington") & ~regexm(cbsaname,"Indianapolis") & ~regexm(cbsaname,"Milwaukee")
replace choice50 = 37 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice50 = 38 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale") & ~regexm(cbsaname,"Kansas City")
replace choice50 = 39 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice50 = 40 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport")
replace choice50 = 41 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice50 = 42 if pop_cat>=2 &  inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville")
replace choice50 = 43 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice50 = 44 if pop_cat==2 & (inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")) & ~regexm(cbsaname,"Austin-Round Rock")
replace choice50 = 45 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice50 = 46 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice50 = 47 if pop_cat==1 &  inlist(stabb,"CA","OR","WA")
replace choice50 = 48 if pop_cat==2 &  inlist(stabb,"CA","OR","WA") & ~regexm(cbsaname,"Portland-Vancouver")
replace choice50 = 49 if inlist(stabb,"AK")
replace choice50 = 50 if inlist(stabb,"HI")
lab def vlchoice50 1 "Atlanta" 2 "Austin" 3 "Baltimore" 4 "Boston" 5 "Chicago" 6 "Cincinnati" 7 "Cleveland" 8 "Dallas" 9 "Denver" 10 "Detroit" 11 "Houston" 12 "Indianapolis" 13 "Kansas City" 14 "Los Angeles" 15 "Miami" 16 "Milwaukee" 17 "Minneapolis" 18 "New York" 19 "Philadelphia" 20 "Phoenix" 21 "Pittsburgh" 22 "Portland" 23 "Providence" 24 "Riverside" 25 "San Diego" 26 "San Francisco" 27 "Seattle" 28 "St. Louis" 29 "Tampa" 30 "Washington DC" 31 "NE small" 32 "NE medium" 33 "Mid Atl small" 34 "Mid Atl medium" 35 "E N Central small" 36 "E N Central medium" 37 "W N Central small" 38 "W N Central medium" 39 "S Atl small" 40 "S Atl medium" 41 "E S Central small" 42 "E S Central medium" 43 "W S Central small" 44 "W S Central medium" 45 "Mountain small" 46 "Mountain medium" 47 "Pacific small" 48 "Pacific medium" 49 "Alaska" 50 "Hawaii"
lab val choice50 vlchoice50


generat choice53 =  .
replace choice53 =  1 if regexm(cbsaname,"Atlanta")
replace choice53 =  2 if regexm(cbsaname,"Austin-Round Rock")
replace choice53 =  3 if regexm(cbsaname,"Baltimore")
replace choice53 =  4 if regexm(cbsaname,"Boston")
replace choice53 =  5 if regexm(cbsaname,"Chicago")
replace choice53 =  6 if regexm(cbsaname,"Cincinnati")
replace choice53 =  7 if regexm(cbsaname,"Cleveland")
replace choice53 =  8 if regexm(cbsaname,"Columbus, OH")
replace choice53 =  9 if regexm(cbsaname,"Dallas")
replace choice53 = 10 if regexm(cbsaname,"Denver")
replace choice53 = 11 if regexm(cbsaname,"Detroit")
replace choice53 = 12 if regexm(cbsaname,"Houston")
replace choice53 = 13 if regexm(cbsaname,"Indianapolis")
replace choice53 = 14 if regexm(cbsaname,"Kansas City")
replace choice53 = 15 if regexm(cbsaname,"Knoxville")
replace choice53 = 16 if regexm(cbsaname,"Los Angeles")
replace choice53 = 17 if regexm(cbsaname,"Miami-Fort Lauderdale")
replace choice53 = 18 if regexm(cbsaname,"Milwaukee")
replace choice53 = 19 if regexm(cbsaname,"Minneapolis")
replace choice53 = 20 if regexm(cbsaname,"New York-Northern New Jersey")
replace choice53 = 21 if regexm(cbsaname,"Philadelphia")
replace choice53 = 22 if regexm(cbsaname,"Phoenix")
replace choice53 = 23 if regexm(cbsaname,"Pittsburgh, PA")
replace choice53 = 24 if regexm(cbsaname,"Portland-Vancouver")
replace choice53 = 25 if regexm(cbsaname,"Providence-New Bedford")
replace choice53 = 26 if regexm(cbsaname,"Richmond, VA")
replace choice53 = 27 if regexm(cbsaname,"Riverside-San Bernardino")
replace choice53 = 28 if regexm(cbsaname,"Sacramento")
replace choice53 = 29 if regexm(cbsaname,"San Diego")
replace choice53 = 30 if regexm(cbsaname,"San Francisco")
replace choice53 = 31 if regexm(cbsaname,"Seattle-Tacoma")
replace choice53 = 32 if regexm(cbsaname,"St. Louis")
replace choice53 = 33 if regexm(cbsaname,"Tampa-St. Petersburg")
replace choice53 = 34 if regexm(cbsaname,"Virginia Beach")
replace choice53 = 35 if regexm(cbsaname,"Washington-Arlington")
replace choice53 = 36 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice53 = 37 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT") & ~regexm(cbsaname,"Providence-New Bedford")
replace choice53 = 38 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice53 = 39 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice53 = 40 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice53 = 41 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington") & ~regexm(cbsaname,"Indianapolis") & ~regexm(cbsaname,"Milwaukee")  & ~regexm(cbsaname,"Columbus, OH")
replace choice53 = 42 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice53 = 43 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale") & ~regexm(cbsaname,"Kansas City")
replace choice53 = 44 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice53 = 45 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport") & ~regexm(cbsaname,"Richmond, VA") & ~regexm(cbsaname,"Virginia Beach")
replace choice53 = 46 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice53 = 47 if pop_cat>=2 & (inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville"))  & ~regexm(cbsaname,"Knoxville")
replace choice53 = 48 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice53 = 49 if pop_cat==2 & (inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")) & ~regexm(cbsaname,"Austin-Round Rock")
replace choice53 = 50 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice53 = 51 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice53 = 52 if pop_cat==1 &  inlist(stabb,"CA","OR","WA","AK","HI")
replace choice53 = 53 if pop_cat==2 &  inlist(stabb,"CA","OR","WA","AK","HI") & ~regexm(cbsaname,"Portland-Vancouver") & ~regexm(cbsaname,"Sacramento")
lab def vlchoice53 1 "Atlanta" 2 "Austin" 3 "Baltimore" 4 "Boston" 5 "Chicago" 6 "Cincinnati" 7 "Cleveland" 8 "Columbus" 9 "Dallas" 10 "Denver" 11 "Detroit" 12 "Houston" 13 "Indianapolis" 14 "Kansas City" 15 "Knoxville" 16 "Los Angeles" 17 "Miami" 18 "Milwaukee" 19 "Minneapolis" 20 "New York" 21 "Philadelphia" 22 "Phoenix" 23 "Pittsburgh" 24 "Portland" 25 "Providence" 26 "Richmond" 27 "Riverside" 28 "Sacramento" 29 "San Diego" 30 "San Francisco" 31 "Seattle" 32 "St. Louis" 33 "Tampa" 34 "Virginia Beach" 35 "Washington DC" 36 "NE small" 37 "NE medium" 38 "Mid Atl small" 39 "Mid Atl medium" 40 "E N Central small" 41 "E N Central medium" 42 "W N Central small" 43 "W N Central medium" 44 "S Atl small" 45 "S Atl medium" 46 "E S Central small" 47 "E S Central medium" 48 "W S Central small" 49 "W S Central medium" 50 "Mountain small" 51 "Mountain medium" 52 "Pacific small" 53 "Pacific medium"
lab val choice53 vlchoice53


generat choice54b =  .
replace choice54b =  1 if regexm(cbsaname,"Atlanta")
replace choice54b =  2 if regexm(cbsaname,"Austin-Round Rock")
replace choice54b =  3 if regexm(cbsaname,"Baltimore")
replace choice54b =  4 if regexm(cbsaname,"Boston")
replace choice54b =  5 if regexm(cbsaname,"Chicago")
replace choice54b =  6 if regexm(cbsaname,"Cincinnati")
replace choice54b =  7 if regexm(cbsaname,"Cleveland")
replace choice54b =  8 if regexm(cbsaname,"Columbus, OH")
replace choice54b =  9 if regexm(cbsaname,"Dallas")
replace choice54b = 10 if regexm(cbsaname,"Denver")
replace choice54b = 11 if regexm(cbsaname,"Detroit")
replace choice54b = 12 if regexm(cbsaname,"Houston")
replace choice54b = 13 if regexm(cbsaname,"Indianapolis")
replace choice54b = 14 if regexm(cbsaname,"Kansas City")
replace choice54b = 15 if regexm(cbsaname,"Knoxville")
replace choice54b = 16 if regexm(cbsaname,"Los Angeles")
replace choice54b = 17 if regexm(cbsaname,"Miami-Fort Lauderdale")
replace choice54b = 18 if regexm(cbsaname,"Milwaukee")
replace choice54b = 19 if regexm(cbsaname,"Minneapolis")
replace choice54b = 20 if regexm(cbsaname,"New York-Northern New Jersey")
replace choice54b = 21 if regexm(cbsaname,"Philadelphia")
replace choice54b = 22 if regexm(cbsaname,"Phoenix")
replace choice54b = 23 if regexm(cbsaname,"Pittsburgh, PA")
replace choice54b = 24 if regexm(cbsaname,"Portland-Vancouver")
replace choice54b = 25 if regexm(cbsaname,"Providence-New Bedford")
replace choice54b = 26 if regexm(cbsaname,"Richmond, VA")
replace choice54b = 27 if regexm(cbsaname,"Riverside-San Bernardino")
replace choice54b = 28 if regexm(cbsaname,"Sacramento")
replace choice54b = 29 if regexm(cbsaname,"San Diego")
replace choice54b = 30 if regexm(cbsaname,"San Francisco")
replace choice54b = 31 if regexm(cbsaname,"Seattle-Tacoma")
replace choice54b = 32 if regexm(cbsaname,"St. Louis")
replace choice54b = 33 if regexm(cbsaname,"Tampa-St. Petersburg")
replace choice54b = 34 if regexm(cbsaname,"Virginia Beach")
replace choice54b = 35 if regexm(cbsaname,"Washington-Arlington")
replace choice54b = 36 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice54b = 37 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT") & ~regexm(cbsaname,"Providence-New Bedford")
replace choice54b = 38 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice54b = 39 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice54b = 40 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice54b = 41 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington") & ~regexm(cbsaname,"Indianapolis") & ~regexm(cbsaname,"Milwaukee")  & ~regexm(cbsaname,"Columbus, OH")
replace choice54b = 42 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice54b = 43 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale") & ~regexm(cbsaname,"Kansas City")
replace choice54b = 44 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice54b = 45 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport") & ~regexm(cbsaname,"Richmond, VA") & ~regexm(cbsaname,"Virginia Beach")
replace choice54b = 46 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice54b = 47 if pop_cat>=2 & (inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville"))  & ~regexm(cbsaname,"Knoxville")
replace choice54b = 48 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice54b = 49 if pop_cat==2 & (inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")) & ~regexm(cbsaname,"Austin-Round Rock")
replace choice54b = 50 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice54b = 51 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice54b = 52 if pop_cat==1 &  inlist(stabb,"CA","OR","WA","AK","HI")
replace choice54b = 53 if pop_cat==2 &  inlist(stabb,"CA","OR","WA") & ~regexm(cbsaname,"Portland-Vancouver") & ~regexm(cbsaname,"Sacramento")
replace choice54b = 54 if inlist(stabb,"AK","HI")
lab def vlchoice54b 1 "Atlanta" 2 "Austin" 3 "Baltimore" 4 "Boston" 5 "Chicago" 6 "Cincinnati" 7 "Cleveland" 8 "Columbus" 9 "Dallas" 10 "Denver" 11 "Detroit" 12 "Houston" 13 "Indianapolis" 14 "Kansas City" 15 "Knoxville" 16 "Los Angeles" 17 "Miami" 18 "Milwaukee" 19 "Minneapolis" 20 "New York" 21 "Philadelphia" 22 "Phoenix" 23 "Pittsburgh" 24 "Portland" 25 "Providence" 26 "Richmond" 27 "Riverside" 28 "Sacramento" 29 "San Diego" 30 "San Francisco" 31 "Seattle" 32 "St. Louis" 33 "Tampa" 34 "Virginia Beach" 35 "Washington DC" 36 "NE small" 37 "NE medium" 38 "Mid Atl small" 39 "Mid Atl medium" 40 "E N Central small" 41 "E N Central medium" 42 "W N Central small" 43 "W N Central medium" 44 "S Atl small" 45 "S Atl medium" 46 "E S Central small" 47 "E S Central medium" 48 "W S Central small" 49 "W S Central medium" 50 "Mountain small" 51 "Mountain medium" 52 "Pacific small" 53 "Pacific medium" 54 "Alaska/Hawaii"
lab val choice54b vlchoice54b


generat choice55 =  .
replace choice55 =  1 if regexm(cbsaname,"Atlanta")
replace choice55 =  2 if regexm(cbsaname,"Austin-Round Rock")
replace choice55 =  3 if regexm(cbsaname,"Baltimore")
replace choice55 =  4 if regexm(cbsaname,"Boston")
replace choice55 =  5 if regexm(cbsaname,"Chicago")
replace choice55 =  6 if regexm(cbsaname,"Cincinnati")
replace choice55 =  7 if regexm(cbsaname,"Cleveland")
replace choice55 =  8 if regexm(cbsaname,"Columbus, OH")
replace choice55 =  9 if regexm(cbsaname,"Dallas")
replace choice55 = 10 if regexm(cbsaname,"Denver")
replace choice55 = 11 if regexm(cbsaname,"Detroit")
replace choice55 = 12 if regexm(cbsaname,"Houston")
replace choice55 = 13 if regexm(cbsaname,"Indianapolis")
replace choice55 = 14 if regexm(cbsaname,"Kansas City")
replace choice55 = 15 if regexm(cbsaname,"Knoxville")
replace choice55 = 16 if regexm(cbsaname,"Los Angeles")
replace choice55 = 17 if regexm(cbsaname,"Miami-Fort Lauderdale")
replace choice55 = 18 if regexm(cbsaname,"Milwaukee")
replace choice55 = 19 if regexm(cbsaname,"Minneapolis")
replace choice55 = 20 if regexm(cbsaname,"New York-Northern New Jersey")
replace choice55 = 21 if regexm(cbsaname,"Philadelphia")
replace choice55 = 22 if regexm(cbsaname,"Phoenix")
replace choice55 = 23 if regexm(cbsaname,"Pittsburgh, PA")
replace choice55 = 24 if regexm(cbsaname,"Portland-Vancouver")
replace choice55 = 25 if regexm(cbsaname,"Providence-New Bedford")
replace choice55 = 26 if regexm(cbsaname,"Richmond, VA")
replace choice55 = 27 if regexm(cbsaname,"Riverside-San Bernardino")
replace choice55 = 28 if regexm(cbsaname,"Sacramento")
replace choice55 = 29 if regexm(cbsaname,"San Diego")
replace choice55 = 30 if regexm(cbsaname,"San Francisco")
replace choice55 = 31 if regexm(cbsaname,"Seattle-Tacoma")
replace choice55 = 32 if regexm(cbsaname,"St. Louis")
replace choice55 = 33 if regexm(cbsaname,"Tampa-St. Petersburg")
replace choice55 = 34 if regexm(cbsaname,"Virginia Beach")
replace choice55 = 35 if regexm(cbsaname,"Washington-Arlington")
replace choice55 = 36 if pop_cat==1 &  inlist(stabb,"CT","ME","MA","NH","RI","VT")
replace choice55 = 37 if pop_cat==2 &  inlist(stabb,"CT","ME","MA","NH","RI","VT") & ~regexm(cbsaname,"Providence-New Bedford")
replace choice55 = 38 if pop_cat==1 &  inlist(stabb,"NJ","NY","PA")
replace choice55 = 39 if pop_cat==2 &  inlist(stabb,"NJ","NY","PA") & ~regexm(cbsaname,"Youngstown")
replace choice55 = 40 if pop_cat==1 &  inlist(stabb,"IN","IL","MI","OH","WI")
replace choice55 = 41 if pop_cat==2 & (inlist(stabb,"IN","IL","MI","OH","WI")| regexm(cbsaname,"Youngstown")) & ~regexm(cbsaname,"Louisville") & ~regexm(cbsaname,"Duluth") & ~regexm(cbsaname,"Huntington") & ~regexm(cbsaname,"Indianapolis") & ~regexm(cbsaname,"Milwaukee")  & ~regexm(cbsaname,"Columbus, OH")
replace choice55 = 42 if pop_cat==1 &  inlist(stabb,"IA","KS","MN","MO","NE","ND","SD")
replace choice55 = 43 if pop_cat==2 & (inlist(stabb,"IA","KS","MN","MO","NE","ND","SD") | regexm(cbsaname,"Duluth")) & ~regexm(cbsaname,"Fayetteville-Springdale") & ~regexm(cbsaname,"Kansas City")
replace choice55 = 44 if pop_cat==1 &  inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")
replace choice55 = 45 if pop_cat==2 & (inlist(stabb,"DE","DC","FL","GA","MD","NC","SC","VA","WV")| regexm(cbsaname,"Huntington")) & ~regexm(cbsaname,"Chattanooga") & ~regexm(cbsaname,"Kingsport") & ~regexm(cbsaname,"Richmond, VA") & ~regexm(cbsaname,"Virginia Beach")
replace choice55 = 46 if pop_cat==1 &  inlist(stabb,"AL","KY","MS","TN")
replace choice55 = 47 if pop_cat>=2 & (inlist(stabb,"AL","KY","MS","TN")| regexm(cbsaname,"Chattanooga") | regexm(cbsaname,"Kingsport") | regexm(cbsaname,"Louisville"))  & ~regexm(cbsaname,"Knoxville")
replace choice55 = 48 if pop_cat==1 &  inlist(stabb,"AR","LA","OK","TX")
replace choice55 = 49 if pop_cat==2 & (inlist(stabb,"AR","LA","OK","TX")| regexm(cbsaname,"Fayetteville-Springdale")) & ~regexm(cbsaname,"Austin-Round Rock")
replace choice55 = 50 if pop_cat==1 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice55 = 51 if pop_cat==2 &  inlist(stabb,"AZ","CO","ID","NM","MT","UT","NV","WY")
replace choice55 = 52 if pop_cat==1 &  inlist(stabb,"CA","OR","WA","AK","HI")
replace choice55 = 53 if pop_cat==2 &  inlist(stabb,"CA","OR","WA","AK","HI") & ~regexm(cbsaname,"Portland-Vancouver") & ~regexm(cbsaname,"Sacramento")
replace choice55 = 54 if inlist(stabb,"AK")
replace choice55 = 55 if inlist(stabb,"HI")
lab def vlchoice55 1 "Atlanta" 2 "Austin" 3 "Baltimore" 4 "Boston" 5 "Chicago" 6 "Cincinnati" 7 "Cleveland" 8 "Columbus" 9 "Dallas" 10 "Denver" 11 "Detroit" 12 "Houston" 13 "Indianapolis" 14 "Kansas City" 15 "Knoxville" 16 "Los Angeles" 17 "Miami" 18 "Milwaukee" 19 "Minneapolis" 20 "New York" 21 "Philadelphia" 22 "Phoenix" 23 "Pittsburgh" 24 "Portland" 25 "Providence" 26 "Richmond" 27 "Riverside" 28 "Sacramento" 29 "San Diego" 30 "San Francisco" 31 "Seattle" 32 "St. Louis" 33 "Tampa" 34 "Virginia Beach" 35 "Washington DC" 36 "NE small" 37 "NE medium" 38 "Mid Atl small" 39 "Mid Atl medium" 40 "E N Central small" 41 "E N Central medium" 42 "W N Central small" 43 "W N Central medium" 44 "S Atl small" 45 "S Atl medium" 46 "E S Central small" 47 "E S Central medium" 48 "W S Central small" 49 "W S Central medium" 50 "Mountain small" 51 "Mountain medium" 52 "Pacific small" 53 "Pacific medium" 54 "Alaska" 55 "Hawaii"
lab val choice55 vlchoice55

*-----------------------------------------
* generate the choice set with work & home
*-----------------------------------------
generat choice54 = choice27    if inlf==1
replace choice54 = choice27+27 if inlf==0
lab def vlchoice54 1 "NE small work" 2 "NE medium work" 3 "NE big work" 4 "Mid Atl small work" 5 "Mid Atl medium work" 6 "Mid Atl big work" 7 "E N Central small work" 8 "E N Central medium work" 9 "E N Central big work" 10 "W N Central small work" 11 "W N Central medium work" 12 "W N Central big work" 13 "S Atl small work" 14 "S Atl medium work" 15 "S Atl big work" 16 "E S Central small work" 17 "E S Central medium work" 18 "E S Central big work" 19 "W S Central small work" 20 "W S Central medium work" 21 "W S Central big work" 22 "Mountain small work" 23 "Mountain medium work" 24 "Mountain big work" 25 "Pacific small work" 26 "Pacific medium work" 27 "Pacific big work"28 "NE small home" 29 "NE medium home" 30 "NE big home" 31 "Mid Atl small home" 32 "Mid Atl medium home" 33 "Mid Atl big home" 34 "E N Central small home" 35 "E N Central medium home" 36 "E N Central big home" 37 "W N Central small home" 38 "W N Central medium home" 39 "W N Central big home" 40 "S Atl small home" 41 "S Atl medium home" 42 "S Atl big home" 43 "E S Central small home" 44 "E S Central medium home" 45 "E S Central big home" 46 "W S Central small home" 47 "W S Central medium home" 48 "W S Central big home" 49 "Mountain small home" 50 "Mountain medium home" 51 "Mountain big home" 52 "Pacific small home" 53 "Pacific medium home" 54 "Pacific big home"
lab val choice54 vlchoice54

generat choice56  = choice28    if inlf==1
replace choice56  = choice28+28 if inlf==0
lab def vlchoice56 1 "NE small work" 2 "NE medium work" 3 "NE big work" 4 "Mid Atl small work" 5 "Mid Atl medium work" 6 "Mid Atl big work" 7 "E N Central small work" 8 "E N Central medium work" 9 "E N Central big work" 10 "W N Central small work" 11 "W N Central medium work" 12 "W N Central big work" 13 "S Atl small work" 14 "S Atl medium work" 15 "S Atl big work" 16 "E S Central small work" 17 "E S Central medium work" 18 "E S Central big work" 19 "W S Central small work" 20 "W S Central medium work" 21 "W S Central big work" 22 "Mountain small work" 23 "Mountain medium work" 24 "Mountain big work" 25 "Pacific small work" 26 "Pacific medium work" 27 "Pacific big work" 28 "Alaska/Hawaii work" 29 "NE small home" 30 "NE medium home" 31 "NE big home" 32 "Mid Atl small home" 33 "Mid Atl medium home" 34 "Mid Atl big home" 35 "E N Central small home" 36 "E N Central medium home" 37 "E N Central big home" 38 "W N Central small home" 39 "W N Central medium home" 40 "W N Central big home" 41 "S Atl small home" 42 "S Atl medium home" 43 "S Atl big home" 44 "E S Central small home" 45 "E S Central medium home" 46 "E S Central big home" 47 "W S Central small home" 48 "W S Central medium home" 49 "W S Central big home" 50 "Mountain small home" 51 "Mountain medium home" 52 "Mountain big home" 53 "Pacific small home" 54 "Pacific medium home" 55 "Pacific big home" 56 "Alaska/Hawaii home"
lab val choice56 vlchoice56

generat choice58  = choice29    if inlf==1
replace choice58  = choice29+29 if inlf==0
lab def vlchoice58 1 "NE small work" 2 "NE medium work" 3 "NE big work" 4 "Mid Atl small work" 5 "Mid Atl medium work" 6 "Mid Atl big work" 7 "E N Central small work" 8 "E N Central medium work" 9 "E N Central big work" 10 "W N Central small work" 11 "W N Central medium work" 12 "W N Central big work" 13 "S Atl small work" 14 "S Atl medium work" 15 "S Atl big work" 16 "E S Central small work" 17 "E S Central medium work" 18 "E S Central big work" 19 "W S Central small work" 20 "W S Central medium work" 21 "W S Central big work" 22 "Mountain small work" 23 "Mountain medium work" 24 "Mountain big work" 25 "Pacific small work" 26 "Pacific medium work" 27 "Pacific big work" 28 "Alaska work" 29 "Hawaii work" 30 "NE small home" 31 "NE medium home" 32 "NE big home" 33 "Mid Atl small home" 34 "Mid Atl medium home" 35 "Mid Atl big home" 36 "E N Central small home" 37 "E N Central medium home" 38 "E N Central big home" 39 "W N Central small home" 40 "W N Central medium home" 41 "W N Central big home" 42 "S Atl small home" 43 "S Atl medium home" 44 "S Atl big home" 45 "E S Central small home" 46 "E S Central medium home" 47 "E S Central big home" 48 "W S Central small home" 49 "W S Central medium home" 50 "W S Central big home" 51 "Mountain small home" 52 "Mountain medium home" 53 "Mountain big home" 54 "Pacific small home" 55 "Pacific medium home" 56 "Pacific big home" 57 "Alaska home" 58 "Hawaii home"
lab val choice58 vlchoice58

generat choice96  = choice48    if inlf==1
replace choice96  = choice48+48 if inlf==0
lab def vlchoice96 1 "Atlanta work" 2 "Austin work" 3 "Baltimore work" 4 "Boston work" 5 "Chicago work" 6 "Cincinnati work" 7 "Cleveland work" 8 "Dallas work" 9 "Denver work" 10 "Detroit work" 11 "Houston work" 12 "Indianapolis work" 13 "Kansas City work" 14 "Los Angeles work" 15 "Miami work" 16 "Milwaukee work" 17 "Minneapolis work" 18 "New York work" 19 "Philadelphia work" 20 "Phoenix work" 21 "Pittsburgh work" 22 "Portland work" 23 "Providence work" 24 "Riverside work" 25 "San Diego work" 26 "San Francisco work" 27 "Seattle work" 28 "St. Louis work" 29 "Tampa work" 30 "Washington DC work" 31 "NE small work" 32 "NE medium work" 33 "Mid Atl small work" 34 "Mid Atl medium work" 35 "N Central small work" 36 "N Central medium work" 37 "NW Central small work" 38 "NW Central medium work" 39 "S Mid-Atl small work" 40 "S Mid-Atl medium work" 41 "S Central small work" 42 "S Central medium work" 43 "SW Central small work" 44 "SW Central medium work" 45 "Mountain small work" 46 "Mountain medium work" 47 "Pacific small work" 48 "Pacific medium work" 49 "Atlanta home" 50 "Austin home" 51 "Baltimore home" 52 "Boston home" 53 "Chicago home" 54 "Cincinnati home" 55 "Cleveland home" 56 "Dallas home" 57 "Denver home" 58 "Detroit home" 59 "Houston home" 60 "Indianapolis home" 61 "Kansas City home" 62 "Los Angeles home" 63 "Miami home" 64 "Milwaukee home" 65 "Minneapolis home" 66 "New York home" 67 "Philadelphia home" 68 "Phoenix home" 69 "Pittsburgh home" 70 "Portland home" 71 "Providence home" 72 "Riverside home" 73 "San Diego home" 74 "San Francisco home" 75 "Seattle home" 76 "St. Louis home" 77 "Tampa home" 78 "Washington DC home" 79 "NE small home" 80 "NE medium home" 81 "Mid Atl small home" 82 "Mid Atl medium home" 83 "E N Central small home" 84 "E N Central medium home" 85 "W N Central small home" 86 "W N Central medium home" 87 "S Atl small home" 88 "S Atl medium home" 89 "E S Central small home" 90 "E S Central medium home" 91 "W S Central small home" 92 "W S Central medium home" 93 "Mountain small home" 94 "Mountain medium home" 95 "Pacific small home" 96 "Pacific medium home"
lab val choice96 vlchoice96

generat choice98  = choice49    if inlf==1
replace choice98  = choice49+49 if inlf==0
lab def vlchoice98 1 "Atlanta work" 2 "Austin work" 3 "Baltimore work" 4 "Boston work" 5 "Chicago work" 6 "Cincinnati work" 7 "Cleveland work" 8 "Dallas work" 9 "Denver work" 10 "Detroit work" 11 "Houston work" 12 "Indianapolis work" 13 "Kansas City work" 14 "Los Angeles work" 15 "Miami work" 16 "Milwaukee work" 17 "Minneapolis work" 18 "New York work" 19 "Philadelphia work" 20 "Phoenix work" 21 "Pittsburgh work" 22 "Portland work" 23 "Providence work" 24 "Riverside work" 25 "San Diego work" 26 "San Francisco work" 27 "Seattle work" 28 "St. Louis work" 29 "Tampa work" 30 "Washington DC work" 31 "NE small work" 32 "NE medium work" 33 "Mid Atl small work" 34 "Mid Atl medium work" 35 "N Central small work" 36 "N Central medium work" 37 "NW Central small work" 38 "NW Central medium work" 39 "S Mid-Atl small work" 40 "S Mid-Atl medium work" 41 "S Central small work" 42 "S Central medium work" 43 "SW Central small work" 44 "SW Central medium work" 45 "Mountain small work" 46 "Mountain medium work" 47 "Pacific small work" 48 "Pacific medium work" 49 "Alaska/Hawaii work" 50 "Atlanta home" 51 "Austin home" 52 "Baltimore home" 53 "Boston home" 54 "Chicago home" 55 "Cincinnati home" 56 "Cleveland home" 57 "Dallas home" 58 "Denver home" 59 "Detroit home" 60 "Houston home" 61 "Indianapolis home" 62 "Kansas City home" 63 "Los Angeles home" 64 "Miami home" 65 "Milwaukee home" 66 "Minneapolis home" 67 "New York home" 68 "Philadelphia home" 69 "Phoenix home" 70 "Pittsburgh home" 71 "Portland home" 72 "Providence home" 73 "Riverside home" 74 "San Diego home" 75 "San Francisco home" 76 "Seattle home" 77 "St. Louis home" 78 "Tampa home" 79 "Washington DC home" 80 "NE small home" 81 "NE medium home" 82 "Mid Atl small home" 83 "Mid Atl medium home" 84 "E N Central small home" 85 "E N Central medium home" 86 "W N Central small home" 87 "W N Central medium home" 88 "S Atl small home" 89 "S Atl medium home" 90 "E S Central small home" 91 "E S Central medium home" 92 "W S Central small home" 93 "W S Central medium home" 94 "Mountain small home" 95 "Mountain medium home" 96 "Pacific small home" 97 "Pacific medium home" 98 "Alaska/Hawaii home"
lab val choice98 vlchoice98

generat choice100 = choice50    if inlf==1
replace choice100 = choice50+50 if inlf==0
lab def vlchoice100 1 "Atlanta work" 2 "Austin work" 3 "Baltimore work" 4 "Boston work" 5 "Chicago work" 6 "Cincinnati work" 7 "Cleveland work" 8 "Dallas work" 9 "Denver work" 10 "Detroit work" 11 "Houston work" 12 "Indianapolis work" 13 "Kansas City work" 14 "Los Angeles work" 15 "Miami work" 16 "Milwaukee work" 17 "Minneapolis work" 18 "New York work" 19 "Philadelphia work" 20 "Phoenix work" 21 "Pittsburgh work" 22 "Portland work" 23 "Providence work" 24 "Riverside work" 25 "San Diego work" 26 "San Francisco work" 27 "Seattle work" 28 "St. Louis work" 29 "Tampa work" 30 "Washington DC work" 31 "NE small work" 32 "NE medium work" 33 "Mid Atl small work" 34 "Mid Atl medium work" 35 "N Central small work" 36 "N Central medium work" 37 "NW Central small work" 38 "NW Central medium work" 39 "S Mid-Atl small work" 40 "S Mid-Atl medium work" 41 "S Central small work" 42 "S Central medium work" 43 "SW Central small work" 44 "SW Central medium work" 45 "Mountain small work" 46 "Mountain medium work" 47 "Pacific small work" 48 "Pacific medium work" 49 "Alaska work" 50 "Hawaii work" 51 "Atlanta home" 52 "Austin home" 53 "Baltimore home" 54 "Boston home" 55 "Chicago home" 56 "Cincinnati home" 57 "Cleveland home" 58 "Dallas home" 59 "Denver home" 60 "Detroit home" 61 "Houston home" 62 "Indianapolis home" 63 "Kansas City home" 64 "Los Angeles home" 65 "Miami home" 66 "Milwaukee home" 67 "Minneapolis home" 68 "New York home" 69 "Philadelphia home" 70 "Phoenix home" 71 "Pittsburgh home" 72 "Portland home" 73 "Providence home" 74 "Riverside home" 75 "San Diego home" 76 "San Francisco home" 77 "Seattle home" 78 "St. Louis home" 79 "Tampa home" 80 "Washington DC home" 81 "NE small home" 82 "NE medium home" 83 "Mid Atl small home" 84 "Mid Atl medium home" 85 "E N Central small home" 86 "E N Central medium home" 87 "W N Central small home" 88 "W N Central medium home" 89 "S Atl small home" 90 "S Atl medium home" 91 "E S Central small home" 92 "E S Central medium home" 93 "W S Central small home" 94 "W S Central medium home" 95 "Mountain small home" 96 "Mountain medium home" 97 "Pacific small home" 98 "Pacific medium home" 99 "Alaska home" 100 "Hawaii home"
lab val choice100 vlchoice100

generat choice106 = choice53    if inlf==1
replace choice106 = choice53+53 if inlf==0
lab def vlchoice106 1 "Atlanta work" 2 "Austin work" 3 "Baltimore work" 4 "Boston work" 5 "Chicago work" 6 "Cincinnati work" 7 "Cleveland work" 8 "Columbus work" 9 "Dallas work" 10 "Denver work" 11 "Detroit work" 12 "Houston work" 13 "Indianapolis work" 14 "Kansas City work" 15 "Knoxville work" 16 "Los Angeles work" 17 "Miami work" 18 "Milwaukee work" 19 "Minneapolis work" 20 "New York work" 21 "Philadelphia work" 22 "Phoenix work" 23 "Pittsburgh work" 24 "Portland work" 25 "Providence"26 "Richmond work" 27 "Riverside work" 28 "Sacramento work" 29 "San Diego work" 30 "San Francisco work" 31 "Seattle work" 32 "St. Louis work" 33 "Tampa work" 34 "Virginia Beach" 35 "Washington DC work" 36 "NE small work" 37 "NE medium work" 38 "Mid Atl small work" 39 "Mid Atl medium work" 40 "E N Central small work" 41 "E N Central medium work" 42 "W N Central small work" 43 "W N Central medium work" 44 "S Atl small work" 45 "S Atl medium work" 46 "E S Central small work" 47 "E S Central medium work" 48 "W S Central small work" 49 "W S Central medium work" 50 "Mountain small work" 51 "Mountain medium work" 52 "Pacific small work" 53 "Pacific medium work"54 "Atlanta home" 55 "Austin home" 56 "Baltimore home" 57 "Boston home" 58 "Chicago home" 59 "Cincinnati home" 60 "Cleveland home" 61 "Columbus home" 62 "Dallas home" 63 "Denver home" 64 "Detroit home" 65 "Houston home" 66 "Indianapolis home" 67 "Kansas City home" 68 "Knoxville home" 69 "Los Angeles home" 70 "Miami home" 71 "Milwaukee home" 72 "Minneapolis home" 73 "New York home" 74 "Philadelphia home" 75 "Phoenix home" 76 "Pittsburgh home" 77 "Portland home" 78 "Providence"79 "Richmond home" 80 "Riverside home" 81 "Sacramento home" 82 "San Diego home" 83 "San Francisco home" 84 "Seattle home" 85 "St. Louis home" 86 "Tampa home" 87 "Virginia Beach"88 "Washington DC home" 89 "NE small home" 90 "NE medium home" 91 "Mid Atl small home" 92 "Mid Atl medium home" 93 "E N Central small home" 94 "E N Central medium home" 95 "W N Central small home" 96 "W N Central medium home" 97 "S Atl small home" 98 "S Atl medium home" 99 "E S Central small home" 100 "E S Central medium home" 101 "W S Central small home" 102 "W S Central medium home" 103 "Mountain small home" 104 "Mountain medium home" 105 "Pacific small home" 106 "Pacific medium home"
lab val choice106 vlchoice106

generat choice108 = choice54b    if inlf==1
replace choice108 = choice54b+54 if inlf==0
lab def vlchoice108 1 "Atlanta work" 2 "Austin work" 3 "Baltimore work" 4 "Boston work" 5 "Chicago work" 6 "Cincinnati work" 7 "Cleveland work" 8 "Columbus work" 9 "Dallas work" 10 "Denver work" 11 "Detroit work" 12 "Houston work" 13 "Indianapolis work" 14 "Kansas City work" 15 "Knoxville work" 16 "Los Angeles work" 17 "Miami work" 18 "Milwaukee work" 19 "Minneapolis work" 20 "New York work" 21 "Philadelphia work" 22 "Phoenix work" 23 "Pittsburgh work" 24 "Portland work" 25 "Providence"26 "Richmond work" 27 "Riverside work" 28 "Sacramento work" 29 "San Diego work" 30 "San Francisco work" 31 "Seattle work" 32 "St. Louis work" 33 "Tampa work" 34 "Virginia Beach" 35 "Washington DC work" 36 "NE small work" 37 "NE medium work" 38 "Mid Atl small work" 39 "Mid Atl medium work" 40 "E N Central small work" 41 "E N Central medium work" 42 "W N Central small work" 43 "W N Central medium work" 44 "S Atl small work" 45 "S Atl medium work" 46 "E S Central small work" 47 "E S Central medium work" 48 "W S Central small work" 49 "W S Central medium work" 50 "Mountain small work" 51 "Mountain medium work" 52 "Pacific small work" 53 "Pacific medium work" 54 "Alaska/Hawaii work" 55 "Atlanta home" 56 "Austin home" 57 "Baltimore home" 58 "Boston home" 59 "Chicago home" 60 "Cincinnati home" 61 "Cleveland home" 62 "Columbus home" 63 "Dallas home" 64 "Denver home" 65 "Detroit home" 66 "Houston home" 67 "Indianapolis home" 68 "Kansas City home" 69 "Knoxville home" 70 "Los Angeles home" 71 "Miami home" 72 "Milwaukee home" 73 "Minneapolis home" 74 "New York home" 75 "Philadelphia home" 76 "Phoenix home" 77 "Pittsburgh home" 78 "Portland home" 79 "Providence"80 "Richmond home" 81 "Riverside home" 82 "Sacramento home" 83 "San Diego home" 84 "San Francisco home" 85 "Seattle home" 86 "St. Louis home" 87 "Tampa home" 88 "Virginia Beach"89 "Washington DC home" 90 "NE small home" 91 "NE medium home" 92 "Mid Atl small home" 93 "Mid Atl medium home" 94 "E N Central small home" 95 "E N Central medium home" 96 "W N Central small home" 97 "W N Central medium home" 98 "S Atl small home" 99 "S Atl medium home" 100 "E S Central small home" 101 "E S Central medium home" 102 "W S Central small home" 103 "W S Central medium home" 104 "Mountain small home" 105 "Mountain medium home" 106 "Pacific small home" 107 "Pacific medium home" 108 "Alaska/Hawaii home"
lab val choice108 vlchoice108

generat choice110 = choice55    if inlf==1
replace choice110 = choice55+55 if inlf==0
lab def vlchoice110 1 "Atlanta work" 2 "Austin work" 3 "Baltimore work" 4 "Boston work" 5 "Chicago work" 6 "Cincinnati work" 7 "Cleveland work" 8 "Columbus work" 9 "Dallas work" 10 "Denver work" 11 "Detroit work" 12 "Houston work" 13 "Indianapolis work" 14 "Kansas City work" 15 "Knoxville work" 16 "Los Angeles work" 17 "Miami work" 18 "Milwaukee work" 19 "Minneapolis work" 20 "New York work" 21 "Philadelphia work" 22 "Phoenix work" 23 "Pittsburgh work" 24 "Portland work" 25 "Providence"26 "Richmond work" 27 "Riverside work" 28 "Sacramento work" 29 "San Diego work" 30 "San Francisco work" 31 "Seattle work" 32 "St. Louis work" 33 "Tampa work" 34 "Virginia Beach" 35 "Washington DC work" 36 "NE small work" 37 "NE medium work" 38 "Mid Atl small work" 39 "Mid Atl medium work" 40 "E N Central small work" 41 "E N Central medium work" 42 "W N Central small work" 43 "W N Central medium work" 44 "S Atl small work" 45 "S Atl medium work" 46 "E S Central small work" 47 "E S Central medium work" 48 "W S Central small work" 49 "W S Central medium work" 50 "Mountain small work" 51 "Mountain medium work" 52 "Pacific small work" 53 "Pacific medium work" 54 "Alaska work" 55 "Hawaii work" 56 "Atlanta home" 57 "Austin home" 58 "Baltimore home" 59 "Boston home" 60 "Chicago home" 61 "Cincinnati home" 62 "Cleveland home" 63 "Columbus home" 64 "Dallas home" 65 "Denver home" 66 "Detroit home" 67 "Houston home" 68 "Indianapolis home" 69 "Kansas City home" 70 "Knoxville home" 71 "Los Angeles home" 72 "Miami home" 73 "Milwaukee home" 74 "Minneapolis home" 75 "New York home" 76 "Philadelphia home" 77 "Phoenix home" 78 "Pittsburgh home" 79 "Portland home" 80 "Providence"81 "Richmond home" 82 "Riverside home" 83 "Sacramento home" 84 "San Diego home" 85 "San Francisco home" 86 "Seattle home" 87 "St. Louis home" 88 "Tampa home" 89 "Virginia Beach"90 "Washington DC home" 91 "NE small home" 92 "NE medium home" 93 "Mid Atl small home" 94 "Mid Atl medium home" 95 "E N Central small home" 96 "E N Central medium home" 97 "W N Central small home" 98 "W N Central medium home" 99 "S Atl small home" 100 "S Atl medium home" 101 "E S Central small home" 102 "E S Central medium home" 103 "W S Central small home" 104 "W S Central medium home" 105 "Mountain small home" 106 "Mountain medium home" 107 "Pacific small home" 108 "Pacific medium home" 109 "Alaska home" 110 "Hawaii home"
lab val choice110 vlchoice110


*----------------------------
* longitudinal survey weights
*----------------------------
bys ID (timmonth): gen weightlong = weight[1]



*******************************************************
*********** ADD ATTRITION ANALYSIS HERE ***************
*******************************************************
gen anyFlag    = schoolFlag | ageFlag | young | old | noCalDateFlag | gapFlag | misPopFlag | ~in_w1 | flag_in_der~=1 | mi(age) // | ihgcFlag
gen attritFlag = (ageFlag | noCalDateFlag | gapFlag | misPopFlag) & ~young & ~old & ~schoolFlag // & ~ihgcFlag
replace attritFlag = . if young | old | schoolFlag | ~in_w1 | mi(flag_in_ser) | mi(age) // | ihgcFlag
gen agefloor = floor(age)

// see how different attriters are from rest of sample
tab attritFlag if firstper, sum(hgc)
tab attritFlag if firstper, sum(age)
tab agefloor attritFlag if firstper, sum(pop_cat)
tab attritFlag if firstper, sum(exper)
tab attritFlag if firstper, sum(tenure)

// examine initial conditions
preserve
keep if firstper & ~anyFlag
tab pop_cat, sum(hgc)
tab pop_cat, sum(exper)
tab pop_cat, sum(tenure)
restore

// make sure 18-20 year olds in the sample are NOT in school
sum hgc if inrange(age,18,20) & ~schoolFlag & ~ageFlag & ~noCalDateFlag & ~gapFlag & ~misPopFlag & ~young & ~old // & ~ihgcFlag

preserve
  keep if ~anyFlag
  xtsum ID
  xtdescribe, patterns(100)
restore

sum *Flag
sum *Flag if ~young & ~old

replace cbsaname = cntyname if cbsaname=="99999" & regexm(substr(cntyname,-2,.),"[A-Z][A-Z]$")
replace cbsaname = cntyname+stabb if cbsaname=="99999"

compress
//! gunzip -f [REDACTED]sipp2004NHWmaleTEST.dta.gz
save "[REDACTED]sipp2004NHWmaleTEST.dta", replace
! gzip -f [REDACTED]sipp2004NHWmaleTEST.dta

// drop if ageFlag | young | old | noCalDateFlag | gapFlag | misPopFlag // | ihgcFlag
// schoolFlag | ageFlag | young | old | noCalDateFlag | gapFlag | misPopFlag | ~in_w1 | mi(flag_in_ser) | mi(age)
preserve
  xtsum ID
  drop if noCalDateFlag
  xtsum ID
  drop if gapFlag
  xtsum ID
  drop if misPopFlag
  xtsum ID
  drop if young | old
  xtsum ID
  drop if ageFlag
  xtsum ID
  xtsum ID if ~schoolFlag

  compress
  save "[REDACTED]sipp2004NHWmale.dta", replace
  ! gzip -f [REDACTED]sipp2004NHWmale.dta
restore

preserve
  tempfile annual
  sum switch_loc
  keep if inlist(timmonth,4,8,20,32,44)
  xtset ID timmonth
  drop switch_state switch_MSA switch_loc switch_county pop_cat_lag
  bys ID (timmonth): generat byte switch_state  =  state[_n-1] ~=state[_n] & _n>1
  bys ID (timmonth): generat byte switch_county =  switch_state | (county[_n-1]~=county[_n] & state[_n-1]==state[_n] & _n>1)
  bys ID (timmonth): generat byte switch_MSA    = (cbsaname[_n-1]~=cbsaname[_n] & _n>1)
  bys ID (timmonth): generat byte switch_pop_cat=  pop_cat[_n-1]~=pop_cat[_n] & _n>1
  bys ID (timmonth): generat byte switch_house  = mover==2
  bys ID (timmonth): generat pop_cat_lag        =  pop_cat[_n-1]
  generat switch_division = 0
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"CT","ME","MA","NH","RI","VT")                & ~inlist(stabb[_n],"CT","ME","MA","NH","RI","VT")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"NJ","NY","PA")                               & ~inlist(stabb[_n],"NJ","NY","PA")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"IN","IL","MI","OH","WI")                     & ~inlist(stabb[_n],"IN","IL","MI","OH","WI")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"IA","KS","MN","MO","NE","ND","SD")           & ~inlist(stabb[_n],"IA","KS","MN","MO","NE","ND","SD")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"DE","DC","FL","GA","MD","NC","SC","VA","WV") & ~inlist(stabb[_n],"DE","DC","FL","GA","MD","NC","SC","VA","WV")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AL","KY","MS","TN")                          & ~inlist(stabb[_n],"AL","KY","MS","TN")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AR","LA","OK","TX")                          & ~inlist(stabb[_n],"AR","LA","OK","TX")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AZ","CO","ID","NM","MT","UT","NV","WY")      & ~inlist(stabb[_n],"AZ","CO","ID","NM","MT","UT","NV","WY")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AK","CA","HI","OR","WA")                     & ~inlist(stabb[_n],"AK","CA","HI","OR","WA")

  bys ID (timmonth): gen choice27_lag  = choice27[_n-1]
  bys ID (timmonth): gen choice28_lag  = choice28[_n-1]
  bys ID (timmonth): gen choice29_lag  = choice29[_n-1]
  bys ID (timmonth): gen choice48_lag  = choice48[_n-1]
  bys ID (timmonth): gen choice49_lag  = choice49[_n-1]
  bys ID (timmonth): gen choice50_lag  = choice50[_n-1]
  bys ID (timmonth): gen choice53_lag  = choice53[_n-1]
  bys ID (timmonth): gen choice54b_lag = choice54b[_n-1]
  bys ID (timmonth): gen choice55_lag  = choice55[_n-1]
  
  bys ID (timmonth): gen choice54_lag  = choice54[_n-1]
  bys ID (timmonth): gen choice56_lag  = choice56[_n-1]
  bys ID (timmonth): gen choice58_lag  = choice58[_n-1]
  bys ID (timmonth): gen choice96_lag  = choice96[_n-1]
  bys ID (timmonth): gen choice98_lag  = choice98[_n-1]
  bys ID (timmonth): gen choice100_lag = choice100[_n-1]
  bys ID (timmonth): gen choice106_lag = choice106[_n-1]
  bys ID (timmonth): gen choice108_lag = choice108[_n-1]
  bys ID (timmonth): gen choice110_lag = choice110[_n-1]
  bys ID (timmonth): gen empFT_lag     = empFT[_n-1]
  bys ID (timmonth): gen emp_lag       = emp[_n-1]
  bys ID (timmonth): gen inlf_lag      = inlf[_n-1]
  drop if wavemap==1
  xtset ID year
  sum switch_loc
  compress
  save `annual', replace
  save "[REDACTED].dta", replace
  ! gzip -f [REDACTED]sipp2004NHWmaleAnnual.dta
restore

preserve
  tempfile trimesterly
  sum switch_loc
  keep if inlist(timmonth,4,8,12,16,20,24,28,32,36,40,44,48)
  xtset ID wavemap
  drop switch_state switch_MSA switch_loc switch_county pop_cat_lag
  bys ID (timmonth): generat byte switch_state  =  state[_n-1] ~=state[_n] & _n>1
  bys ID (timmonth): generat byte switch_county =  switch_state | (county[_n-1]~=county[_n] & state[_n-1]==state[_n] & _n>1)
  bys ID (timmonth): generat byte switch_MSA    = (cbsaname[_n-1]~=cbsaname[_n] & _n>1)
  bys ID (timmonth): generat byte switch_pop_cat=  pop_cat[_n-1]~=pop_cat[_n] & _n>1
  bys ID (timmonth): generat byte switch_house  = mover==2
  bys ID (timmonth): generat pop_cat_lag        =  pop_cat[_n-1]
  generat switch_division = 0
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"CT","ME","MA","NH","RI","VT")                & ~inlist(stabb[_n],"CT","ME","MA","NH","RI","VT")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"NJ","NY","PA")                               & ~inlist(stabb[_n],"NJ","NY","PA")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"IN","IL","MI","OH","WI")                     & ~inlist(stabb[_n],"IN","IL","MI","OH","WI")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"IA","KS","MN","MO","NE","ND","SD")           & ~inlist(stabb[_n],"IA","KS","MN","MO","NE","ND","SD")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"DE","DC","FL","GA","MD","NC","SC","VA","WV") & ~inlist(stabb[_n],"DE","DC","FL","GA","MD","NC","SC","VA","WV")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AL","KY","MS","TN")                          & ~inlist(stabb[_n],"AL","KY","MS","TN")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AR","LA","OK","TX")                          & ~inlist(stabb[_n],"AR","LA","OK","TX")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AZ","CO","ID","NM","MT","UT","NV","WY")      & ~inlist(stabb[_n-1],"AZ","CO","ID","NM","MT","UT","NV","WY")
  bys ID (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AK","CA","HI","OR","WA")                     & ~inlist(stabb[_n-1],"AK","CA","HI","OR","WA")

  bys ID (timmonth): gen choice27_lag  = choice27[_n-1]
  bys ID (timmonth): gen choice28_lag  = choice28[_n-1]
  bys ID (timmonth): gen choice29_lag  = choice29[_n-1]
  bys ID (timmonth): gen choice48_lag  = choice48[_n-1]
  bys ID (timmonth): gen choice49_lag  = choice49[_n-1]
  bys ID (timmonth): gen choice50_lag  = choice50[_n-1]
  bys ID (timmonth): gen choice53_lag  = choice53[_n-1]
  bys ID (timmonth): gen choice54b_lag = choice54b[_n-1]
  bys ID (timmonth): gen choice55_lag  = choice55[_n-1]
  
  bys ID (timmonth): gen choice54_lag  = choice54[_n-1]
  bys ID (timmonth): gen choice56_lag  = choice56[_n-1]
  bys ID (timmonth): gen choice58_lag  = choice58[_n-1]
  bys ID (timmonth): gen choice96_lag  = choice96[_n-1]
  bys ID (timmonth): gen choice98_lag  = choice98[_n-1]
  bys ID (timmonth): gen choice100_lag = choice100[_n-1]
  bys ID (timmonth): gen choice106_lag = choice106[_n-1]
  bys ID (timmonth): gen choice108_lag = choice108[_n-1]
  bys ID (timmonth): gen choice110_lag = choice110[_n-1]
  bys ID (timmonth): gen empFT_lag     = empFT[_n-1]
  bys ID (timmonth): gen emp_lag       = emp[_n-1]
  bys ID (timmonth): gen inlf_lag      = inlf[_n-1]
  
  drop exper experAlt
  bys ID (timmonth): replace empFT_lag = 0 if timmonth==4
  gen exper_temp = (1/3)*lag_empFT
  bys ID (timmonth): gen exper_s = sum(exper_temp)
  gen exper     = exper0      + exper_s
  gen experAlt  = exper0_alt  + exper_s
  drop exper_temp exper_s
  
  drop if timmonth==4
  xtset ID wavemap
  sum switch_loc
  compress
  save `trimesterly', replace
  save "[REDACTED]sipp2004NHWmaleTrimester.dta", replace
  ! gzip -f [REDACTED]sipp2004NHWmaleTrimester.dta
restore

*******************************************************
*********** EXPORT (ANNUAL) WIDE PANEL TO MATLAB ******
*******************************************************
use `annual', clear
gen constant = 1
drop anyFlag
gen  anyFlag   = schoolFlag | ageFlag | young | old | noCalDateFlag | gapFlag | misPopFlag | ~in_w1 | flag_in_der~=1 | mi(age)
local outcomes choice* pop_cat pop_cat_lag lnWageHr lnWageHrJ lnWageHrJb lnearnfinal lnearnfinalJ lnearnfinalJb inlf inlf_lag empFT empFT_lag emp emp_lag rhrsweek anyFlag earnflag
local Svars hgc educlevel
local Gvars urate
local Evars age exper experFTorPT experAlt

preserve
  gen earnflag   = (empFT & inrange(lnearnfinal,6,10) & ~anyFlag)
  gen earnflagJ  = (empFT & inrange(lnearnfinal,6,10) & ~anyFlag)
  gen earnflagJb = (empFT & inrange(lnearnfinal,6,10) & ~anyFlag)
  * keep if ~anyFlag
  * drop if schoolFlag
  sum lnWage* lnearnfinal* constant pop_cat hgc educlevel age exper tenure
  replace lnWageHr      = -999 if ~earnflag
  replace lnWageHrJ     = -999 if ~earnflagJ
  replace lnWageHrJb    = -999 if ~earnflagJb
  replace lnearnfinal   = -999 if ~earnflag
  replace lnearnfinalJ  = -999 if ~earnflagJ
  replace lnearnfinalJb = -999 if ~earnflagJb
  sum lnWage* lnearn* constant pop_cat hgc educlevel age exper tenure
  keep ID timmonth constant weightlong `outcomes' `Evars' `Gvars' `Svars'
  reshape wide `outcomes' `Evars' `Gvars' `Svars', i(ID) j(timmonth)
  order ID constant weightlong
  outsheet using "[REDACTED]sipp2004NHWmale_matlab_wide_annual.csv", comma replace nolabel
  ! gzip -f [REDACTED]sipp2004NHWmale_matlab_wide_annual.csv
restore

log close
exit

