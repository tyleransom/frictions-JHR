version 14.1
clear all
set more off
capture log close
log using an_sipp_nhww.log, replace

local pufcntrl ssuid swave srotaton srefmon spanel eentaid epppnum lgtkey epopstat
local pufpanel aorigin eorigin // ebornus abornus
local pufmnth rhtype ehhnumpp tage ems eenrlm rwkesr1 rwkesr2 rwkesr3 rwkesr4 rwkesr5 rmwkwjb rmwksab ersend1 ersend2 ersnowrk eabre elayoff wpfinwgt tpmsum1 tpmsum2 eemploc1 eemploc2 tempsiz1 tempsiz2 estlemp1 estlemp2 tsjdate1 tsjdate2 tejdate1 tejdate2 eclwrk1 eclwrk2 swave rmhrswk rmesr ?jbhrs? tmlmsum
local pufwave esex erace rhcalyr rhcalmn eppintvw eeducate renrlma eenlevel ejbind1 tjbocc1 ejbind2 tjbocc2 eeno1 eeno2 ebflag ecflag ejobcntr epayhr1 apayhr1 epayhr2 apayhr2 elkwrk alkwrk epropb1 epropb2 eptwrk aptwrk eptresn aptresn eoutcome renroll // rged ehrsall tmovrflg 
local iuwavegeo tfipsst
local pufimpwave arace aeducate aenlevel ajbind1 ajbocc1 ajbind2 ajbocc2 // aged ahrsall 
local pufimpmnth aenrlm arsend1 arsend2 arsnowrk aabre alayoff apmsum1 apmsum2 aemploc1 aemploc2 aempsiz1 aempsiz2 astlemp1 astlemp2 asjdate1 asjdate2 aejdate1 aejdate2 aclwrk1 aclwrk2

local myvars96 `pufcntrl' `pufpanel' `pufmnth' `pufwave' `iuwavegeo' `pufimpwave' `pufimpmnth'
local myvars `pufcntrl' `pufpanel' `pufmnth' `pufwave' `iuwavegeo' `pufimpwave' `pufimpmnth' ebornus abornus rged ehrsall tmovrflg aged ahrsall

local tm1vars spanel ssuid swave srotaton eentaid epppnum tage esex epopstat erace eorigin tmakmnyr tprvjbyr tlstwrky tfrmryr etimeoff eoff6mtn amakmnyr aprvjbyr alstwrky afrmryr atimeoff aoff6mtn
local tm2vars spanel ssuid swave srotaton eentaid epppnum tmovyryr toutinyr tmovest egedtm ebachfld tprstate tbrstate eprevres amovyryr aoutinyr amovest agedtm abachfld aprstate abrstate aprevres
local tm1varsNew yrenterlf yrprevlf yrexitlf yrstartpvjob mosoutlf exitlfcaregiver 
local tm2varsNew yr_move_this_house yr_move_pvs_house yr_move_this_state gedever college_major prev_state_country birth_state_country prev_home_loc

* Two changes to the file:
*  1: drop all obs for which "hours vary"; this will only impact later 
*     years (2003-2013), but it should help with the drop off
*     - This is the first thing done now after saving ../Data/Temp/full`foo'_data.dta
*     - Requires 8 lines of code or so for each survey
*  2: change hour requirement from 35 to 30; this will help with all years
*     - global find and replace of: 
*       'hours_job1_month`k' >= 35' to 'hours_job1_month`k' >= 30'
*       'hours_job2_month`k' >= 35' to 'hours_job2_month`k' >= 30'

* * 1996
* forvalues k=1(1)12 {
*     append using "../Data/Temp/sipp96_core`k'.dta", keep(`myvars96')
*     keep if srefmon==4 & esex==1 & erace==1 & inrange(tage,18,55) & !inrange(eorigin,20,28)
*     tab1 swave srefmon
* }

* * 2001
* forvalues k=1(1)9 {
*     append using "../Data/Temp/sipp01_core`k'.dta", keep(`myvars96')
*     keep if srefmon==4 & esex==1 & erace==1 & inrange(tage,18,55) & !inrange(eorigin,20,28)
*     tab1 swave srefmon
* }

tempfile tm1 tm2
append using "../Data/Temp/sipp04wave1.dta", keep(`tm1vars')
append using "../Data/Temp/sipp08wave1.dta", keep(`tm1vars')
ren tmakmnyr yrenterlf
ren tprvjbyr yrprevlf
ren tlstwrky yrexitlf
ren tfrmryr  yrstartpvjob
ren etimeoff mosoutlf
ren eoff6mtn exitlfcaregiver
ren amakmnyr qyrenterlf
ren aprvjbyr qyrprevlf
ren alstwrky qyrexitlf
ren afrmryr  qyrstartpvjob
ren atimeoff qmosoutlf
ren aoff6mtn qexitlfcaregiver
save `tm1', replace
clear

append using "../Data/Temp/sipp04wave2.dta", keep(`tm2vars')
append using "../Data/Temp/sipp08wave2.dta", keep(`tm2vars')
ren tmovyryr yr_move_this_house
ren toutinyr yr_move_pvs_house
ren tmovest yr_move_this_state
ren egedtm gedever
ren ebachfld college_major
ren tprstate prev_state_country
ren tbrstate birth_state_country
ren eprevres prev_home_loc
ren amovyryr qyr_move_this_house
ren aoutinyr qyr_move_pvs_house
ren amovest qyr_move_this_state
ren agedtm qgedever
ren abachfld qcollege_major
ren aprstate qprev_state_country
ren abrstate qbirth_state_country
save `tm2', replace
clear

* 2004
forvalues k=1(1)12 {
    append using "../Data/Temp/sipp04_core`k'.dta", keep(`myvars')
    keep if srefmon==4 & esex==2 & erace==1 & eeducate<44 & !inrange(eorigin,1,1)
    tab1 swave srefmon
}

* 2008
forvalues k=1(1)16 {
    append using "../Data/Temp/sipp08_core`k'.dta", keep(`myvars')
    keep if srefmon==4 & esex==2 & erace==1 & eeducate<44 & !inrange(eorigin,1,1)
    tab1 swave srefmon
}

merge m:1 spanel ssuid eentaid epppnum using `tm1', nogen keep(match master)
gen missW1 = mi(yrenterlf)
merge m:1 spanel ssuid eentaid epppnum using `tm2', nogen keep(match master)
gen missW2 = mi(birth_state_country)

ds `myvars' `tm1varsNew' `tm2varsNew', has(type numeric)
local mynumericvars `r(varlist)'
local qvars a* q*
recode `mynumericvars' (-1 = .v)
recode `qvars' (2 3 4 = 1)

egen ID = group(spanel ssuid eentaid epppnum)
gen timmonth = swave*4-(4-srefmon+1)+1

*------------------------------------------------------------------
* generate mover flag
*------------------------------------------------------------------
bys ID (timmonth): egen moverMax1 = max(tmovrflg) if inlist(timmonth,12,16,20)
bys ID (timmonth): egen moverMax2 = max(tmovrflg) if inlist(timmonth,24,28,32)
bys ID (timmonth): egen moverMax3 = max(tmovrflg) if inlist(timmonth,36,40,44)
bys ID (timmonth): egen moverMax4 = max(tmovrflg) if inlist(timmonth,48,52,56)
bys ID (timmonth): egen everMovedA = max(tmovrflg) if inrange(tmovrflg,0,4)
replace everMovedA = 1 if mi(everMovedA)
gen everMoved = everMovedA>1
clonevar tmovrflgA = tmovrflg
replace  tmovrflgA = moverMax1 if inlist(timmonth,12,16,20)
replace  tmovrflgA = moverMax2 if inlist(timmonth,24,28,32)
replace  tmovrflgA = moverMax3 if inlist(timmonth,36,40,44)
replace  tmovrflgA = moverMax4 if inlist(timmonth,48,52,56)
l ID timmonth everMovedA tmovrflg tmovrflgA moverMax? if everMoved in 1/100000, sepby(ID)

*------------------------------------------------------------------
* Keep only one wave per year and rename variables
*------------------------------------------------------------------
keep if inlist(timmonth,4,8,20,32,44,56,68)
xtset ID timmonth
xtsum ID if swave==1

ren rhtype hhtype
ren ehhnumpp hhsize
ren tage age
ren ems marst
ren eenrlm enrlsch
ren rwkesr1 empstatwk1
ren rwkesr2 empstatwk2
ren rwkesr3 empstatwk3
ren rwkesr4 empstatwk4
ren rwkesr5 empstatwk5
ren rmwkwjb weekswjob
ren rmwksab weeksabjob
ren ersend1 rsn_stop_j1
ren ersend2 rsn_stop_j2
ren ersnowrk rsn_no_work
ren eabre rsn_absent
ren elayoff on_layoff
ren wpfinwgt weight
ren tpmsum1 monthearnj1
ren tpmsum2 monthearnj2
ren eemploc1 multlocj1
ren eemploc2 multlocj2
ren tempsiz1 firmsizej1
ren tempsiz2 firmsizej2
ren estlemp1 stillatj1
ren estlemp2 stillatj2
ren tsjdate1 sdatej1
ren tsjdate2 sdatej2
ren tejdate1 edatej1
ren tejdate2 edatej2
ren eclwrk1 classj1
ren eclwrk2 classj2
ren swave wavemap
ren rhcalyr year
ren rhcalmn calmonth
ren rmhrswk rhrsweek
ren esex sex
ren erace race
ren eppintvw intstatus
ren eeducate hgc_cat
ren rged ged
ren renrlma enrsch
ren eenlevel enrwhich
ren ehrsall hrswork
ren ejbind1 ind1
ren tjbocc1 occ1
ren ejbind2 ind2
ren tjbocc2 occ2
ren eeno1 job1id
ren eeno2 job2id
ren ebflag qjobflag
ren ecflag contingent
ren ejobcntr numjobs
ren tmovrflgA mover
ren epayhr1 hourly1
ren apayhr1 qhourly1
ren epayhr2 hourly2
ren apayhr2 qhourly2
ren elkwrk look_work
ren alkwrk qlook_work
ren epropb1 prop1
ren epropb2 prop2
ren eptwrk workedpt1plusweeks
ren aptwrk qworkedpt1plusweeks
ren eptresn reasonworkpt
ren aptresn qreasonworkpt
ren arace qrace
ren aeducate qhgc
ren aged qged
ren aenlevel qenrwhich
ren ahrsall qhrswork
ren ajbind1 qind1
ren ajbocc1 qocc1
ren ajbind2 qind2
ren ajbocc2 qocc2
ren aenrlm qenrlsch
ren arsend1 qrsn_stop_j1
ren arsend2 qrsn_stop_j2
ren arsnowrk qrsn_no_work
ren aabre qrsn_absent
ren alayoff qon_layoff
ren apmsum1 qmonthearnj1
ren apmsum2 qmonthearnj2
ren aemploc1 qmultlocj1
ren aemploc2 qmultlocj2
ren aempsiz1 qfirmsizej1
ren aempsiz2 qfirmsizej2
ren astlemp1 qstillatj1
ren astlemp2 qstillatj2
ren asjdate1 qsdatej1
ren asjdate2 qsdatej2
ren aejdate1 qedatej1
ren aejdate2 qedatej2
ren aclwrk1 qclassj1
ren aclwrk2 qclassj2
ren tfipsst state

bys ID (timmonth): gen  agebeg = age[1] // generate age at beginning of survey

*------------------------------------------------------------------
* get rid of people not present in wave 1
*------------------------------------------------------------------
* fillin ID wave
* gen notInW1 = wave==1 & _fillin==1
* bys ID (timmonth): egen noW1 = mean(notInW1)
* gen noW1flag = noW1>0
* drop if _fillin==1

clonevar noW1flag = missW1

xtsum ID if wave==1
xtsum ID if wave==1 & !noW1flag
xtsum ID
xtsum ID if !noW1flag

drop if noW1flag
xtsum ID if wave==1

*gen noW1flagA  = inrange(eoutcome,208,.) & wave==1
*bys ID (wave): egen noW1flagB   = mean(noW1flagA)
*gen noW1flag   = noW1flagB>0
*l ID noW1flag* if inlist(ID,1,10,100,1000,10000), sepby(ID)
*
*xtsum ID if wave==1
*xtsum ID if wave==1 & !noW1flag
*xtsum ID 
*xtsum ID if !noW1flag

drop if noW1flag
xtsum ID if wave==1

*------------------------------------------------------------------
* get rid of people in school
*------------------------------------------------------------------
recode enrsch (2 .v = 0)
* generate variable for last time enrolled in school
gen last_enr_schA = timmonth*enrsch
bys ID (timmonth): egen last_enr_sch  = max(last_enr_schA)
bys ID (timmonth): egen ever_enr_schA = mean(last_enr_schA)
* generate schoolFlag variable to mark that a person was in school
gen schoolFlagOld = timmonth<=last_enr_sch
gen schoolFlag    = ever_enr_schA>0

drop if schoolFlag
xtsum ID if wave==1

*------------------------------------------------------------------
* only keep people in prime working ages
*------------------------------------------------------------------
gen young = agebeg<18
gen old   = agebeg>55

xtsum ID if wave==1   & inrange(agebeg,18,55)
xtsum ID if timmonth>5 & inrange(agebeg,18,55)

gen weirdage = age<18 | age>60
bys ID (timmonth): egen everWeirdAge = max(weirdage)

drop if !inrange(agebeg,18,55) | everWeirdAge

*------------------------------------------------------------------
* Drop immigrants
*------------------------------------------------------------------
drop if birth_state_country>56

*------------------------------------------------------------------
* Drop observations after missed interview
*------------------------------------------------------------------
bys ID (timmonth): gen gap = wave[_n]-wave[_n-1]>3 if _n>1
gen t_gap = timmonth*gap
recode t_gap (0 = .)
bys ID (timmonth): egen first_gap = min(t_gap)
gen post_gap = (timmonth>=first_gap)
drop gap t_gap first_gap
gen gapFlag = post_gap

drop if gapFlag

*----------------------------------------------------
* experience
*----------------------------------------------------
replace yrexitlf = 2004 if (mi(yrexitlf) | yrexitlf==0) & spanel==2004
replace yrexitlf = 2008 if (mi(yrexitlf) | yrexitlf==0) & spanel==2008
replace yrenterlf = 2004 if (mi(yrenterlf) | yrenterlf==0) & spanel==2004
replace yrenterlf = 2008 if (mi(yrenterlf) | yrenterlf==0) & spanel==2008
replace mosoutlf = 0    if  mi(mosoutlf)
gen exper0temp = yrexitlf-yrenterlf-(mosoutlf/12) if wave==1
bys ID (timmonth): egen exper0 = mean(exper0temp)
drop exper0temp
replace exper0 = 0 if exper0<0

*----------------------------------------------------
* generate labor force participation/outcomes
*----------------------------------------------------
gen selfEmp =  inlist(prop1,1,2) | inlist(prop2,1,2) | numjobs==0 | qjobflag==1 | contingent==1
gen wantPT  =  reasonworkpt==2
gen inlf    = (inrange(rhrsweek,1,6) | (rhrsweek==0 & look_work==1)) & ~selfEmp & ~wantPT
gen empFT   =  inrange(rhrsweek,1,1) & inlf & inlist(reasonworkpt,3,4,10,.v) 
gen emp     =  inrange(rhrsweek,1,2) & inlf & inlist(reasonworkpt,3,4,10,.v) 
gen unemp   =  inlf & ~empFT
bys ID (timmonth): gen empFT_prev = empFT[_n-1]
bys ID (timmonth): gen  inlf_prev =  inlf[_n-1]
bys ID (timmonth): gen state_prev = state[_n-1]

bys ID (timmonth): gen exper = exper0+sum(empFT_prev)
replace exper = 0 if exper<0

l ID timmonth empFT empFT_prev yrexitlf yrenterlf mosoutlf exper0 exper age in 1/50, sepby(ID)
l ID timmonth empFT empFT_prev yrexitlf yrenterlf mosoutlf exper0 exper age if inlist(ID,64,159,533,536,542,585,765), sepby(ID)

drop if wave==1

*----------------------------------------------------
* State unemployment rates
*----------------------------------------------------
tempfile countyurate stateurate
preserve
    use ../../Input/CountyUnemp/county_unemp_annual.dta, clear
    egen cid = group(state county)
    xtset cid year
    bys cid (year): gen lurate = l.urate
    bys cid (year): gen llf    = l.labor_force
    bys cid (year): gen lunemp = l.unemployed
    ren state statefip
    ren county countyfips
    keep year statefip countyfips urate lurate labor_force llf unemployed lunemp
    keep if inrange(year,2001,2016)
    mdesc _all
    l in 1/10, sep(0)
    save `countyurate', replace
restore

preserve
    use `countyurate', clear
    drop lurate urate
    egen ltotlf = total(llf)   , by(year statefip)
    egen ltotun = total(lunemp), by(year statefip)
    gen  lurate = ltotun/ltotlf*100
    egen totlf = total(labor_force), by(year statefip)
    egen totun = total(unemployed) , by(year statefip)
    gen  urate = totun/totlf*100
    bys year statefip (countyfips): gen firstObs = _n==1
    keep if firstObs
    drop firstObs countyfips
    ren statefip state
    l, sepby(year)
    isid state year
    save `stateurate', replace
restore

merge m:1 state year using `stateurate', keep(match master) nogen

*----------------------------------------------------
* Summary stats
*----------------------------------------------------
gen mover_dummy = inlist(mover,3,4)
sum empFT if empFT_prev==0
sum empFT if empFT_prev==1
sum empFT if empFT_prev==0 & mover_dummy==1
sum empFT if empFT_prev==0 & mover_dummy==0
sum empFT if empFT_prev==1 & mover_dummy==1
sum empFT if empFT_prev==1 & mover_dummy==0

*----------------------------------------------------
* Logits -- lagged urate
*----------------------------------------------------
logit empFT c.exper##c.exper        mover_dummy if empFT_prev==1
logit empFT c.exper##c.exper        mover_dummy if empFT_prev==0
logit empFT c.exper##c.exper lurate mover_dummy if empFT_prev==1
logit empFT c.exper##c.exper lurate mover_dummy if empFT_prev==0

*----------------------------------------------------
* Logits -- current urate
*----------------------------------------------------
logit empFT c.exper##c.exper  urate mover_dummy if empFT_prev==1
logit empFT c.exper##c.exper  urate mover_dummy if empFT_prev==0

*----------------------------------------------------
* LPMs -- lagged urate
*----------------------------------------------------
reg empFT c.exper##c.exper        mover_dummy if empFT_prev==1
reg empFT c.exper##c.exper        mover_dummy if empFT_prev==0
reg empFT c.exper##c.exper lurate mover_dummy if empFT_prev==1
reg empFT c.exper##c.exper lurate mover_dummy if empFT_prev==0

reg empFT c.exper##c.exper lurate mover_dummy i.state##i.year if empFT_prev==1
reg empFT c.exper##c.exper lurate mover_dummy i.state##i.year if empFT_prev==0


*----------------------------------------------------
* LPMs -- lagged urate
*----------------------------------------------------
clonevar moverd = mover_dummy
replace  moverd = 3*mover_dummy if wave==2
tab year empFT_prev, sum(moverd) mean nofreq


*----------------------------------------------------
* Graphs
*----------------------------------------------------
generat MRcity  = mover==3
generat MRstate = mover==4
replace MRcity  = 3*(mover==3) if wave==2
replace MRstate = 3*(mover==4) if wave==2

* Graph mobility rates by distance, time, and prev employment status (cond'l non-college grad)
preserve
	collapse MR* , by(empFT_prev year)
	set obs `=_N+1'
	replace year = 2009.5 if _n==_N
	set obs `=_N+1'
	replace year = 2007.9 if _n==_N
	sort year empFT_prev
	local cheight1 5
	gen c1 = `cheight1' if inrange(year,2007.8,2009.6)
	local cheight 5
	gen c2 = `cheight' if inrange(year,2007.8,2009.6)
	l, sepby(year)
	generat MRSameState    = 100*MRcity
	generat MRDiffState    = 100*MRstate   

	graph set eps fontface "Palatino Linotype"
	graph twoway ///(area c2             year if inrange(year,2004,2013)                ,              color(gs15)  )  /// recession
				 (line MRSameState    year if inrange(year,2004,2013) & empFT_prev==1, lpattern(.) lcolor(black) )  /// between counties
				 (line MRDiffState    year if inrange(year,2004,2013) & empFT_prev==1, lpattern(-) lcolor(black) ), /// between states
				 legend(label(1 "Recession") label(2 "Within state") label(3 "Between states") cols(2)) ytitle("Annual Migration Rate (%)") xtitle("Year") ylabel(0(.5)`cheight') xlabel(2004(1)2013) graphregion(color(white))
	graph export mobility_nhww_employed.eps, replace

	graph twoway ///(area c2             year if inrange(year,2004,2013)                ,              color(gs15)  )  /// recession
				 (line MRSameState    year if inrange(year,2004,2013) & empFT_prev==0, lpattern(.) lcolor(black) )  /// between counties
				 (line MRDiffState    year if inrange(year,2004,2013) & empFT_prev==0, lpattern(-) lcolor(black) ), /// between states
				 legend(label(1 "Recession") label(2 "Within state") label(3 "Between states") cols(2)) ytitle("Annual Migration Rate (%)") xtitle("Year") ylabel(0(.5)`cheight') xlabel(2004(1)2013) graphregion(color(white))
	graph export mobility_nhww_nonemployed.eps, replace
restore

* Graph mobility rates by distance, time, and curr employment status (cond'l non-college grad)
preserve
	collapse MR* , by(empFT year)
	set obs `=_N+1'
	replace year = 2009.5 if _n==_N
	set obs `=_N+1'
	replace year = 2007.9 if _n==_N
	sort year empFT
	local cheight1 3.5
	gen c1 = `cheight1' if inrange(year,2007.8,2009.6)
	local cheight 3.5
	gen c2 = `cheight' if inrange(year,2007.8,2009.6)
	l, sepby(year)
	generat MRSameState    = 100*MRcity
	generat MRDiffState    = 100*MRstate   

	graph set eps fontface "Palatino Linotype"
	graph twoway (area c2             year if inrange(year,2004,2013)           ,              color(gs15)  )  /// recession
				 (line MRSameState    year if inrange(year,2004,2013) & empFT==1, lpattern(.) lcolor(black) )  /// between counties
				 (line MRDiffState    year if inrange(year,2004,2013) & empFT==1, lpattern(-) lcolor(black) ), /// between states
				 legend(label(1 "Recession") label(2 "Within state") label(3 "Between states") cols(2)) ytitle("Annual Migration Rate (%)") xtitle("Year") ylabel(0(.5)`cheight') xlabel(2004(1)2013) graphregion(color(white))
	graph export mobility_nhww_curr_employed.eps, replace

	graph twoway (area c2             year if inrange(year,2004,2013)           ,              color(gs15)  )  /// recession
				 (line MRSameState    year if inrange(year,2004,2013) & empFT==0, lpattern(.) lcolor(black) )  /// between counties
				 (line MRDiffState    year if inrange(year,2004,2013) & empFT==0, lpattern(-) lcolor(black) ), /// between states
				 legend(label(1 "Recession") label(2 "Within state") label(3 "Between states") cols(2)) ytitle("Annual Migration Rate (%)") xtitle("Year") ylabel(0(.5)`cheight') xlabel(2004(1)2013) graphregion(color(white))
	graph export mobility_nhww_curr_nonemployed.eps, replace
restore



log close
 
