clear all
set maxvar 32000
version 13.0
capture cd "[REDACTED]"
capture log close
set more off
log using panel_stacker.log, replace

tempfile p2008
	! gunzip -f [REDACTED]sipp2008NHWmaleTEST.dta.gz
	use "[REDACTED]sipp2008NHWmaleTEST.dta", clear
	! gzip -f [REDACTED]sipp2008NHWmaleTEST.dta
	capture confirm variable SPANEL
	if !_rc assert SPANEL==2008
	else gen SPANEL=2008
save `p2008', replace

! gunzip -f [REDACTED]sipp2004NHWmaleTEST.dta.gz
use "[REDACTED]sipp2004NHWmaleTEST.dta", clear
! gzip -f [REDACTED]sipp2004NHWmaleTEST.dta

append using `p2008'
qui sum ID if SPANEL==2004
local max04 = `r(max)'
qui sum ID if SPANEL==2008
local max08 = `r(max)'
local maxoverall = max(`max04',`max08')
local numdigits = length("`maxoverall'")
gen double IDtilde = SPANEL*1e`numdigits'+ID
tostring IDtilde, gen(IDstring)
codebook IDtilde IDstring
egen puid_new = concat(SPANEL puid)

// ! gunzip -f [REDACTED]sippCombinedNHWmaleTEST.dta.gz
save "[REDACTED]sippCombinedNHWmaleTEST.dta", replace
! gzip -f [REDACTED]sippCombinedNHWmaleTEST.dta

generat birthDiv = 0
replace birthDiv = 1 if inlist(birth_state_country,9, 23, 25, 33, 44, 50)
replace birthDiv = 2 if inlist(birth_state_country,34, 36, 42)
replace birthDiv = 3 if inlist(birth_state_country,18, 17, 26, 39, 55)
replace birthDiv = 4 if inlist(birth_state_country,19, 20, 27, 29, 31, 38, 46)
replace birthDiv = 5 if inlist(birth_state_country,10, 11, 12, 13, 24, 37, 45, 51, 54)
replace birthDiv = 6 if inlist(birth_state_country,1, 21, 28, 47)
replace birthDiv = 7 if inlist(birth_state_country,5, 22, 40, 48)
replace birthDiv = 8 if inlist(birth_state_country,4, 8, 16, 35, 30, 49, 32, 56)
replace birthDiv = 9 if inlist(birth_state_country,6, 41, 53, 2, 15)
forvalues x=1/110 {
	qui gen born_here_state`x'_ = 0
	qui gen born_here_div`x'_ = 0
}

qui replace born_here_state1_   = 1 if inlist(birth_state_country,13)
qui replace born_here_state2_   = 1 if inlist(birth_state_country,48)
qui replace born_here_state3_   = 1 if inlist(birth_state_country,24)
qui replace born_here_state4_   = 1 if inlist(birth_state_country,25,33)
qui replace born_here_state5_   = 1 if inlist(birth_state_country,17,18,55)
qui replace born_here_state6_   = 1 if inlist(birth_state_country,18,21,39)
qui replace born_here_state7_   = 1 if inlist(birth_state_country,39)
qui replace born_here_state8_   = 1 if inlist(birth_state_country,39)
qui replace born_here_state9_   = 1 if inlist(birth_state_country,48)
qui replace born_here_state10_  = 1 if inlist(birth_state_country,8)
qui replace born_here_state11_  = 1 if inlist(birth_state_country,26)
qui replace born_here_state12_  = 1 if inlist(birth_state_country,48)
qui replace born_here_state13_  = 1 if inlist(birth_state_country,18)
qui replace born_here_state14_  = 1 if inlist(birth_state_country,29,20)
qui replace born_here_state15_  = 1 if inlist(birth_state_country,47)
qui replace born_here_state16_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state17_  = 1 if inlist(birth_state_country,12)
qui replace born_here_state18_  = 1 if inlist(birth_state_country,55)
qui replace born_here_state19_  = 1 if inlist(birth_state_country,27,55)
qui replace born_here_state20_  = 1 if inlist(birth_state_country,36,9,34,42)
qui replace born_here_state21_  = 1 if inlist(birth_state_country,42,34,10)
qui replace born_here_state22_  = 1 if inlist(birth_state_country,4)
qui replace born_here_state23_  = 1 if inlist(birth_state_country,42)
qui replace born_here_state24_  = 1 if inlist(birth_state_country,41,53)
qui replace born_here_state25_  = 1 if inlist(birth_state_country,44,25)
qui replace born_here_state26_  = 1 if inlist(birth_state_country,51)
qui replace born_here_state27_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state28_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state29_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state30_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state31_  = 1 if inlist(birth_state_country,53)
qui replace born_here_state32_  = 1 if inlist(birth_state_country,29,17)
qui replace born_here_state33_  = 1 if inlist(birth_state_country,12)
qui replace born_here_state34_  = 1 if inlist(birth_state_country,51)
qui replace born_here_state35_  = 1 if inlist(birth_state_country,51,11,24)
qui replace born_here_state36_  = 1 if inlist(birth_state_country,9,23,25,33,44,50)
qui replace born_here_state37_  = 1 if inlist(birth_state_country,9,23,25,33,44,50)
qui replace born_here_state38_  = 1 if inlist(birth_state_country,34,36,42)
qui replace born_here_state39_  = 1 if inlist(birth_state_country,34,36,42)
qui replace born_here_state40_  = 1 if inlist(birth_state_country,18,17,26,39,55)
qui replace born_here_state41_  = 1 if inlist(birth_state_country,18,17,26,39,55)
qui replace born_here_state42_  = 1 if inlist(birth_state_country,19,20,27,29,31,38,46)
qui replace born_here_state43_  = 1 if inlist(birth_state_country,19,20,27,29,31,38,46)
qui replace born_here_state44_  = 1 if inlist(birth_state_country,10,11,12,13,24,37,45,51,54)
qui replace born_here_state45_  = 1 if inlist(birth_state_country,10,11,12,13,24,37,45,51,54)
qui replace born_here_state46_  = 1 if inlist(birth_state_country,1,21,28,47)
qui replace born_here_state47_  = 1 if inlist(birth_state_country,1,21,28,47)
qui replace born_here_state48_  = 1 if inlist(birth_state_country,5,22,40,48)
qui replace born_here_state49_  = 1 if inlist(birth_state_country,5,22,40,48)
qui replace born_here_state50_  = 1 if inlist(birth_state_country,4,8,16,35,30,49,32,56)
qui replace born_here_state51_  = 1 if inlist(birth_state_country,4,8,16,35,30,49,32,56)
qui replace born_here_state52_  = 1 if inlist(birth_state_country,6,41,53)
qui replace born_here_state53_  = 1 if inlist(birth_state_country,6,41,53)
qui replace born_here_state54_  = 1 if inlist(birth_state_country,2)
qui replace born_here_state55_  = 1 if inlist(birth_state_country,15)
qui replace born_here_state56_  = 1 if inlist(birth_state_country,13)
qui replace born_here_state57_  = 1 if inlist(birth_state_country,48)
qui replace born_here_state58_  = 1 if inlist(birth_state_country,24)
qui replace born_here_state59_  = 1 if inlist(birth_state_country,25,33)
qui replace born_here_state60_  = 1 if inlist(birth_state_country,17,18,55)
qui replace born_here_state61_  = 1 if inlist(birth_state_country,18,21,39)
qui replace born_here_state62_  = 1 if inlist(birth_state_country,39)
qui replace born_here_state63_  = 1 if inlist(birth_state_country,39)
qui replace born_here_state64_  = 1 if inlist(birth_state_country,48)
qui replace born_here_state65_  = 1 if inlist(birth_state_country,8)
qui replace born_here_state66_  = 1 if inlist(birth_state_country,26)
qui replace born_here_state67_  = 1 if inlist(birth_state_country,48)
qui replace born_here_state68_  = 1 if inlist(birth_state_country,18)
qui replace born_here_state69_  = 1 if inlist(birth_state_country,29,20)
qui replace born_here_state70_  = 1 if inlist(birth_state_country,47)
qui replace born_here_state71_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state72_  = 1 if inlist(birth_state_country,12)
qui replace born_here_state73_  = 1 if inlist(birth_state_country,55)
qui replace born_here_state74_  = 1 if inlist(birth_state_country,27,55)
qui replace born_here_state75_  = 1 if inlist(birth_state_country,36,9,34,42)
qui replace born_here_state76_  = 1 if inlist(birth_state_country,42,34,10)
qui replace born_here_state77_  = 1 if inlist(birth_state_country,4)
qui replace born_here_state78_  = 1 if inlist(birth_state_country,42)
qui replace born_here_state79_  = 1 if inlist(birth_state_country,41,53)
qui replace born_here_state80_  = 1 if inlist(birth_state_country,44,25)
qui replace born_here_state81_  = 1 if inlist(birth_state_country,51)
qui replace born_here_state82_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state83_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state84_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state85_  = 1 if inlist(birth_state_country,6)
qui replace born_here_state86_  = 1 if inlist(birth_state_country,53)
qui replace born_here_state87_  = 1 if inlist(birth_state_country,29,17)
qui replace born_here_state88_  = 1 if inlist(birth_state_country,12)
qui replace born_here_state89_  = 1 if inlist(birth_state_country,51)
qui replace born_here_state90_  = 1 if inlist(birth_state_country,51,11,24)
qui replace born_here_state91_  = 1 if inlist(birth_state_country,9,23,25,33,44,50)
qui replace born_here_state92_  = 1 if inlist(birth_state_country,9,23,25,33,44,50)
qui replace born_here_state93_  = 1 if inlist(birth_state_country,34,36,42)
qui replace born_here_state94_  = 1 if inlist(birth_state_country,34,36,42)
qui replace born_here_state95_  = 1 if inlist(birth_state_country,18,17,26,39,55)
qui replace born_here_state96_  = 1 if inlist(birth_state_country,18,17,26,39,55)
qui replace born_here_state97_  = 1 if inlist(birth_state_country,19,20,27,29,31,38,46)
qui replace born_here_state98_  = 1 if inlist(birth_state_country,19,20,27,29,31,38,46)
qui replace born_here_state99_  = 1 if inlist(birth_state_country,10,11,12,13,24,37,45,51,54)
qui replace born_here_state100_ = 1 if inlist(birth_state_country,10,11,12,13,24,37,45,51,54)
qui replace born_here_state101_ = 1 if inlist(birth_state_country,1,21,28,47)
qui replace born_here_state102_ = 1 if inlist(birth_state_country,1,21,28,47)
qui replace born_here_state103_ = 1 if inlist(birth_state_country,5,22,40,48)
qui replace born_here_state104_ = 1 if inlist(birth_state_country,5,22,40,48)
qui replace born_here_state105_ = 1 if inlist(birth_state_country,4,8,16,35,30,49,32,56)
qui replace born_here_state106_ = 1 if inlist(birth_state_country,4,8,16,35,30,49,32,56)
qui replace born_here_state107_ = 1 if inlist(birth_state_country,6,41,53)
qui replace born_here_state108_ = 1 if inlist(birth_state_country,6,41,53)
qui replace born_here_state109_ = 1 if inlist(birth_state_country,2)
qui replace born_here_state110_ = 1 if inlist(birth_state_country,15)

qui replace born_here_div1_ = 1 if inlist(birthDiv,5)
qui replace born_here_div2_ = 1 if inlist(birthDiv,7)
qui replace born_here_div3_ = 1 if inlist(birthDiv,5)
qui replace born_here_div4_ = 1 if inlist(birthDiv,1)
qui replace born_here_div5_ = 1 if inlist(birthDiv,3)
qui replace born_here_div6_ = 1 if inlist(birthDiv,3,6)
qui replace born_here_div7_ = 1 if inlist(birthDiv,3)
qui replace born_here_div8_ = 1 if inlist(birthDiv,3)
qui replace born_here_div9_ = 1 if inlist(birthDiv,7)
qui replace born_here_div10_ = 1 if inlist(birthDiv,8)
qui replace born_here_div11_ = 1 if inlist(birthDiv,3)
qui replace born_here_div12_ = 1 if inlist(birthDiv,7)
qui replace born_here_div13_ = 1 if inlist(birthDiv,3)
qui replace born_here_div14_ = 1 if inlist(birthDiv,4)
qui replace born_here_div15_ = 1 if inlist(birthDiv,6)
qui replace born_here_div16_ = 1 if inlist(birthDiv,9)
qui replace born_here_div17_ = 1 if inlist(birthDiv,5)
qui replace born_here_div18_ = 1 if inlist(birthDiv,3)
qui replace born_here_div19_ = 1 if inlist(birthDiv,3,4)
qui replace born_here_div20_ = 1 if inlist(birthDiv,1,2)
qui replace born_here_div21_ = 1 if inlist(birthDiv,2,5)
qui replace born_here_div22_ = 1 if inlist(birthDiv,8)
qui replace born_here_div23_ = 1 if inlist(birthDiv,2)
qui replace born_here_div24_ = 1 if inlist(birthDiv,9)
qui replace born_here_div25_ = 1 if inlist(birthDiv,1)
qui replace born_here_div26_ = 1 if inlist(birthDiv,5)
qui replace born_here_div27_ = 1 if inlist(birthDiv,9)
qui replace born_here_div28_ = 1 if inlist(birthDiv,9)
qui replace born_here_div29_ = 1 if inlist(birthDiv,9)
qui replace born_here_div30_ = 1 if inlist(birthDiv,9)
qui replace born_here_div31_ = 1 if inlist(birthDiv,9)
qui replace born_here_div32_ = 1 if inlist(birthDiv,3,4)
qui replace born_here_div33_ = 1 if inlist(birthDiv,5)
qui replace born_here_div34_ = 1 if inlist(birthDiv,5)
qui replace born_here_div35_ = 1 if inlist(birthDiv,5)
qui replace born_here_div36_ = 1 if inlist(birthDiv,1)
qui replace born_here_div37_ = 1 if inlist(birthDiv,1)
qui replace born_here_div38_ = 1 if inlist(birthDiv,2)
qui replace born_here_div39_ = 1 if inlist(birthDiv,2)
qui replace born_here_div40_ = 1 if inlist(birthDiv,3)
qui replace born_here_div41_ = 1 if inlist(birthDiv,3)
qui replace born_here_div42_ = 1 if inlist(birthDiv,4)
qui replace born_here_div43_ = 1 if inlist(birthDiv,4)
qui replace born_here_div44_ = 1 if inlist(birthDiv,5)
qui replace born_here_div45_ = 1 if inlist(birthDiv,5)
qui replace born_here_div46_ = 1 if inlist(birthDiv,6)
qui replace born_here_div47_ = 1 if inlist(birthDiv,6)
qui replace born_here_div48_ = 1 if inlist(birthDiv,7)
qui replace born_here_div49_ = 1 if inlist(birthDiv,7)
qui replace born_here_div50_ = 1 if inlist(birthDiv,8)
qui replace born_here_div51_ = 1 if inlist(birthDiv,8)
qui replace born_here_div52_ = 1 if inlist(birthDiv,9)
qui replace born_here_div53_ = 1 if inlist(birthDiv,9)
qui replace born_here_div54_ = 1 if inlist(birthDiv,9)
qui replace born_here_div55_ = 1 if inlist(birthDiv,9)
qui replace born_here_div56_ = 1 if inlist(birthDiv,5)
qui replace born_here_div57_ = 1 if inlist(birthDiv,7)
qui replace born_here_div58_ = 1 if inlist(birthDiv,5)
qui replace born_here_div59_ = 1 if inlist(birthDiv,1)
qui replace born_here_div60_ = 1 if inlist(birthDiv,3)
qui replace born_here_div61_ = 1 if inlist(birthDiv,3,6)
qui replace born_here_div62_ = 1 if inlist(birthDiv,3)
qui replace born_here_div63_ = 1 if inlist(birthDiv,3)
qui replace born_here_div64_ = 1 if inlist(birthDiv,7)
qui replace born_here_div65_ = 1 if inlist(birthDiv,8)
qui replace born_here_div66_ = 1 if inlist(birthDiv,3)
qui replace born_here_div67_ = 1 if inlist(birthDiv,7)
qui replace born_here_div68_ = 1 if inlist(birthDiv,3)
qui replace born_here_div69_ = 1 if inlist(birthDiv,4)
qui replace born_here_div70_ = 1 if inlist(birthDiv,6)
qui replace born_here_div71_ = 1 if inlist(birthDiv,9)
qui replace born_here_div72_ = 1 if inlist(birthDiv,5)
qui replace born_here_div73_ = 1 if inlist(birthDiv,3)
qui replace born_here_div74_ = 1 if inlist(birthDiv,3,4)
qui replace born_here_div75_ = 1 if inlist(birthDiv,1,2)
qui replace born_here_div76_ = 1 if inlist(birthDiv,2,5)
qui replace born_here_div77_ = 1 if inlist(birthDiv,8)
qui replace born_here_div78_ = 1 if inlist(birthDiv,2)
qui replace born_here_div79_ = 1 if inlist(birthDiv,9)
qui replace born_here_div80_ = 1 if inlist(birthDiv,1)
qui replace born_here_div81_ = 1 if inlist(birthDiv,5)
qui replace born_here_div82_ = 1 if inlist(birthDiv,9)
qui replace born_here_div83_ = 1 if inlist(birthDiv,9)
qui replace born_here_div84_ = 1 if inlist(birthDiv,9)
qui replace born_here_div85_ = 1 if inlist(birthDiv,9)
qui replace born_here_div86_ = 1 if inlist(birthDiv,9)
qui replace born_here_div87_ = 1 if inlist(birthDiv,3,4)
qui replace born_here_div88_ = 1 if inlist(birthDiv,5)
qui replace born_here_div89_ = 1 if inlist(birthDiv,5)
qui replace born_here_div90_ = 1 if inlist(birthDiv,5)
qui replace born_here_div91_ = 1 if inlist(birthDiv,1)
qui replace born_here_div92_ = 1 if inlist(birthDiv,1)
qui replace born_here_div93_ = 1 if inlist(birthDiv,2)
qui replace born_here_div94_ = 1 if inlist(birthDiv,2)
qui replace born_here_div95_ = 1 if inlist(birthDiv,3)
qui replace born_here_div96_ = 1 if inlist(birthDiv,3)
qui replace born_here_div97_ = 1 if inlist(birthDiv,4)
qui replace born_here_div98_ = 1 if inlist(birthDiv,4)
qui replace born_here_div99_ = 1 if inlist(birthDiv,5)
qui replace born_here_div100_ = 1 if inlist(birthDiv,5)
qui replace born_here_div101_ = 1 if inlist(birthDiv,6)
qui replace born_here_div102_ = 1 if inlist(birthDiv,6)
qui replace born_here_div103_ = 1 if inlist(birthDiv,7)
qui replace born_here_div104_ = 1 if inlist(birthDiv,7)
qui replace born_here_div105_ = 1 if inlist(birthDiv,8)
qui replace born_here_div106_ = 1 if inlist(birthDiv,8)
qui replace born_here_div107_ = 1 if inlist(birthDiv,9)
qui replace born_here_div108_ = 1 if inlist(birthDiv,9)
qui replace born_here_div109_ = 1 if inlist(birthDiv,9)
qui replace born_here_div110_ = 1 if inlist(birthDiv,9)

preserve
	tempfile annual
	sum switch_loc
	keep if inlist(timmonth,4,8,20,32,44,56)
	xtset IDtilde timmonth
	drop switch_state switch_MSA switch_loc switch_county pop_cat_lag
	bys IDtilde (timmonth): generat byte switch_state  =  state[_n-1] ~=state[_n] & _n>1
	bys IDtilde (timmonth): generat byte switch_county =  switch_state | (county[_n-1]~=county[_n] & state[_n-1]==state[_n] & _n>1)
	bys IDtilde (timmonth): generat byte switch_MSA    = (cbsaname[_n-1]~=cbsaname[_n] & _n>1)
	bys IDtilde (timmonth): generat byte switch_pop_cat=  pop_cat[_n-1]~=pop_cat[_n] & _n>1
	bys IDtilde (timmonth): generat byte switch_house  = mover==2
	bys IDtilde (timmonth): generat pop_cat_lag        =  pop_cat[_n-1]
	generat switch_division = 0
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"CT","ME","MA","NH","RI","VT")                & ~inlist(stabb[_n],"CT","ME","MA","NH","RI","VT")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"NJ","NY","PA")                               & ~inlist(stabb[_n],"NJ","NY","PA")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"IN","IL","MI","OH","WI")                     & ~inlist(stabb[_n],"IN","IL","MI","OH","WI")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"IA","KS","MN","MO","NE","ND","SD")           & ~inlist(stabb[_n],"IA","KS","MN","MO","NE","ND","SD")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"DE","DC","FL","GA","MD","NC","SC","VA","WV") & ~inlist(stabb[_n],"DE","DC","FL","GA","MD","NC","SC","VA","WV")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AL","KY","MS","TN")                          & ~inlist(stabb[_n],"AL","KY","MS","TN")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AR","LA","OK","TX")                          & ~inlist(stabb[_n],"AR","LA","OK","TX")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AZ","CO","ID","NM","MT","UT","NV","WY")      & ~inlist(stabb[_n],"AZ","CO","ID","NM","MT","UT","NV","WY")
	bys IDtilde (timmonth): replace switch_division = 1 if inlist(stabb[_n-1],"AK","CA","HI","OR","WA")                     & ~inlist(stabb[_n],"AK","CA","HI","OR","WA")

	bys IDtilde (timmonth): gen choice27_lag  = choice27[_n-1]
	bys IDtilde (timmonth): gen choice28_lag  = choice28[_n-1]
	bys IDtilde (timmonth): gen choice29_lag  = choice29[_n-1]
	bys IDtilde (timmonth): gen choice48_lag  = choice48[_n-1]
	bys IDtilde (timmonth): gen choice49_lag  = choice49[_n-1]
	bys IDtilde (timmonth): gen choice50_lag  = choice50[_n-1]
	bys IDtilde (timmonth): gen choice53_lag  = choice53[_n-1]
	bys IDtilde (timmonth): gen choice54b_lag = choice54b[_n-1]
	bys IDtilde (timmonth): gen choice55_lag  = choice55[_n-1]

	bys IDtilde (timmonth): gen choice54_lag  = choice54[_n-1]
	bys IDtilde (timmonth): gen choice56_lag  = choice56[_n-1]
	bys IDtilde (timmonth): gen choice58_lag  = choice58[_n-1]
	bys IDtilde (timmonth): gen choice96_lag  = choice96[_n-1]
	bys IDtilde (timmonth): gen choice98_lag  = choice98[_n-1]
	bys IDtilde (timmonth): gen choice100_lag = choice100[_n-1]
	bys IDtilde (timmonth): gen choice106_lag = choice106[_n-1]
	bys IDtilde (timmonth): gen choice108_lag = choice108[_n-1]
	bys IDtilde (timmonth): gen choice110_lag = choice110[_n-1]
	bys IDtilde (timmonth): gen empFT_lag     = empFT[_n-1]
	bys IDtilde (timmonth): gen emp_lag       = emp[_n-1]
	bys IDtilde (timmonth): gen inlf_lag      = inlf[_n-1]
	drop if wavemap==1
	drop born_here_state* born_here_div*
	xtset IDtilde year
	sum switch_loc
	compress
	save `annual', replace
	save "[REDACTED]sippCombinedNHWmaleAnnual.dta", replace
	! gzip -f [REDACTED]sippCombinedNHWmaleAnnual.dta
restore

preserve
	keep if inlist(timmonth,8,20,32,44,56)
	xtset IDtilde timmonth
	local bstates born_here_state1_ born_here_state2_ born_here_state3_ born_here_state4_ born_here_state5_ born_here_state6_ born_here_state7_ born_here_state8_ born_here_state9_ born_here_state10_ born_here_state11_ born_here_state12_ born_here_state13_ born_here_state14_ born_here_state15_ born_here_state16_ born_here_state17_  born_here_state18_  born_here_state19_  born_here_state20_  born_here_state21_  born_here_state22_  born_here_state23_  born_here_state24_  born_here_state25_  born_here_state26_  born_here_state27_  born_here_state28_  born_here_state29_  born_here_state30_  born_here_state31_  born_here_state32_  born_here_state33_  born_here_state34_  born_here_state35_  born_here_state36_  born_here_state37_  born_here_state38_  born_here_state39_  born_here_state40_  born_here_state41_  born_here_state42_  born_here_state43_  born_here_state44_  born_here_state45_  born_here_state46_  born_here_state47_  born_here_state48_  born_here_state49_  born_here_state50_  born_here_state51_  born_here_state52_ born_here_state53_  born_here_state54_  born_here_state55_  born_here_state56_  born_here_state57_  born_here_state58_  born_here_state59_  born_here_state60_  born_here_state61_  born_here_state62_  born_here_state63_  born_here_state64_  born_here_state65_  born_here_state66_  born_here_state67_  born_here_state68_  born_here_state69_  born_here_state70_  born_here_state71_  born_here_state72_  born_here_state73_  born_here_state74_  born_here_state75_  born_here_state76_  born_here_state77_  born_here_state78_  born_here_state79_  born_here_state80_  born_here_state81_  born_here_state82_  born_here_state83_  born_here_state84_  born_here_state85_  born_here_state86_  born_here_state87_  born_here_state88_  born_here_state89_  born_here_state90_  born_here_state91_  born_here_state92_  born_here_state93_  born_here_state94_  born_here_state95_  born_here_state96_  born_here_state97_  born_here_state98_  born_here_state99_  born_here_state100_ born_here_state101_ born_here_state102_ born_here_state103_ born_here_state104_ born_here_state105_ born_here_state106_ born_here_state107_ born_here_state108_ born_here_state109_ born_here_state110_
	local bdivs   born_here_div1_ born_here_div2_ born_here_div3_ born_here_div4_ born_here_div5_ born_here_div6_ born_here_div7_ born_here_div8_ born_here_div9_ born_here_div10_ born_here_div11_ born_here_div12_ born_here_div13_ born_here_div14_ born_here_div15_ born_here_div16_ born_here_div17_  born_here_div18_  born_here_div19_  born_here_div20_  born_here_div21_  born_here_div22_  born_here_div23_  born_here_div24_  born_here_div25_  born_here_div26_  born_here_div27_  born_here_div28_  born_here_div29_  born_here_div30_  born_here_div31_  born_here_div32_  born_here_div33_  born_here_div34_  born_here_div35_  born_here_div36_  born_here_div37_  born_here_div38_  born_here_div39_  born_here_div40_  born_here_div41_  born_here_div42_  born_here_div43_  born_here_div44_  born_here_div45_  born_here_div46_  born_here_div47_  born_here_div48_  born_here_div49_  born_here_div50_  born_here_div51_  born_here_div52_ born_here_div53_  born_here_div54_  born_here_div55_  born_here_div56_  born_here_div57_  born_here_div58_  born_here_div59_  born_here_div60_  born_here_div61_  born_here_div62_  born_here_div63_  born_here_div64_  born_here_div65_  born_here_div66_  born_here_div67_  born_here_div68_  born_here_div69_  born_here_div70_  born_here_div71_  born_here_div72_  born_here_div73_  born_here_div74_  born_here_div75_  born_here_div76_  born_here_div77_  born_here_div78_  born_here_div79_  born_here_div80_  born_here_div81_  born_here_div82_  born_here_div83_  born_here_div84_  born_here_div85_  born_here_div86_  born_here_div87_  born_here_div88_  born_here_div89_  born_here_div90_  born_here_div91_  born_here_div92_  born_here_div93_  born_here_div94_  born_here_div95_  born_here_div96_  born_here_div97_  born_here_div98_  born_here_div99_  born_here_div100_ born_here_div101_ born_here_div102_ born_here_div103_ born_here_div104_ born_here_div105_ born_here_div106_ born_here_div107_ born_here_div108_ born_here_div109_ born_here_div110_
	keep  IDtilde timmonth SPANEL `bstates' `bdivs'
	order IDtilde timmonth SPANEL `bstates' `bdivs'
	reshape wide `bstates' `bdivs', i(IDtilde SPANEL) j(timmonth)
	order IDtilde SPANEL
	d
	ds
	ds, alpha
	outsheet using "[REDACTED]sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv", comma replace nolabel
	! gzip -f [REDACTED]sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv
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
	gen calmo = calmonth
	gen calyr = year
	gen earnflag   = (empFT & inrange(lnearnfinal,6,10) & ~anyFlag)
	gen earnflagJ  = (empFT & inrange(lnearnfinal,6,10) & ~anyFlag)
	gen earnflagJb = (empFT & inrange(lnearnfinal,6,10) & ~anyFlag)
	* keep if ~anyFlag
	* drop if schoolFlag
	sum lnWage* lnearnfinal* constant pop_cat hgc educlevel age exper
	replace lnWageHr      = -999 if ~earnflag
	replace lnWageHrJ     = -999 if ~earnflagJ
	replace lnWageHrJb    = -999 if ~earnflagJb
	replace lnearnfinal   = -999 if ~earnflag
	replace lnearnfinalJ  = -999 if ~earnflagJ
	replace lnearnfinalJb = -999 if ~earnflagJb
	sum lnWage* lnearn* constant pop_cat hgc educlevel age exper
	keep IDtilde timmonth SPANEL constant weightlong `outcomes' `Evars' `Gvars' `Svars' calmo calyr
	reshape wide `outcomes' `Evars' `Gvars' `Svars' calmo calyr, i(IDtilde SPANEL) j(timmonth)
	order IDtilde SPANEL constant weightlong
	outsheet using "[REDACTED]sippCombinedNHWmale_matlab_wide_annual.csv", comma replace nolabel
	! gzip -f [REDACTED]sippCombinedNHWmale_matlab_wide_annual.csv
restore

log close

