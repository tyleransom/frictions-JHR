clear all
version 14.1
set more off
capture log close

log using GrapherCfls.log, replace

import excel using ../Paper/Fall2018/allCflsLong.xls, clear firstrow

*------------------------------------------------------------------------------
* Add labels to the data
*------------------------------------------------------------------------------
* convert outcome from string to categorical
generat outcome = .
replace outcome = 1 if strpos(outc,"Mig")
replace outcome = 2 if strpos(outc,"Unemp")
replace outcome = 3 if strpos(outc,"LFP")
lab def vloutcome 1 "Migration" 2 "Unemployment" 3 "LFP"
lab val outcome vloutcome
drop outc
order outcome 

* add value labels for counterfactuals
lab def vlcfl 1 "indep w decline" 2 "corr w decline" 3 "indep u increase" 4 "corr u increase" 5 "zero search costs" 6 "$10k moving subsidy" 7 "baseline move prob"
lab val Cfl vlcfl

* add value labels to employment
recode employed (0 = 2)
lab def vlemp 1 "Employed" 2 "Unemployed"
lab val employed vlemp

* add more info about locations
gen amenities = 2*inlist(loc,1)+0*inlist(loc,2)+1*inlist(loc,3,4,5,6)
gen earnings  = 2*inlist(loc,3)+0*inlist(loc,4)+1*inlist(loc,1,2,5,6)
gen emp_prob  = 2*inlist(loc,5)+0*inlist(loc,6)+1*inlist(loc,1,2,3,4)

egen id = group(year outcome type loc employed)
reshape wide delta_p, i(id) j(Cfl)

set scheme s2color
* graphs
forval x=1/6 {
    forval t=1/2 {
        forval y=1/3 {
            if "`y'"=="1" local name = "Mig"
            if "`y'"=="1" local namelong = "migration probability"
            if "`y'"=="1" local nameshort = "mig. prob."
            if "`y'"=="1" local yrange = " -.04(.02).06"
            if "`y'"=="1" local yscale = " -.03(.02).05"
            if "`y'"=="2" local name = "Unemp"
            if "`y'"=="2" local namelong = "unemployment rate"
            if "`y'"=="2" local nameshort = "unemp. rate"
            if "`y'"=="2" local yrange = " -.04(.02).06"
            if "`y'"=="2" local yscale = " -.03(.02).05"
            if "`y'"=="3" local name = "LFP"
            if "`y'"=="3" local namelong = "LFP rate"
            if "`y'"=="3" local nameshort = "LFP rate"
            if "`y'"=="3" local yrange = " -.04(.02).06"
            if "`y'"=="3" local yscale = " -.03(.02).05"
            qui sum delta_p7 if year==2007 & outcome==`y' & type==`t' & loc==`x' & employed==1
            local baseE = `r(mean)'
            qui sum delta_p7 if year==2007 & outcome==`y' & type==`t' & loc==`x' & employed==2
            local baseU = `r(mean)'
            graph bar delta_p1 delta_p2 delta_p3 delta_p4 delta_p6 if year==2007 & outcome==`y' & type==`t' & loc==`x', over(employed) ///
            graphregion(color(white)) ytitle("Change in `namelong'") yscale(range(`yrange')) ylabel(`yscale') legend(label(1 "Independent decrease in w") label(2 "Correlated decrease in w") label(3 "Independent increase in UR") label(4 "Correlated increase in UR") label(5 "Move subsidy") cols(2) symxsize(10) keygap(1) ) bar(1, lcolor(black) fcolor(gs16)) bar(2, lcolor(black) fcolor(gs12)) bar(3 , lcolor(black) fcolor(gs8)) bar(4, lcolor(black) fcolor(gs4)) bar(5, lcolor(black) fcolor(gs1)) text(-.04 8 "baseline `nameshort': `:di %4.3f `baseE''", place(e)) text(-.04 58 "baseline `nameshort': `:di %4.3f `baseU''", place(e))
            graph export ../Paper/Fall2018/Cfl`name'Loc`x'T`t'.eps, replace
        }
    }
}

* graphs
forval x=1/6 {
    forval t=1/2 {
        forval y=1/3 {
            if "`y'"=="1" local name = "Mig"
            if "`y'"=="1" local namelong = "migration probability"
            if "`y'"=="1" local nameshort = "mig. prob."
            if "`y'"=="1" local yrange = " -.04(.02).06"
            if "`y'"=="1" local yscale = " -.03(.02).05"
            if "`y'"=="2" local name = "Unemp"
            if "`y'"=="2" local namelong = "unemployment rate"
            if "`y'"=="2" local nameshort = "unemp. rate"
            if "`y'"=="2" local yrange = " -.04(.02).06"
            if "`y'"=="2" local yscale = " -.03(.02).05"
            if "`y'"=="3" local name = "LFP"
            if "`y'"=="3" local namelong = "LFP rate"
            if "`y'"=="3" local nameshort = "LFP rate"
            if "`y'"=="3" local yrange = " -.04(.02).06"
            if "`y'"=="3" local yscale = " -.03(.02).05"
            qui sum delta_p7 if year==2007 & outcome==`y' & type==`t' & loc==`x' & employed==1
            local baseE = `r(mean)'
            qui sum delta_p7 if year==2007 & outcome==`y' & type==`t' & loc==`x' & employed==2
            local baseU = `r(mean)'
            graph bar delta_p1 delta_p2 delta_p3 delta_p4 delta_p6 if year==2007 & outcome==`y' & type==`t' & loc==`x', over(employed) ///
            graphregion(color(oucream)) plotregion(color(oucream)) ytitle("Change in `namelong'") yscale(range(`yrange')) ylabel(`yscale') legend(label(1 "Independent decrease in w") label(2 "Correlated decrease in w") label(3 "Independent increase in UR") label(4 "Correlated increase in UR") label(5 "Move subsidy") cols(2) symxsize(10) keygap(1) region(fcolor(oucream) lcolor(oucream)) ) bar(1, lcolor(oucrimson) fcolor(oucrimson) fintensity(inten0)) bar(2, lcolor(oucrimson) fcolor(oucrimson) fintensity(inten30)) bar(3 , lcolor(oucrimson) fcolor(oucrimson) fintensity(inten50)) bar(4, lcolor(oucrimson) fcolor(oucrimson) fintensity(inten80)) bar(5, lcolor(oucrimson) fcolor(oucrimson) fintensity(inten100)) text(-.04 8 "baseline `nameshort': `:di %4.3f `baseE''", place(e)) text(-.04 58 "baseline `nameshort': `:di %4.3f `baseU''", place(e))
            graph export ../Paper/Fall2018/creamCfl`name'Loc`x'T`t'.eps, replace
        }
    }
}

log close

