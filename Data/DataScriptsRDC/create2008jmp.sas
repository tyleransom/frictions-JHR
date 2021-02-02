/*************************************************************************/
/* Program: create2004jmp.sas                                            */
/* Original Authors:  Martha Stinson and Lisa Dragoset                   */
/* Updated by:  Melissa Bjelland, Chen Zhao, Tyler Ransom                */
/* Date: <2013-03-04>                                                    */
/* Purpose:  This program makes a person-month-level dataset of select   */
/* variables in the 2004 panel of the SIPP                               */
/*************************************************************************/

options fullstimer;
options nocenter ls=100 obs=max;

%macro createjmp(year);

%include "../macros/macro.sipp_defs_jmp.sas";
%define(&year.);
%include "../macros/sippcorevars&year._jmp.sas";
%include "../macros/sipp_new_corevars&year._jmp.sas";
%include "../macros/labelme.sas";

%do wave=1 %to &maxwave;

proc sort data=	
	  %if &year. eq 2001 %then %do;
	  m&year.w&wave..ppmpuw_with_labels_&wave.
	  %end;
	  %else %if &year eq 2004 or &year eq 2008 %then %do;
	  m&year.w&wave..ppmpuw&wave. 
	  %end;
	  (keep=&pufcntrl &pufpanel &pufwave &iuwavegeo &pufimpwave 
	  &pufmnth &pufimpmnth)
	  out=temp&wave;
   by ssuid eentaid epppnum;
run;

data pufpanel&year.w&wave./*(drop=ssuid eentaid epppnum 
     srefmon swave srotaton)*/;
   set temp&wave./*(keep=&pufcntrl &pufpanel)*/;
   by ssuid eentaid epppnum;
   length puid $19;
   puid=ssuid||put(eentaid,z3.)||put(epppnum,z4.);
run;

proc sort data=pufpanel&year.w&wave.;
  by puid;
run;

%end;

data finalstacked&year.;
      set pufpanel&year.w1;
run;

%do wave = 2 %to &maxwave.;
	data finalstacked&year.;
		set finalstacked&year. pufpanel&year.w&wave.;
	run;
%end;


data finalstacked&year.;
      set finalstacked&year.;
	timmonth = swave*4-(4-srefmon+1)+1;
run;

/******** rename variables *******/
data finalstacked&year.;
     set finalstacked&year.;
		/* do wave vars */
		%let m = 1;
		%do %until (%scan(&pufwave,&m)= );
			/*%put &m;*/
			%let lastvar =%scan(&pufwave,-1);
			%let varname1=%scan(&pufwave,&m);
			%let varname2=%scan(&newpufwave,&m);
			rename &varname1 = &varname2;
			/*
			%if %scan(&pufwave,&m) eq &&lastvar %then %do;
				output;
			%end;
			*/
			%let m=%eval(&m+1);
		%end;
		
		/* now do wave geo vars */
		%let m = 1;
		%do %until (%scan(&iuwavegeo,&m)= );
			/*%put &m;*/
			%let lastvar =%scan(&iuwavegeo,-1);
			%let varname1=%scan(&iuwavegeo,&m);
			%let varname2=%scan(&newiuwavegeo,&m);
			rename &varname1 = &varname2;
			/*
			%if %scan(&iuwavegeo,&m) eq &&lastvar %then %do;
				output;
			%end;
			*/
			%let m=%eval(&m+1);
		%end;
		
		/*now do wave imp vars */
		%let m = 1;
		%do %until (%scan(&pufimpwave,&m)= );
			/*%put &m;*/
			%let lastvar =%scan(&pufimpwave,-1);
			%let varname1=%scan(&pufimpwave,&m);
			%let varname2=%scan(&newpufimpwave,&m);
			rename &varname1 = &varname2;
			/*
			%if %scan(&pufimpwave,&m) eq &&lastvar %then %do;
				output;
			%end;
			*/
			%let m=%eval(&m+1);
		%end;
		/* do month vars */
		%let m = 1;
		%do %until (%scan(&pufmnth,&m)= );
			/*%put &m;*/
			%let lastvar =%scan(&pufmnth,-1);
			%let varname1=%scan(&pufmnth,&m);
			%let varname2=%scan(&newpufmnth,&m);
			rename &varname1 = &varname2;
			/*
			%if %scan(&pufmnth,&m) eq &&lastvar %then %do;
				output;
			%end;
			*/
			%let m=%eval(&m+1);
		%end;
		
		/*now do month imp vars */
		%let m = 1;
		%do %until (%scan(&pufimpmnth,&m)= );
			/*%put &m;*/
			%let lastvar =%scan(&pufimpmnth,-1);
			%let varname1=%scan(&pufimpmnth,&m);
			%let varname2=%scan(&newpufimpmnth,&m);
			rename &varname1 = &varname2;
			/*
			%if %scan(&pufimpmnth,&m) eq &&lastvar %then %do;
				output;
			%end;
			*/
			%let m=%eval(&m+1);
		%end;
run;

data holder;
     set finalstacked&year.;
     if puid~='[REDACTED]' then delete;
run;

/*proc print data=holder;*/
/*var puid swave srefmon srotaton year age county state;*/
/*var puid swave srefmon srotaton rhcalyr tage gcounty gfipsst pop00;*/
/*run;*/

/******* merge wave 1 and wave 2 topical modules *******/
data tmw1vars&year.;
	set jmpdata.tmw1vars&year.;
run;

proc sort data = tmw1vars&year. tagsort;
	by puid;
run;

proc sort data = finalstacked&year. tagsort;
	by puid;
run;

data finalstacked&year.;
   merge finalstacked&year. (in=e) 
	 tmw1vars&year (in=f); 
	    by puid;
            if e and f then output;
run;

data tmw2vars&year.;
	set jmpdata.tmw2vars&year.;
run;

proc sort data = tmw2vars&year. tagsort;
	by puid;
run;

data finalstacked&year.;
   merge finalstacked&year. (in=g) 
	 tmw2vars&year (in=h); 
	    by puid;
	    if g and h then output;
run;


data holder;
     set finalstacked&year.;
     if puid~='[REDACTED]' then delete;
run;

proc print data=holder;
var puid wavemap srefmon srotaton year age county state yrenterlf yrexitlf mosoutlf exper0 birth_state_country yr_move_this_state yr_move_this_house;
run;

proc freq data=finalstacked&year.;
	tables year calmonth*year;
	title 'after concatenating all panels';
run;

/******** sort city data *******/
data citysize_long;
  set jmpdata.citysize_long;
length state year 6 county 4 default=4;
run;

proc sort data=citysize_long;
  by state county year;
run;

/******** sort main data *******/
data finalstacked&year.;
  set finalstacked&year.;
length state year 6 county 4 default=4;
run;

proc sort data=finalstacked&year.;
  by state county year;
run;

/******** merge city data *******/
data finalstacked&year.;
   merge finalstacked&year. (in=a) 
	 citysize_long (in=b); /* keep=county state pop pop_large year LandArea); */
	    by state county year;
            if a and b then output;
run;

/*
data finalstacked&year.;
	set finalstacked&year.;
	if puid = . then delete;
run;
*/

/******** merge longitudinal weights *
data lgtwgt&year.;
	set m&year.w12.lgtwgt&year.w12;
run;

proc sort data = lgtwgt&year. tagsort;
	by ssuid epppnum;
run;

proc sort data = finalstacked&year. tagsort;
	by ssuid epppnum;
run;

data finalstacked&year.;
   merge finalstacked&year. (in=k) 
	 lgtwgt&year. (in=j); 
	    by ssuid epppnum;
            if k and j then output;
run;
******/

proc freq data=finalstacked&year.;
	tables year;
	title 'after merging city data';
run;


/******** merge GSF data*******/
/* NOTES ABOUT ADMIN VARIABLES
Post 1978, wqc variables are SSI earnings equivalents
Pre-1978, the definition matches more closely with 
actual work history.
But you also need to be careful about FICA-covered
jobs
Best definition to use is sum(flag_pos_earn)
Martha suggests: 
for pre-1978 work histories are OK
for post-1978, use sum(flag_pos_earn)

For earnings: use sum of DER FICA and DER NON-FICA
*/
data gsf&year. (drop=panel);
	set sasdata.gsf_ssa_v6_0 (keep=pik puid panel flag_valid_SSN flag_earn1937_to_1951 flag_in_der flag_in_ser earn1937_to_1951 totearn_ser_: total_der_: pos_der_: wqc: brthmn_final brthyr_final);
	if panel~=2008 then delete;
run;

/*proc contents data=gsf&year.;
     title "CONTENTS of pared-down GSF";
run;*/
/* endsas; */
proc sort data = gsf&year. tagsort;
	by puid;
run;

proc sort data = finalstacked&year. tagsort;
	by puid;
run;

data finalstacked&year.;
   merge finalstacked&year. (in=r) 
	 gsf&year. (in=q);
	    by puid;
            if q and r then output;
run;

proc freq data=finalstacked&year.;
	tables year;
	title 'after merging GSF';
run;

/* reshape SER wage data */
data finalstacked&year.;
     set finalstacked&year.;
     
     totearn_ser = .;
     %do i=2008 %to 2013;
	  if year=&i. then totearn_ser = totearn_ser_&i.; 
     %end;
     
     quarter = .;
     if calmonth = 1  or calmonth=2  or calmonth=3  then quarter = 1;
     if calmonth = 4  or calmonth=5  or calmonth=6  then quarter = 2;
     if calmonth = 7  or calmonth=8  or calmonth=9  then quarter = 3;
     if calmonth = 10 or calmonth=11 or calmonth=12 then quarter = 4;
     
     /* initialize variables that will later be replaced */
     quarters_worked = .;
     working_this_qtr = .;
     /* now replace them */
     %do i=2008 %to 2013;
	  if year=&i. then quarters_worked = wqc_yrtot_&i.; 
          %do j=1 %to 4;
	      if year=&i. and quarter=&j. then working_this_qtr = wqc_&i.q&j.;
	  %end;
     %end;
run;

/******** print out test cases *******/
proc sort data=finalstacked&year.;
  by puid wavemap srefmon;
run;

data holder;
     set finalstacked&year.;
     if puid~='[REDACTED]' then delete;
run;
proc print data=holder;
var puid wavemap srefmon timmonth srotaton year age county state totearn_ser monthearnj1 monthearnj2 pop00 exper0 birth_state_country yr_move_this_state yr_move_this_house;
run;
proc print data=holder;
var puid wavemap srefmon timmonth age brthyr_final working_this_qtr;
run;

/* save full dataset */
data jmpdata.fulllongcombinedpuf&year.;
	set finalstacked&year.;
run;

/* keep only non-Hispanic white males (for now) */
data jmpdata.longcombinedpuf&year.;
	set finalstacked&year.;
	if race ~= 1 or sex ~= 1 or eorigin ~= 2 or epopstat~=1 then delete;
run;

/* create some variables for SAS analysis */
data finalstacked&year.2;
  set finalstacked&year.;
  if race ~= 1 or sex ~= 1 or eorigin ~= 2 or epopstat~=1 or wavemap~=1 or srefmon~=4 or age<18 or age>55 then delete;
  empFT = (rhrsweek=1);
  lpop = log(pop00);
  lden = log(dens00);
  learn = log(monthearnj1+monthearnj2);
  learn_ser = log(totearn_ser);
  lwage = log((monthearnj1+monthearnj2)/( weekswjob*hrswork));
  exper2 = exper0**2;
run;

/* summary statistics */
proc means data=finalstacked&year.2;
var wavemap srefmon exper0 empFT;
title 'Means for 18-55 NHW males in first wave reference month';
run;

proc logistic data = finalstacked&year.2;
model empFT (event='1') = exper0;
title 'Logit of full-time work on work experience';
run;

proc reg data = finalstacked&year.2;
model lwage = exper0 exper2 lpop;
title 'OLS of log wage on log population';
run;

proc reg data = finalstacked&year.2;
model lwage = exper0 exper2 lden;
title 'OLS of log wage on log density';
run;

proc reg data = finalstacked&year.2;
model learn = exper0 exper2 lpop;
title 'OLS of log monthly earnings on log population';
run;

proc reg data = finalstacked&year.2;
model learn_ser = exper0 exper2 lpop;
title 'OLS of log annual SER earnings on log population';
run;

proc contents data=jmpdata.longcombinedpuf&year.;
     title "CONTENTS of puf&year.";
run;

%mend;

%createjmp(2008);
