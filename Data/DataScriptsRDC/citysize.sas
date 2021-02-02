/*************************************************************************/
/* Program: citysize.sas                                                 */
/* Original Author:  Tyler Ransom                                        */
/* Date: <2012-08-28>                                                    */
/* Purpose:  This program imports city population data for merging with  */
/* the SIPP geographical variables (2001, 2004, 2008 panels only)        */
/*************************************************************************/

%include "../macros/macro.sipp_defs_jmp.sas";

options fullstimer;
options nocenter ls=100 obs=max;

/*
filename myfile "[REDACTED]FIPS_MSA_populations2000_2011.csv";
data varlabels;
	infile myfile dlm=',' dsd missover firstobs=2;
	length MSAname1 $100;
	length MSAname $100;
	length CBSAandDivisionTitlesandComponen $100;
	input county $ state $ SSA FIPS MSAname1 $ CBSAcode MetroDivisionCode CBSAandDivisionTitlesandComponen $ Population LandArea MSAname $
	pop2009 pop2008 pop2007 pop2006 pop2005 pop2004 pop2003 pop2002 pop2001 pop2000 pop2010 pop2011 pop_large2000 pop_large2001 pop_large2002
	pop_large2003 pop_large2004 pop_large2005 pop_large2006 pop_large2007 pop_large2008 pop_large2009 pop_large2010 pop_large2011 FIPSstring $
	FIPScounty $ FIPSstate $;
run;
*/
/* x gunzip [REDACTED]county_city_chars.csv.gz */
/* proc import datafile="[REDACTED]/county_city_chars.csv" */
proc import datafile="[REDACTED]county_city_chars.csv"
	out = citysize
	dbms=csv
	replace;
	getnames=yes;
run;
data citysize;
	set citysize;
	FIPScounty_new = input(county, best3.);
	FIPSstate_new  = input(state , best2.);
	FIPScounty     = FIPScounty_new;
	FIPSstate      = FIPSstate_new;
	/*gdpnew         = input(gdp, best32.);*/
run;

/*data citysize (drop = FIPScounty FIPSstate);
	set citysize;
run;
data citysize (rename = (FIPScounty_new = FIPScounty));
	set citysize;
run;
data citysize (rename = (FIPSstate_new = FIPSstate));
	set citysize;
run;
proc contents data=citysize;
run;
proc means data=citysize;
run;*/

data jmpdata.citysize_long;
	set citysize;
run;

/*
data jmpdata.citysize_long;
	set citysize;
	array apop(2000:2011) pop2000-pop2011;
	array apopbig(2000:2011) pop_large2000-pop_large2011;
	
	do year = 2000 to 2011;
		pop = apop(year);
		pop_large = apopbig(year);
		output;
	end;
	
	drop pop2000-pop2011 pop_large2000-pop_large2011;
run;
*/

proc contents data=jmpdata.citysize_long;
run;
proc means data=jmpdata.citysize_long;
run;
proc print data=jmpdata.citysize_long (firstobs= 1 obs = 500);
run;
* x gzip [REDACTED]county_city_chars.csv
