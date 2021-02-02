/*************************************************************************/
/* Program:  macro.sipp_defs_jmp.sas                                     */
/* Author:  Tyler Ransom (taken from Melissa Bjelland)                   */
/* Date: <2012-08-24 13:35>                                              */
/* Purpose:  This program defines key macro variables that are pulled    */
/* into the SIPP code to facilitate the renaming of directory structures */
/* when files are ported to a new location for revision.                 */
/*                                                                       */
/* Tyler made the following changes to this file:                        */
/* added a macro variable "impdata" to line 47                           */
/* added this macro variable to the list of globals on line 31           */
/* added a libname statement for "impdata" to line 74                    */
/* added a macro variable "jmpdata" to line 48                           */
/* added this macro variable to the list of globals on line 31           */
/* added a libname statement for "impdata" to line 75                    */
/*************************************************************************/
/* To use these definitions add the following   */
/* lines to your SAS program                    */
/*  %include '../../macros/macro.sipp_defs.sas' */
/* and call the macro defined here as           */
/*  %define(1990) (f.i.)                        */
/* Editing this file:                           */
/*   Feel free to add macro variables, but      */
/*   to work, they need to be referenced by     */
/*   the first command within the macro         */
/*   (the global statement)                     */

%put ::+++ It defines the following macro variables: ;


%let globvars=stem dbase pbase version impdata sasdata jmpdata
        pmacros pcf ser der f831 mbr ssr phus xwalk hcef
	iuxwalk ssaendyr aime cpibase intdob
        panel1 panel2 panel3 panel4 panel5 panel6 panel7
        first_mo last_mo first_yr last_yr sipp1mb sipp1me
        sipp2mb sipp2me sipp3mb sipp3me sipp4mb sipp4me 
	maxwave maxnum m_first_yr m_last_yr;

%put ::+++ &globvars.;
%put ::+++++++++++++++++++++++++++++++++++++++++++++++++++++;

%global &globvars;
%let stem=[REDACTED];
%let dbase=[REDACTED];
%let pbase=[REDACTED];
%let version=[REDACTED];
%let impdata=&stem.&dbase.&version./[REDACTED];
%let jmpdata=&stem.&dbase.&version./[REDACTED];
%let sasdata=&stem.&dbase.&version./[REDACTED];
%let pmacros=&stem.&pbase.&version./[REDACTED];
%let ssaendyr=[REDACTED];; /* SSA file last year - earnings files */
%let mbaendyr=[REDACTED];; /* SSA file MBR last year*/
%let cpibase=[REDACTED];; /* Base year for CPI deflation */

/*------- SSA File Names -------*/
%let pcf=[REDACTED];
%let ser=[REDACTED];
%let der=[REDACTED];
%let phus=[REDACTED];
%let f831=[REDACTED];
%let mbr=[REDACTED];
%let ssr=[REDACTED];
%let xwalk=[REDACTED];
%let iuxwalk=[REDACTED];
%let aime=[REDACTED];
**%let hcef=[REDACTED];
%let intdob=[REDACTED];


/*------- Library Names -------*/

libname sasdata "&sasdata.";
libname impdata "&impdata.";
libname jmpdata "&jmpdata.";
libname pcf [REDACTED];
libname ser [REDACTED];
libname der [REDACTED];
libname phus [REDACTED];
libname f831 [REDACTED];
libname mbr [REDACTED];
libname ssr [REDACTED];
libname aime [REDACTED];
libname xwalk [REDACTED];
libname iuxwalk [REDACTED];
libname isu1984 [REDACTED];
libname isu1990 [REDACTED];
libname isu1991 [REDACTED];
libname isu1992 [REDACTED];
libname isu1993 [REDACTED];
libname isu1996 [REDACTED];
libname i2u1996 [REDACTED];
libname i2u2001 [REDACTED];
libname i2u2004 [REDACTED];
libname i2u2008 [REDACTED];
libname intbday [REDACTED];
libname mkjobs [REDACTED];
**libname hcef "&stem.[REDACTED]";

%macro define(year);

%let panel1=1990;
%let panel2=1991;
%let panel3=1992;
%let panel4=1993;
%let panel5=1996;
%let panel6=2001;
%let panel7=2004;
%let panel8=2008;
%let panel9=1984;


*%let m_first_yr=1983; /* master first year */
*%let m_last_yr=2013;  /* master last year */


%if &year = 1990 %then %goto def1990;
%else %if &year = 1991 %then %goto def1991;
%else %if &year = 1992 %then %goto def1992;
%else %if &year = 1993 %then %goto def1993;
%else %if &year = 1996 %then %goto def1996;
%else %if &year = 2001 %then %goto def2001;
%else %if &year = 2004 %then %goto def2004;
%else %if &year = 2008 %then %goto def2008;
%else %if &year = 1984 %then %goto def1984;

/*-------------------- definitions: 1984 -------------------- */
/*Oddities of the 1984 panel:  2 rotation groups miss a wave*/
/*Rotation group 4:  misses wave 2 but questions in wave 3 are asked about wave 2 reference period*/
/*This has the effect of making rotation group 4 seem to be missing the last wave*/
/*Rotation group 3:  misses wave 8 but no shifting of following months is done*/
/*Data for this wave are missing for this entire rotation group*/
%def1984: ;
libname mast1984 "[REDACTED]";
%let maxwave = 9;
%let maxnum = %eval(&maxwave.*4);
%let first_mo=6;         /* first month in SIPP */
%let last_mo=8;           /* last month in SIPP */
%let first_yr=1983;       /* first year in SIPP */
%let last_yr=1986;        /* last year in SIPP */
%let sipp1mb=6;           /* rotation 1 month begin */
%let sipp1me=5;           /* rotation 1 month end */
%let sipp2mb=7;           /* rotation 2 month begin */
%let sipp2me=6;           /* rotation 2 month end */
%let sipp3mb=8;           /* rotation 3 month begin */
%let sipp3me=7;           /* rotation 3 month end */
%let sipp4mb=9;           /* rotation 4 month begin */
%let sipp4me=8;           /* rotation 4 month end */
/***************************************************
***************************************************/

%goto enddef;

/*-------------------- definitions: 1990 -------------------- */
%def1990: ; 
libname mast1990 "[REDACTED]"; 
%let maxwave = 8;
%let maxnum = %eval(&maxwave.*4);
%let first_mo=10;         /* first month in SIPP */
%let last_mo=8;           /* last month in SIPP */
%let first_yr=1989;       /* first year in SIPP */
%let last_yr=1992;        /* last year in SIPP */
%let sipp1mb=1;           /* rotation 1 month begin */
%let sipp1me=&last_mo.;   /* rotation 1 month end */
%let sipp2mb=&first_mo.;  /* rotation 2 month begin */
%let sipp2me=5;           /* rotation 2 month end */
%let sipp3mb=11;          /* rotation 3 month begin */
%let sipp3me=6;           /* rotation 3 month end */
%let sipp4mb=12;          /* rotation 4 month begin */
%let sipp4me=7;           /* rotation 4 month end */
/***************************************************
Note that rotation 1 month end is last_mo and
rotation 2 month begin is first_mo.  The rotation
groups are actually temporally ordered as:
 2, 3, 4, 1. 
***************************************************/

%goto enddef;

/*-------------------- definitions: 1991 -------------------- */
%def1991: ;
libname mast1991 "[REDACTED]"; 
%let maxwave =8;
%let maxnum =%eval(&maxwave.*4);
%let first_mo=10;         /* first month in SIPP */
%let last_mo=8;           /* last month in SIPP */
%let first_yr=1990;       /* first year in SIPP */
%let last_yr=1993;        /* last year in SIPP */
%let sipp1mb=1;           /* rotation 1 month begin */
%let sipp1me=&last_mo.;   /* rotation 1 month end */
%let sipp2mb=&first_mo.;  /* rotation 2 month begin */
%let sipp2me=5;           /* rotation 2 month end */
%let sipp3mb=11;          /* rotation 3 month begin */
%let sipp3me=6;           /* rotation 3 month end */
%let sipp4mb=12;          /* rotation 4 month begin */
%let sipp4me=7;           /* rotation 4 month end */
/***************************************************
Note that rotation 1 month end is last_mo and
rotation 2 month begin is first_mo.  The rotation
groups are actually temporally ordered as:
 2, 3, 4, 1. 
***************************************************/

%goto enddef;

/*-------------------- definitions: 1992 -------------------- */
%def1992: ;
libname mast1992 "[REDACTED]"; 
%let maxwave =9; /* really 10 */
%let maxnum =%eval(&maxwave.*4);
%let first_mo=10;         /* first month in SIPP */
%let last_mo=12;           /* last month in SIPP */
%let first_yr=1991;       /* first year in SIPP */
%let last_yr=1994;        /* last year in SIPP */
%let sipp1mb=1;           /* rotation 1 month begin */
%let sipp1me=&last_mo.;   /* rotation 1 month end */
%let sipp2mb=&first_mo.;  /* rotation 2 month begin */
%let sipp2me=9;           /* rotation 2 month end */
%let sipp3mb=11;          /* rotation 3 month begin */
%let sipp3me=10;           /* rotation 3 month end */
%let sipp4mb=12;          /* rotation 4 month begin */
%let sipp4me=11;           /* rotation 4 month end */
/***************************************************
Note that rotation 1 month end is last_mo and
rotation 2 month begin is first_mo.  The rotation
groups are actually temporally ordered as:
 2, 3, 4, 1. 

Also note that while the 1992 panel was slated to 
have 10 waves, only 9 wave files exist in the data
repository.
***************************************************/

%goto enddef;

/*-------------------- definitions: 1993 -------------------- */
%def1993: ;
libname mast1993 "[REDACTED]"; 
%let maxwave =9;
%let maxnum =%eval(&maxwave.*4);
%let first_mo=10;         /* first month in SIPP */
%let last_mo=12;           /* last month in SIPP */
%let first_yr=1992;       /* first year in SIPP */
%let last_yr=1995;        /* last year in SIPP */
%let sipp1mb=1;           /* rotation 1 month begin */
%let sipp1me=&last_mo.;   /* rotation 1 month end */
%let sipp2mb=&first_mo.;  /* rotation 2 month begin */
%let sipp2me=9;           /* rotation 2 month end */
%let sipp3mb=11;          /* rotation 3 month begin */
%let sipp3me=10;           /* rotation 3 month end */
%let sipp4mb=12;          /* rotation 4 month begin */
%let sipp4me=11;           /* rotation 4 month end */
/***************************************************
Note that rotation 1 month end is last_mo and
rotation 2 month begin is first_mo.  The rotation
groups are actually temporally ordered as:
 2, 3, 4, 1. 
***************************************************/

%goto enddef;


/*-------------------- definitions: 1996 -------------------- */
%def1996: ;
libname mast1996 "[REDACTED]";
libname m1996w1 "[REDACTED]"; 
libname m1996w2 "[REDACTED]"; 
libname m1996w3 "[REDACTED]"; 
libname m1996w4 "[REDACTED]"; 
libname m1996w5 "[REDACTED]"; 
libname m1996w6 "[REDACTED]"; 
libname m1996w7 "[REDACTED]"; 
libname m1996w8 "[REDACTED]"; 
libname m1996w9 "[REDACTED]"; 
libname m1996w10 "[REDACTED]"; 
libname m1996w11 "[REDACTED]"; 
libname m1996w12 "[REDACTED]"; 
%let maxwave =12;
%let maxnum =%eval(&maxwave.*4);
%let first_mo=12;         /* first month in SIPP */
%let last_mo=2;           /* last month in SIPP */
%let first_yr=1995;       /* first year in SIPP */
%let last_yr=2000;        /* last year in SIPP */
%let sipp4mb=3;           /* rotation 4 month begin */
%let sipp4me=&last_mo.;   /* rotation 4 month end */
%let sipp1mb=&first_mo.;  /* rotation 1 month begin */
%let sipp1me=11;           /* rotation 1 month end */
%let sipp2mb=1;          /* rotation 2 month begin */
%let sipp2me=12;           /* rotation 2 month end */
%let sipp3mb=2;          /* rotation 3 month begin */
%let sipp3me=1;           /* rotation 3 month end */
/***************************************************
The rotation groups are temporally ordered as:
 1, 2, 3, 4. 
***************************************************/

%goto enddef;


/*-------------------- definitions: 2001 -------------------- */
%def2001: ;
libname m2001w1 "[REDACTED]"; 
libname m2001w2 "[REDACTED]"; 
libname m2001w3 "[REDACTED]"; 
libname m2001w4 "[REDACTED]"; 
libname m2001w5 "[REDACTED]"; 
libname m2001w6 "[REDACTED]"; 
libname m2001w7 "[REDACTED]"; 
libname m2001w8 "[REDACTED]"; 
libname m2001w9 "[REDACTED]"; 
%let maxwave =9;
%let maxnum =%eval(&maxwave.*4);
%let first_mo=10;         /* first month in SIPP */
%let last_mo=12;           /* last month in SIPP */
%let first_yr=2000;       /* first year in SIPP */
%let last_yr=2003;        /* last year in SIPP */
%let sipp4mb=1;           /* rotation 4 month begin */
%let sipp4me=&last_mo.;   /* rotation 4 month end */
%let sipp1mb=&first_mo.;  /* rotation 1 month begin */
%let sipp1me=9;           /* rotation 1 month end */
%let sipp2mb=11;          /* rotation 2 month begin */
%let sipp2me=10;          /* rotation 2 month end */
%let sipp3mb=12;          /* rotation 3 month begin */
%let sipp3me=11;          /* rotation 3 month end */
/***************************************************
The rotation groups are temporally ordered as:
 1, 2, 3, 4. 
***************************************************/

%goto enddef;



/*-------------------- definitions: 2004 -------------------- */
%def2004: ;
libname m2004w1  "[REDACTED]";
libname m2004w2  "[REDACTED]"; 
libname m2004w3  "[REDACTED]"; 
libname m2004w4  "[REDACTED]"; 
libname m2004w5  "[REDACTED]"; 
libname m2004w6  "[REDACTED]"; 
libname m2004w7  "[REDACTED]"; 
libname m2004w8  "[REDACTED]"; 
libname m2004w9  "[REDACTED]"; 
libname m2004w10 "[REDACTED]"; 
libname m2004w11 "[REDACTED]"; 
libname m2004w12 "[REDACTED]"; 
libname w2004w12 "[REDACTED]"; 
%let maxwave =12; 
%let maxnum =%eval(&maxwave.*4);
%let first_mo=10;         /* first month in SIPP */
%let last_mo=12;           /* last month in SIPP */
%let first_yr=2003;       /* first year in SIPP */
%let last_yr=2007;        /* last year in SIPP */
%let sipp4mb=1;           /* rotation 4 month begin */
%let sipp4me=&last_mo.;   /* rotation 4 month end */
%let sipp1mb=&first_mo.;  /* rotation 1 month begin */
%let sipp1me=9;           /* rotation 1 month end */
%let sipp2mb=11;          /* rotation 2 month begin */
%let sipp2me=10;           /* rotation 2 month end */
%let sipp3mb=12;          /* rotation 3 month begin */
%let sipp3me=11;           /* rotation 3 month end */
/***************************************************
The rotation groups are temporally ordered as:
 1, 2, 3, 4. 

For this version (Nov 2009) have included all 12 waves
***************************************************/

%goto enddef;

/*-------------------- definitions: 2008 -------------------- */
%def2008: ;
libname m2008w1  "[REDACTED]"; 
libname m2008w2  "[REDACTED]"; 
libname m2008w3  "[REDACTED]"; 
libname m2008w4  "[REDACTED]"; 
libname m2008w5  "[REDACTED]"; 
libname m2008w6  "[REDACTED]"; 
libname m2008w7  "[REDACTED]"; 
libname m2008w8  "[REDACTED]"; 
libname m2008w9  "[REDACTED]"; 
libname m2008w10 "[REDACTED]"; 
libname m2008w11 "[REDACTED]"; 
libname m2008w12 "[REDACTED]"; 
libname m2008w13 "[REDACTED]"; 
libname m2008w14 "[REDACTED]"; 
 
%let maxwave =14; 
%let maxnum =%eval(&maxwave.*4);
%let first_mo=5;         /* first month in SIPP */
%let last_mo=3;           /* last month in SIPP */
%let first_yr=2008;       /* first year in SIPP */
%let last_yr=2013;        /* last year in SIPP */
%let sipp4mb=8;           /* rotation 4 month begin */
%let sipp4me=&last_mo.;   /* rotation 4 month end */
%let sipp1mb=&first_mo.;  /* rotation 1 month begin */
%let sipp1me=12;           /* rotation 1 month end */
%let sipp2mb=6;          /* rotation 2 month begin */
%let sipp2me=1;           /* rotation 2 month end */
%let sipp3mb=7;          /* rotation 3 month begin */
%let sipp3me=2;           /* rotation 3 month end */
/***************************************************
The rotation groups are temporally ordered as:
 1, 2, 3, 4. 

For this version (Mar 2014) have included 14 waves
***************************************************/

%goto enddef;

%enddef: ;

%mend;
