options fullstimer;
options nocenter ls=100 obs=max;

%macro export_stata(year);

%include "../macros/macro.sipp_defs_jmp.sas";
%define(&year.);
%include "../macros/sippcorevars&year._jmp.sas";
%include "../macros/sipp_new_corevars&year._jmp.sas";
%include "../macros/labelme.sas";

/*--------------------------------------->
libname misc "~/work";
libname xportout xport "[REDACTED]S&year.NHWM.xpt";
libname xportout xport "[REDACTED]";


* x gunzip -f longcombinedpuf&year..sas7bdat.gz

data S&year.NHWM
	set jmpdata.longcombinedpuf&year.;
run;


proc copy in=misc out=xportout memtype=data;
	select S&year.NHWM;
run;



data xportout.S&year.NHWM.xpt;
	set S&year.NHWM;
run;


<---------------------------------------*/

data tempy;
	set jmpdata.longcombinedpuf&year.;
	/* %if &year eq 2004 %then %do */
	if puid   = "[REDACTED]" and &year. eq 2004 then puid = "[REDACTED]";
	/* %end
	* %if &year eq 2008 %then %do */
	if puid   = "[REDACTED]" and &year. eq 2008 then puid = "[REDACTED]";
	/* %end	
	*if lgtkey = [REDACTED] then lgtkey = "[REDACTED]";*/
run;

proc print data=tempy (firstobs=1 obs=50);
var puid lgtkey;
run;


proc freq data=tempy;
tables calmonth*srefmon year*srefmon;
run;

proc export data=tempy
	outfile = "[REDACTED]sipp&year.NHWmale.csv"
	dbms = csv
	replace;
run;

* x gzip -f [REDACTED]sipp&year.NHWmale.csv

%mend;
%export_stata(2004);
%export_stata(2008);
