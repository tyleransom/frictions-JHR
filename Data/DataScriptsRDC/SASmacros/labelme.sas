%macro labelme(lbl_file = /* SAS dataset that x-walk var names to descriptions */,
               lbl_var  = /* Column listing the variables                      */,
	       lbl_desc = /* Column listing the labels (description)           */,
	       pds_lib  = /* Libref of the input file, WORK if file is in WORK */,
	       pds_file = /* SAS dataset needing labels (w/o libref prefix)    */
	       );

/*****************************************************************/
/**** [1] check that all parameters are entered into the macro ***/
/*****************************************************************/

%if %length(&lbl_file) = 0 %then %do; %goto ENDIT1; %end;
%if %length(&lbl_var ) = 0 %then %do; %goto ENDIT1; %end;
%if %length(&lbl_desc) = 0 %then %do; %goto ENDIT1; %end;
%if %length(&pds_lib ) = 0 %then %do; %goto ENDIT1; %end;
%if %length(&pds_file) = 0 %then %do; %goto ENDIT1; %end;

/**************************************************************************/
/**** [2] check if label file exists, then check if the variables exist ***/
/**************************************************************************/

%if %sysfunc(exist(&lbl_file)) eq 0 %then %do;
	%put ERROR: &lbl_file CROSSWALK FILE DOES NOT EXIST, CHECK SPELLLING;
	%goto ENDIT2;
%end;
%else %do;
	%let dsid = %sysfunc(open(&lbl_file));
	%let val1 = %sysfunc(varnum(&dsid, &lbl_var));
	%let val2 = %sysfunc(varnum(&dsid, &lbl_desc));
	%let rc   = %sysfunc(close(&dsid));
	
	%if &val1 eq 0 %then %do;
		%put ERROR: &lbl_var DOES NOT EXIST IN THE LABEL DATASET, CHECK SPELLING;
		%goto ENDIT2;
	%end;
	%if &val2 eq 0 %then %do;
		%put ERROR: &lbl_desc DOES NOT EXIST IN THE LABEL DATASET, CHECK SPELLING;
		%goto ENDIT2;
	%end;
%end;

/***********************************/
/**** [3] check for valid libref ***/
/***********************************/
%if %sysfunc(libref(&pds_lib)) eq 0 %then %do;
 /* do nothing - libref assigned */
%end;
%else %do;
	%put ERROR: LIBREF &pds_lib HAS NOT BEEN ASSIGNED, CHECK SPELLING;
	%goto ENDIT2;
%end;

/***************************************/
/**** [4] check if input file exists ***/
/***************************************/
%if %sysfunc(exist(&pds_lib..&pds_file)) eq 0 %then %do;
	%put ERROR: &pds_file DATASET DOES NOT EXIST IN LIBREF, CHECK SPELLING;
	%goto ENDIT2;
%end;

/***************************************************************************************/
/**** [5] if there is a tilde (~) in desc, change to (-), fix amperstand and percent ***/
/***************************************************************************************/
data temp0001;
	set &lbl_file (where=(&lbl_desc ne "") keep=&lbl_var &lbl_desc);
	temp1_    = compress(translate(&lbl_desc,"-","~"),"""");
	temp2_    = tranwrd(temp1_,"&","& ");
	temp3_    = tranwrd(temp2_,"%","% ");
	&lbl_desc = trim(left(compbl(temp3_)));
run;

/******************************************/
/**** [6] make two lists using PROC SQL ***/
/******************************************/
proc sql noprint;
	select &lbl_var into :nm separated by '~' from temp0001
	;
	select &lbl_desc into :lbl separated by '~' from temp0001
	;
quit;


/**************************************************/
/**** [7] write to LOG to visually verify lists ***/
/**************************************************/
%put nm list  = %bquote(&nm);
%put lbl list = %bquote(&lbl);

/***************************************************************/
/**** [8] finally, add labels to dataset using PROC DATASETS ***/
/***************************************************************/
proc datasets library=&pds_lib nolist;
	modify &pds_file;
	label
		%let i=1;
		%do %while(%length(%qscan(&nm,&i,%str(~))) > 0);
			%qscan(&nm,&i,%str(~)) = "%qscan(%bquote(&lbl),&i,%str(~))"
			%let i = %eval(&i+1);
		%end;
	;
run;
quit;

/*************************************/
/**** [9] delete data set temp0001 ***/
/*************************************/
proc datasets nolist;
	delete temp0001;
run;
quit;

%goto ENDIT2;

%ENDIT1:
	%PUT ERROR: MACRO STOPPED - MUST ENTER ALL PARAMETERS;

%ENDIT2:

%mend labelme;

