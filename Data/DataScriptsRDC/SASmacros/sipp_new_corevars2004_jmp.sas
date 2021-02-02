/*
%let pufcntrl=
SSUID
SROTATON
SREFMON
EENTAID
EPPPNUM
LGTKEY
;

%let newpufpanel=
qhispanic
hispanic
born_us
qborn_us
;
*/

%let newpufmnth=
hhtype
hhsize
age
marst
enrlsch
empstatwk1
empstatwk2
empstatwk3
empstatwk4
empstatwk5
weekswjob
weeksabjob
rsn_stop_j1
rsn_stop_j2
rsn_no_work
rsn_absent
on_layoff
weight
monthearnj1
monthearnj2
multlocj1
multlocj2
firmsizej1
firmsizej2
stillatj1
stillatj2
sdatej1
sdatej2
edatej1
edatej2
classj1
classj2
wavemap
year
calmonth
rhrsweek
;

%let newiuwavegeo=
county
state
;

%let newpufwave =
sex
race
year
calmonth
intstatus
hgc_cat
ged
enrsch
enrwhich
hrswork
ind1
occ1
ind2
occ2
job1ID
job2ID
qjobflag
contingent
numjobs
mover
hourly1
qhourly1
hourly2
qhourly2
look_work
qlook_work
prop1
prop2
workedpt1plusweeks
qworkedpt1plusweeks
reasonworkpt
qreasonworkpt
;

%let newpufimpwave=
qrace
qhgc
qged
qenrwhich
qhrswork
qind1
qocc1
qind2
qocc2
;

%let newpufimpmnth=
qenrlsch
qrsn_stop_j1
qrsn_stop_j2
qrsn_no_work
qrsn_absent
qon_layoff
qmonthearnj1
qmonthearnj2
qmultlocj1
qmultlocj2
qfirmsizej1
qfirmsizej2
qstillatj1
qstillatj2
qsdatej1
qsdatej2
qedatej1
qedatej2
qclassj1
qclassj2
;

%let pufcntrllabels=
"Sample Unit Identifier"
"Wave of data collection"  
"Rotation within wave" 
"Reference month of this record"     
"Address ID of hhld where person entered sample"    
"Person number"    
"Person longitudinal key"
;

%let pufpanellabels=
"Imputation flag for Hispanic Origin"   
"Hispanic origin"     
"Respondent was born in the U.S."  
"Imputation flag for Respondent was born in the U.S."   
;

%let newpufmnthlabels=
"Household type"     
"Total number of persons in this household in this month" 
"Age as of last birthday (should use ad rec for true age)" 
"Marital Status"                
"Was R enrolled in school in this month?"   
"Emp status week 1" 
"Emp status week 2" 
"Emp status week 3" 
"Emp status week 4" 
"Emp status week 5" 
"Weeks with a job"  
"Weeks absent from job" 
"Reason stopped working for job 1"
"Reason stopped working for job 2"  
"Reason no job" 
"Reason absent from work"  
"On layoff?" 
"Person weight" 
"Monthly earnings at job 1" 
"Monthly earnings at job 2" 
"Employer operations in more than one location" 
"Employer operations in more than one location" 
"Number of employees at workers location" 
"Number of employees at workers location" 
"Still working for this employer?" 
"Still working for this employer?" 
"Job 1 Start Date" 
"Job 2 Start Date" 
"Job 1 End Date" 
"Job 2 End Date" 
"Class of worker at Job 1"  
"Class of worker at Job 2"  
;

%let newiuwavegeolabels=
"FIPS county code" 
"FIPS state code" 
;

%let newpufwavelabels =
"sex (1=male, 2=female)" 
"race (NHW, NHB, NHAsian, Other)" 
"Persons interview status" 
"Highest Degree Received or grade completed" 
"received GED" 
"Was R enrolled in school in all four months?" 
"At what level or grade was...enrolled?" 
"Usual hours worked per week at all jobs during the reference period" 
"Usual hours worked per week recode in month" 
"Industry at job 1 (5 digits)" 
"Occupation at job 1 (5 digits)" 
"Industry at job 2 (5 digits)" 
"Occupation at job 2 (5 digits)" 
"Unique ID of job 1" 
"Unique ID of job 2" 
"Flag indicating worker with unknown job dates" 
"Number of jobs held in reference period" 
"Mover flag" 
"Paid by the hour at job 1"
"Imp flag for paid by the hour at job 1"
"Paid by the hour at job 2"
"Imp flag for paid by the hour at job 2" 
"Looked for work during reference period"
"Imp flag for looked for work during reference period"
;

%let newpufimpwavelabels=
"Imputation flag for race" 
"Imputation flag for hgc"  
"Imputation flag for received GED" 
"Imputation flag for what level or grade was...enrolled?" 
"Imp flag for usual hours worked per week at all jobs during the reference period" 
"Imputation flag for industry" 
"Imputation flag for occupation" 
"Imputation flag for industry" 
"Imputation flag for occupation" 
;

%let newpufimpmnthlabels=
"Imp flag for Was R enrolled in school in this month?" 
"Imp flag for Reason stopped working for job 1" 
"Imp flag for Reason stopped working for job 2" 
"Imp flag for Reason no job" 
"Imp flag for Reason absent from work" 
"Imp flag for On layoff?" 
"Imp flag for Monthly earnings at job 1" 
"Imp flag for Monthly earnings at job 2" 
"Imp flag for Employer operations in more than one location" 
"Imp flag for Employer operations in more than one location" 
"Imp flag for Number of employees at workers location" 
"Imp flag for Number of employees at workers location" 
"Imp flag for Still working for this employer?" 
"Imp flag for Still working for this employer?" 
"Imp flag for Job 1 Start Date" 
"Imp flag for Job 2 Start Date" 
"Imp flag for Job 1 End Date" 
"Imp flag for Job 2 End Date" 
"Imp flag for Class of worker at Job 1" 
"Imp flag for Class of worker at Job 2" 
;
