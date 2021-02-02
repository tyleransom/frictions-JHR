%******************************************************************************
% Matlab code to read in data from Stata for structural estimation
%******************************************************************************
clear all; clc;
delete dataImportCombinedAnnual.diary;
diary  dataImportCombinedAnnual.diary;

tic

! test -f [REDACTED]sippCombinedNHWmale_matlab_wide_annual.csv.gz && gunzip -f [REDACTED]sippCombinedNHWmale_matlab_wide_annual.csv.gz || echo "sippCombinedNHWmale_matlab_wide_annual.csv.gz Not Found"
! head -1 [REDACTED]sippCombinedNHWmale_matlab_wide_annual.csv
A = importdata('[REDACTED]sippCombinedNHWmale_matlab_wide_annual.csv');
! gzip -f [REDACTED]sippCombinedNHWmale_matlab_wide_annual.csv

data = A.data;
clear A
[N KK] = size(data)
T  =  5; % maximum number of time periods in panel
K0 =  4; % number of time-invariant variables
K1 = 62; % number of time-varying variables

%==============================================================================
% Read in the data. Columns 1:K0 are time-invariant variables; the rest are
% time-varying regressors.
% Create a master NxKxT tensor that stores all the K variables for all N people
% in each of the T periods.
%==============================================================================
cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),N,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

ID               = all_vars(:,1 ,:);
panel            = all_vars(:,2 ,:); 
constant         = all_vars(:,3 ,:); 
weightlong       = all_vars(:,4 ,:);
age              = all_vars(:,5 ,:);
rhrsweek         = all_vars(:,6 ,:);
urate            = all_vars(:,7 ,:);
pop_cat          = all_vars(:,8 ,:);
lnWageHr         = all_vars(:,9 ,:); lnWageHr  (lnWageHr  ==-999) = NaN;
lnWageHrJ        = all_vars(:,10,:); lnWageHrJ (lnWageHrJ ==-999) = NaN;
lnWageHrJb       = all_vars(:,11,:); lnWageHrJb(lnWageHrJb==-999) = NaN;
hgc              = all_vars(:,12,:);
educlevel        = all_vars(:,13,:);
inlf             = all_vars(:,14,:);
empFT            = all_vars(:,15,:);
emp              = all_vars(:,16,:);
lnearnfinal      = all_vars(:,17,:); lnearnfinal  (lnearnfinal  ==-999) = NaN;
lnearnfinalJ     = all_vars(:,18,:); lnearnfinalJ (lnearnfinalJ ==-999) = NaN;
lnearnfinalJb    = all_vars(:,19,:); lnearnfinalJb(lnearnfinalJb==-999) = NaN;
exper            = all_vars(:,20,:);
experAlt         = all_vars(:,21,:);
experFTorPT      = all_vars(:,22,:);
choice29         = all_vars(:,23,:);
choice28         = all_vars(:,24,:);
choice27         = all_vars(:,25,:);
choice48         = all_vars(:,26,:);
choice49         = all_vars(:,27,:);
choice50         = all_vars(:,28,:);
choice53         = all_vars(:,29,:);
choice54b        = all_vars(:,30,:);
choice55         = all_vars(:,31,:);
choice54         = all_vars(:,32,:);
choice56         = all_vars(:,33,:);
choice58         = all_vars(:,34,:);
choice96         = all_vars(:,35,:);
choice98         = all_vars(:,36,:);
choice100        = all_vars(:,37,:);
choice106        = all_vars(:,38,:);
choice108        = all_vars(:,39,:);
choice110        = all_vars(:,40,:);
pop_cat_lag      = all_vars(:,41,:);
choice27_lag     = all_vars(:,42,:);
choice28_lag     = all_vars(:,43,:);
choice29_lag     = all_vars(:,44,:);
choice48_lag     = all_vars(:,45,:);
choice49_lag     = all_vars(:,46,:);
choice50_lag     = all_vars(:,47,:);
choice53_lag     = all_vars(:,48,:);
choice54b_lag    = all_vars(:,49,:);
choice55_lag     = all_vars(:,50,:);
choice54_lag     = all_vars(:,51,:);
choice56_lag     = all_vars(:,52,:);
choice58_lag     = all_vars(:,53,:);
choice96_lag     = all_vars(:,54,:);
choice98_lag     = all_vars(:,55,:);
choice100_lag    = all_vars(:,56,:);
choice106_lag    = all_vars(:,57,:);
choice108_lag    = all_vars(:,58,:);
choice110_lag    = all_vars(:,59,:);
empFT_lag        = all_vars(:,60,:);
emp_lag          = all_vars(:,61,:);
inlf_lag         = all_vars(:,62,:);
anyFlag          = all_vars(:,63,:);
calmo            = all_vars(:,64,:);
calyr            = all_vars(:,65,:);
earnflag         = all_vars(:,66,:);
collgrad         = educlevel==3;


%==============================================================================
% Summary Stats to check with Stata
%==============================================================================
all_long = reshape(permute(all_vars, [3 1 2]),N*T,(K0+K1));
stat_mat = all_long(:,[6:8 2 24 23 9:10 20:22 4 14 16]); 
% stata equivalent: sum lnWage* constant pop_cat logpop hgc educlevel lngdp m_lngdp lnipc age exper tenure
sum_stats = cat(2,nanmean(stat_mat)',nanstd(stat_mat)',nanmin(stat_mat)',nanmax(stat_mat)')
sum_N     = sum(~isnan(stat_mat(:,[4 1:3])))
summarize(stat_mat);

%==============================================================================
% Check basic estimation with Stata
%==============================================================================
wageflag = anyFlag==0 & earnflag==1; %& ~isnan(anyFlag) & ~isnan(empFT) & ~isnan(lnWageHr) & ~isnan(hgc) & ~isnan(age) & ~isnan(exper) & ~isnan(educlevel);
summarize(lnearnfinal(wageflag & educlevel==3));
summarize(lnearnfinal(wageflag & educlevel~=3));
% regress(lnWageHr(wageflag==1),[hgc(wageflag==1) age(wageflag==1) exper(wageflag==1) collgrad(wageflag==1) ones(sum(sum(wageflag==1)),1)])

choice6temp = nan(size(pop_cat));
choice6temp(pop_cat==1 & empFT==1) = 1;
choice6temp(pop_cat==2 & empFT==1) = 2;
choice6temp(pop_cat==3 & empFT==1) = 3;
choice6temp(pop_cat==1 & empFT==0) = 4;
choice6temp(pop_cat==2 & empFT==0) = 5;
choice6temp(pop_cat==3 & empFT==0) = 6;

choice6     = nan(size(pop_cat));
choice6_lag = nan(size(pop_cat_lag));
for j=1:3
    choice6    (pop_cat    ==j & empFT    ==1) = j;
    choice6    (pop_cat    ==j & empFT    ==0) = 3+j;
    choice6_lag(pop_cat_lag==j & empFT_lag==1) = j; 
    choice6_lag(pop_cat_lag==j & empFT_lag==0) = 3+j;
end

choice54     = nan(size(choice27));
choice54_lag = nan(size(choice27_lag));
for j=1:27
    choice54    (choice27    ==j & empFT    ==1) = j;
    choice54    (choice27    ==j & empFT    ==0) = 27+j;
    choice54_lag(choice27_lag==j & empFT_lag==1) = j; 
    choice54_lag(choice27_lag==j & empFT_lag==0) = 27+j;
end


choice56     = nan(size(choice28));
choice56_lag = nan(size(choice28_lag));
for j=1:28
    choice56    (choice28    ==j & empFT    ==1) = j;
    choice56    (choice28    ==j & empFT    ==0) = 28+j;
    choice56_lag(choice28_lag==j & empFT_lag==1) = j; 
    choice56_lag(choice28_lag==j & empFT_lag==0) = 28+j;
end

choice58     = nan(size(choice29));
choice58_lag = nan(size(choice29_lag));
for j=1:29
    choice58    (choice29    ==j & empFT    ==1) = j;
    choice58    (choice29    ==j & empFT    ==0) = 29+j;
    choice58_lag(choice29_lag==j & empFT_lag==1) = j; 
    choice58_lag(choice29_lag==j & empFT_lag==0) = 29+j;
end

choice96     = nan(size(choice48));
choice96_lag = nan(size(choice48_lag));
for j=1:48
    choice96    (choice48    ==j & empFT    ==1) = j;
    choice96    (choice48    ==j & empFT    ==0) = 48+j;
    choice96_lag(choice48_lag==j & empFT_lag==1) = j; 
    choice96_lag(choice48_lag==j & empFT_lag==0) = 48+j;
end

choice98     = nan(size(choice49));
choice98_lag = nan(size(choice49_lag));
for j=1:49
    choice98    (choice49    ==j & empFT    ==1) = j;
    choice98    (choice49    ==j & empFT    ==0) = 49+j;
    choice98_lag(choice49_lag==j & empFT_lag==1) = j; 
    choice98_lag(choice49_lag==j & empFT_lag==0) = 49+j;
end

choice100     = nan(size(choice50));
choice100_lag = nan(size(choice50_lag));
for j=1:50
    choice100    (choice50    ==j & empFT    ==1) = j;
    choice100    (choice50    ==j & empFT    ==0) = 50+j;
    choice100_lag(choice50_lag==j & empFT_lag==1) = j; 
    choice100_lag(choice50_lag==j & empFT_lag==0) = 50+j;
end

choice106     = nan(size(choice53));
choice106_lag = nan(size(choice53_lag));
for j=1:53
    choice106    (choice53    ==j & empFT    ==1) = j;
    choice106    (choice53    ==j & empFT    ==0) = 53+j;
    choice106_lag(choice53_lag==j & empFT_lag==1) = j; 
    choice106_lag(choice53_lag==j & empFT_lag==0) = 53+j;
end

choice108     = nan(size(choice54b));
choice108_lag = nan(size(choice54b_lag));
for j=1:54
    choice108    (choice54b    ==j & empFT    ==1) = j;
    choice108    (choice54b    ==j & empFT    ==0) = 54+j;
    choice108_lag(choice54b_lag==j & empFT_lag==1) = j; 
    choice108_lag(choice54b_lag==j & empFT_lag==0) = 54+j;
end

choice110     = nan(size(choice55));
choice110_lag = nan(size(choice55_lag));
for j=1:55
    choice110    (choice55    ==j & empFT    ==1) = j;
    choice110    (choice55    ==j & empFT    ==0) = 55+j;
    choice110_lag(choice55_lag==j & empFT_lag==1) = j; 
    choice110_lag(choice55_lag==j & empFT_lag==0) = 55+j;
end

choicefrict6     = nan(size(pop_cat));
choicefrict6_lag = nan(size(pop_cat_lag));
for j=1:3
    choicefrict6    (pop_cat    ==j & inlf    ==1) = j;
    choicefrict6    (pop_cat    ==j & inlf    ==0) = 3+j;
    choicefrict6_lag(pop_cat_lag==j & inlf_lag==1) = j; 
    choicefrict6_lag(pop_cat_lag==j & inlf_lag==0) = 3+j;
end

choicefrict54     = nan(size(choice27));
choicefrict54_lag = nan(size(choice27_lag));
for j=1:27
    choicefrict54    (choice27    ==j & inlf    ==1) = j;
    choicefrict54    (choice27    ==j & inlf    ==0) = 27+j;
    choicefrict54_lag(choice27_lag==j & inlf_lag==1) = j; 
    choicefrict54_lag(choice27_lag==j & inlf_lag==0) = 27+j;
end

choicefrict56     = nan(size(choice28));
choicefrict56_lag = nan(size(choice28_lag));
for j=1:28
    choicefrict56    (choice28    ==j & inlf    ==1) = j;
    choicefrict56    (choice28    ==j & inlf    ==0) = 28+j;
    choicefrict56_lag(choice28_lag==j & inlf_lag==1) = j; 
    choicefrict56_lag(choice28_lag==j & inlf_lag==0) = 28+j;
end

choicefrict58     = nan(size(choice29));
choicefrict58_lag = nan(size(choice29_lag));
for j=1:29
    choicefrict58    (choice29    ==j & inlf    ==1) = j;
    choicefrict58    (choice29    ==j & inlf    ==0) = 29+j;
    choicefrict58_lag(choice29_lag==j & inlf_lag==1) = j; 
    choicefrict58_lag(choice29_lag==j & inlf_lag==0) = 29+j;
end

choicefrict96     = nan(size(choice48));
choicefrict96_lag = nan(size(choice48_lag));
for j=1:48
    choicefrict96    (choice48    ==j & inlf    ==1) = j;
    choicefrict96    (choice48    ==j & inlf    ==0) = 48+j;
    choicefrict96_lag(choice48_lag==j & inlf_lag==1) = j; 
    choicefrict96_lag(choice48_lag==j & inlf_lag==0) = 48+j;
end

choicefrict98     = nan(size(choice49));
choicefrict98_lag = nan(size(choice49_lag));
for j=1:49
    choicefrict98    (choice49    ==j & inlf    ==1) = j;
    choicefrict98    (choice49    ==j & inlf    ==0) = 49+j;
    choicefrict98_lag(choice49_lag==j & inlf_lag==1) = j; 
    choicefrict98_lag(choice49_lag==j & inlf_lag==0) = 49+j;
end

choicefrict100     = nan(size(choice50));
choicefrict100_lag = nan(size(choice50_lag));
for j=1:50
    choicefrict100    (choice50    ==j & inlf    ==1) = j;
    choicefrict100    (choice50    ==j & inlf    ==0) = 50+j;
    choicefrict100_lag(choice50_lag==j & inlf_lag==1) = j; 
    choicefrict100_lag(choice50_lag==j & inlf_lag==0) = 50+j;
end

choicefrict106     = nan(size(choice53));
choicefrict106_lag = nan(size(choice53_lag));
for j=1:53
    choicefrict106    (choice53    ==j & inlf    ==1) = j;
    choicefrict106    (choice53    ==j & inlf    ==0) = 53+j;
    choicefrict106_lag(choice53_lag==j & inlf_lag==1) = j; 
    choicefrict106_lag(choice53_lag==j & inlf_lag==0) = 53+j;
end

choicefrict108     = nan(size(choice54b));
choicefrict108_lag = nan(size(choice54b_lag));
for j=1:54
    choicefrict108    (choice54b    ==j & inlf    ==1) = j;
    choicefrict108    (choice54b    ==j & inlf    ==0) = 54+j;
    choicefrict108_lag(choice54b_lag==j & inlf_lag==1) = j; 
    choicefrict108_lag(choice54b_lag==j & inlf_lag==0) = 54+j;
end

choicefrict110     = nan(size(choice55));
choicefrict110_lag = nan(size(choice55_lag));
for j=1:55
    choicefrict110    (choice55    ==j & inlf    ==1) = j;
    choicefrict110    (choice55    ==j & inlf    ==0) = 55+j;
    choicefrict110_lag(choice55_lag==j & inlf_lag==1) = j; 
    choicefrict110_lag(choice55_lag==j & inlf_lag==0) = 55+j;
end

choiceflag = anyFlag==0 & ~isnan(choice6) & ~isnan(educlevel) & ~isnan(age);
mnrfit([collgrad(choiceflag==1) age(choiceflag==1) exper(choiceflag==1)],choice6temp(choiceflag==1))
clear choice6temp

%==============================================================================
% Save the data
%==============================================================================
clear K0 K1 KK all_vars all_long chng_vars cons_vars data stat_mat sum_N sum_stats j k t ans
ID                 = squeeze(ID);
panel              = squeeze(panel);
urate              = squeeze(urate);
age                = squeeze(age);
anyFlag            = squeeze(anyFlag);
rhrsweek           = squeeze(rhrsweek);
choice6            = squeeze(choice6);
choice6_lag        = squeeze(choice6_lag);
choicefrict6       = squeeze(choicefrict6);
choicefrict6_lag   = squeeze(choicefrict6_lag);
choice54           = squeeze(choice54);
choice54_lag       = squeeze(choice54_lag);
choicefrict54      = squeeze(choicefrict54);
choicefrict54_lag  = squeeze(choicefrict54_lag);
choice56           = squeeze(choice56);
choice56_lag       = squeeze(choice56_lag);
choicefrict56      = squeeze(choicefrict56);
choicefrict56_lag  = squeeze(choicefrict56_lag);
choice58           = squeeze(choice58);
choice58_lag       = squeeze(choice58_lag);
choicefrict58      = squeeze(choicefrict58);
choicefrict58_lag  = squeeze(choicefrict58_lag);
choice96           = squeeze(choice96);
choice96_lag       = squeeze(choice96_lag);
choicefrict96      = squeeze(choicefrict96);
choicefrict96_lag  = squeeze(choicefrict96_lag);
choice98           = squeeze(choice98);
choice98_lag       = squeeze(choice98_lag);
choicefrict98      = squeeze(choicefrict98);
choicefrict98_lag  = squeeze(choicefrict98_lag);
choice100          = squeeze(choice100);
choice100_lag      = squeeze(choice100_lag);
choicefrict100     = squeeze(choicefrict100);
choicefrict100_lag = squeeze(choicefrict100_lag);
choice106          = squeeze(choice106);
choice106_lag      = squeeze(choice106_lag);
choicefrict106     = squeeze(choicefrict106);
choicefrict106_lag = squeeze(choicefrict106_lag);
choice108          = squeeze(choice108);
choice108_lag      = squeeze(choice108_lag);
choicefrict108     = squeeze(choicefrict108);
choicefrict108_lag = squeeze(choicefrict108_lag);
choice110          = squeeze(choice110);
choice110_lag      = squeeze(choice110_lag);
choicefrict110     = squeeze(choicefrict110);
choicefrict110_lag = squeeze(choicefrict110_lag);
choiceflag         = squeeze(choiceflag);
collgrad           = squeeze(collgrad);
constant           = squeeze(constant);
educlevel          = squeeze(educlevel);
empFT              = squeeze(empFT);
empFT_lag          = squeeze(empFT_lag);
emp                = squeeze(emp);
emp_lag            = squeeze(emp_lag);
exper_survey       = squeeze(exper);     % SIPP experience (# years with at least two quarters of work)
exper_IRS          = squeeze(experAlt);  % IRS experience (in years)
experFTorPT        = squeeze(experFTorPT);
hgc                = squeeze(hgc);
inlf               = squeeze(inlf);
inlf_lag           = squeeze(inlf_lag);
lnWageHr           = squeeze(lnWageHr);
lnWageHrJ          = squeeze(lnWageHrJ);
lnWageHrJb         = squeeze(lnWageHrJb);
lnearnfinal        = squeeze(lnearnfinal);
lnearnfinalJ       = squeeze(lnearnfinalJ);
lnearnfinalJb      = squeeze(lnearnfinalJb);
pop_cat            = squeeze(pop_cat);
pop_cat_lag        = squeeze(pop_cat_lag);
earnflag           = squeeze(earnflag);
wageflag           = squeeze(wageflag);
weightlong         = squeeze(weightlong);
calmo              = squeeze(calmo);
calyr              = squeeze(calyr);

clear choice27* choice28* choice29* choice48* choice49* choice50* choice53* choice54b* choice55* exper experAlt

whos

% read-in birth location data created in panel_stacker.do
! test -f [REDACTED]sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv.gz && gunzip -f [REDACTED]sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv.gz || echo "sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv.gz Not Found"
! head -1 [REDACTED]sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv
B = importdata('[REDACTED]sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv');
! gzip -f [REDACTED]sippCombinedNHWmaleBirthOnly_matlab_wide_annual.csv
data = B.data;
clear B
[N KK] = size(data)
T  =   5; % maximum number of time periods in panel
K0 =   2; % number of time-invariant variables
K1 = 220; % number of time-varying variables

%==============================================================================
% Read in the data. Columns 1:K0 are time-invariant variables; the rest are
% time-varying regressors.
% Create a master NxKxT tensor that stores all the K variables for all N people
% in each of the T periods.
%==============================================================================
cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),N,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

IDo         = all_vars(:,  1    ,:);
panelo      = all_vars(:,  2    ,:); 
birthLoc110 = all_vars(:,  3:112,:); 
birthDiv110 = all_vars(:,113:222,:);
birthLoc = reshape(permute(birthLoc110,[1 3 2]),N*T,110);
birthDiv = reshape(permute(birthDiv110,[1 3 2]),N*T,110);

% read-in city-level data created in city_level_data.do
% 3 locations:
B = importdata('[REDACTED]urate3.csv');
data = B.data;
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 1; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

pop_cat_CL3      = all_vars(:,1 ,:);
urate3           = all_vars(:,2 ,:);
urate_lag3       = all_vars(:,3 ,:);


% 27 locations:
B = importdata('[REDACTED]urate27.csv');
data = B.data;
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex27        = all_vars(:,1,:);
lnpop27           = all_vars(:,2,:);
urate27           = all_vars(:,3,:);
urate_lag27       = all_vars(:,4,:);


% 28 locations:
B = importdata('[REDACTED]urate28.csv');
data = B.data;
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex28        = all_vars(:,1,:);
lnpop28           = all_vars(:,2,:);
urate28           = all_vars(:,3,:);
urate_lag28       = all_vars(:,4,:);


% 29 locations:
B = importdata('[REDACTED]urate29.csv');
data = B.data;
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex29        = all_vars(:,1,:);
lnpop29           = all_vars(:,2,:);
urate29           = all_vars(:,3,:);
urate_lag29       = all_vars(:,4,:);


% 48 locations:
B = importdata('[REDACTED]urate48.csv');
data = B.data;
if size(B.textdata,1)>1
	data = [[1:48]' data];
end
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex48        = all_vars(:,1,:);
lnpop48           = all_vars(:,2,:);
urate48           = all_vars(:,3,:);
urate_lag48       = all_vars(:,4,:);


% 49 locations:
B = importdata('[REDACTED]urate49.csv');
data = B.data;
if size(B.textdata,1)>1
	data = [[1:49]' data];
end
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex49        = all_vars(:,1,:);
lnpop49           = all_vars(:,2,:);
urate49           = all_vars(:,3,:);
urate_lag49       = all_vars(:,4,:);


% 50 locations:
B = importdata('[REDACTED]urate50.csv');
data = B.data;
if size(B.textdata,1)>1
	data = [[1:50]' data];
end
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex50        = all_vars(:,1,:);
lnpop50           = all_vars(:,2,:);
urate50           = all_vars(:,3,:);
urate_lag50       = all_vars(:,4,:);


% 53 locations:
B = importdata('[REDACTED]urate53.csv');
data = B.data;
if size(B.textdata,1)>1
	data = [[1:53]' data];
end
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex53        = all_vars(:,1,:);
lnpop53           = all_vars(:,2,:);
urate53           = all_vars(:,3,:);
urate_lag53       = all_vars(:,4,:);


% 54 locations:
B = importdata('[REDACTED]urate54.csv');
data = B.data;
if size(B.textdata,1)>1
	data = [[1:54]' data];
end
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex54        = all_vars(:,1,:);
lnpop54           = all_vars(:,2,:);
urate54           = all_vars(:,3,:);
urate_lag54       = all_vars(:,4,:);


% 55 locations:
B = importdata('[REDACTED]urate55.csv');
data = B.data;
if size(B.textdata,1)>1
	data = [[1:55]' data];
end
clear B
[J LL] = size(data)
T  = 10; % maximum number of time periods in panel
K0 = 2; % number of time-invariant variables
K1 = 2; % number of time-varying variables

cons_vars = data(:,1:K0);
chng_vars = reshape(data(:,K0+1:end),J,K1,T);
all_vars  = [cons_vars(:,:,ones(1,T)) chng_vars];

locindex55        = all_vars(:,1,:);
lnpop55           = all_vars(:,2,:);
urate55           = all_vars(:,3,:);
urate_lag55       = all_vars(:,4,:);



year = kron(ones(N,1),[2004 2005 2006 2007 2008 2009 2010 2011 2012 2013]);

urate3      = squeeze(urate3);
urate27     = squeeze(urate27);
urate28     = squeeze(urate28);
urate29     = squeeze(urate29);
urate48     = squeeze(urate48);
urate49     = squeeze(urate49);
urate50     = squeeze(urate50);
urate53     = squeeze(urate53);
urate54     = squeeze(urate54);
urate55     = squeeze(urate55);
urate_lag3  = squeeze(urate_lag3);
urate_lag27 = squeeze(urate_lag27);
urate_lag28 = squeeze(urate_lag28);
urate_lag29 = squeeze(urate_lag29);
urate_lag48 = squeeze(urate_lag48);
urate_lag49 = squeeze(urate_lag49);
urate_lag50 = squeeze(urate_lag50);
urate_lag53 = squeeze(urate_lag53);
urate_lag54 = squeeze(urate_lag54);
urate_lag55 = squeeze(urate_lag55);
lnpop27     = squeeze(lnpop27);
lnpop28     = squeeze(lnpop28);
lnpop29     = squeeze(lnpop29);
lnpop48     = squeeze(lnpop48);
lnpop49     = squeeze(lnpop49);
lnpop50     = squeeze(lnpop50);
lnpop53     = squeeze(lnpop53);
lnpop54     = squeeze(lnpop54);
lnpop55     = squeeze(lnpop55);

clear locindex* pop_cat_CL3 data all_vars J K0 K1 LL chng_vars cons_vars

% read-in distance data created in cr_vincentys.do
% 3 locations:
B = importdata('[REDACTED]dist3.csv');
dist3 = B.data;
dist6 = kron(ones(2),dist3);

% 27 locations:
B = importdata('[REDACTED]dist27.csv');
dist27 = B.data;
dist54 = kron(ones(2),dist27);

% 28 locations:
B = importdata('[REDACTED]dist28.csv');
dist28 = B.data;
dist56 = kron(ones(2),dist28);

% 29 locations:
B = importdata('[REDACTED]dist29.csv');
dist29 = B.data;
dist58 = kron(ones(2),dist29);

% 48 locations:
B = importdata('[REDACTED]dist48.csv');
dist48 = B.data;
dist96 = kron(ones(2),dist48);

% 49 locations:
B = importdata('[REDACTED]dist49.csv');
dist49 = B.data;
dist98 = kron(ones(2),dist49);

% 50 locations:
B = importdata('[REDACTED]dist50.csv');
dist50  = B.data;
dist100 = kron(ones(2),dist50);

% 53 locations:
B = importdata('[REDACTED]dist53.csv');
dist53  = B.data;
dist106 = kron(ones(2),dist53);

% 54 locations:
B = importdata('[REDACTED]dist54.csv');
dist54b  = B.data;
dist108 = kron(ones(2),dist54b);

% 55 locations:
B = importdata('[REDACTED]dist55.csv');
dist55  = B.data;
dist110 = kron(ones(2),dist55);

clear B dist3 dist27 dist28 dist29 dist48 dist49 dist50 dist53 dist54b dist55

whos

save [REDACTED]sippCombinedNHWmaleAnnual.mat
% no need to compress .mat files because they already come that way

toc
diary off;
