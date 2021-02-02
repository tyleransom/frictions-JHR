%**************************************************************************
% Program outline:
% 1. Compute correlations between offers, layoffs, wages, and amenities
%**************************************************************************
%==========================================================================
% Import data
%==========================================================================
clear all;
clc;
addpath [REDACTED]

pathstring = pwd;
pathstring = lower(pathstring);
pathstring = pathstring(end-11:end);

% global Beta
nloc   = str2num(pathstring(1:2));
money  = pathstring(7:10); % enter 'wage' or 'earn'
time   = 'annual'; % enter 'annual' or 'trimes'
Beta     =.9;
% Beta     = 0;
sample = upper(pathstring(end-1:end)); % enter 'HS' or 'BA'

if Beta==0
%    delete(['backOutCorrelationsFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['typeFreqChecker_beta',num2str(10*Beta),'.diary']);
else                               
%    delete(['backOutCorrelationsFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['typeFreqChecker_beta',num2str(10*Beta),'.diary']);
end

load allBefDebuggingCCPs wageflagS ChoiceS flagS empFTS empFT_lagS dat PTypeS IDS inlfS S typeS

%---------------------------------------
% Wages
%---------------------------------------
disp('Wages')
wageflagStempl = wageflagS(:)==1 & ChoiceS(:)<=nloc & flagS(:) & empFTS(:)==1;

% Number of person-years in each type
for s=1:S
	sum(PTypeS(wageflagStempl).*(dat.Xw(wageflagStempl,end)==2-s))
end

% Number of persons in each type
for s=1:S
	numel(unique(IDS(wageflagStempl))).*(sum(PTypeS(wageflagStempl).*(dat.Xw(wageflagStempl,end)==2-s))/sum(PTypeS(wageflagStempl)))
end

%---------------------------------------
% Employment
%---------------------------------------
disp('Emp -> Emp')
conditioner = inlfS(:)==1 & empFT_lagS(:)==1 & flagS(:)==1 & empFTS(:)==1;

% Number of person-years in each type
for s=1:S
	sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))
end

% Number of persons in each type
for s=1:S
	numel(unique(IDS(conditioner))).*(sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))/sum(PTypeS(conditioner)))
end

disp('Emp -> Not Emp')
conditioner = inlfS(:)==1 & empFT_lagS(:)==1 & flagS(:)==1 & empFTS(:)==0;

% Number of person-years in each type
for s=1:S
	sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))
end

% Number of persons in each type
for s=1:S
	numel(unique(IDS(conditioner))).*(sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))/sum(PTypeS(conditioner)))
end


disp('Not Emp -> Emp')
conditioner = inlfS(:)==1 & empFT_lagS(:)==0 & flagS(:)==1 & empFTS(:)==1;

% Number of person-years in each type
for s=1:S
	sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))
end

% Number of persons in each type
for s=1:S
	numel(unique(IDS(conditioner))).*(sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))/sum(PTypeS(conditioner)))
end


disp('Not Emp -> Not Emp')
conditioner = inlfS(:)==1 & empFT_lagS(:)==0 & flagS(:)==1 & empFTS(:)==0;

% Number of person-years in each type
for s=1:S
	sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))
end

% Number of persons in each type
for s=1:S
	numel(unique(IDS(conditioner))).*(sum(PTypeS(conditioner).*(typeS(conditioner)==2-s))/sum(PTypeS(conditioner)))
end

diary('off')

