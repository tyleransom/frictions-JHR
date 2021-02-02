if strcmp(time,'trimest')==1
	error('You shouldn''t be using trimesterly data!!');
    Beta = Beta.^(1/3); % convert annual discount factor to trimesterly
end

if Beta==0
    delete(['modelFitFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['modelFitFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
else                               
    delete(['modelFitFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['modelFitFrict_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
end

tic;

%==========================================================================
% Read in the data
%==========================================================================
load [REDACTED]sippCombinedNHWmaleAnnual.mat
%==========================================================================
% Read in the parameter estimates
%==========================================================================
if Beta==0
    load(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'Beta0.mat'],'bstruc','YS','ZS','Adj','pimatS','flagSl');
else
    load('allBefDebuggingCCPs.mat','experS','empFT_lagS','inlf_lagS','empFTS','inlfS','nloc','N','T','S','J','dat','bwage','wageflagStemp','calyrS','calyr','calmo','LYS','ageS','lambda_e_hat','lambda_hat','delta_hat','PTypel','baseAlt','wHat','prior','bARw','rho_hat_urate','sig_hat_urate','bLemp','bLunemp','urate_lag55','pimatS','YS');
    load(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'.mat'],'bstruc','YS','ZS','Adj','pimatS','flagSl');
    calmoS       = repmat(calmo,[1 1 S]);         %kron(ones(S,1),calyrS);
end

%==========================================================================
% Preamble
%==========================================================================
seed = 1234;
rng(seed,'twister');

if strcmp(time,'annual')==1
    annual=1;
elseif strcmp(time,'trimest')==1
    annual=0;
end

% code Jan 2009 as Dec 2008 so that 2008 isn't just one rotation group
calyr(calmo==1 & calyr==2009) = 2008;
calmo(calmo==1 & calyr==2008) = 12;
calyrS(calmoS==1 & calyrS==2009) = 2008;
calmoS(calmoS==1 & calyrS==2008) = 12;

expermeanUtemp = wmean(experS(flagSl & empFT_lagS(:)==0 & ageS(:)==25),PTypel(flagSl & empFT_lagS(:)==0 & ageS(:)==25));
expermeanEtemp = wmean(experS(flagSl & empFT_lagS(:)==1 & ageS(:)==25),PTypel(flagSl & empFT_lagS(:)==0 & ageS(:)==25));

for s=1:S
	for k=1:nloc
		for t=1:size(urate55,2);
		Tmattemp = zeros(1,size(urate55,2));
		Tvectemp = t;
		Tmattemp(1,t) = 1;
			deltabar(k,t,s)     = 1-lambdervecRhoS(nloc,k,Tmattemp,Tvectemp,urate_lag55,1,1,expermeanEtemp,s,rho_hat_urate,0,bLemp,bLunemp);
			lambda_e_bar(k,t,s) =   lambdervecRhoS(nloc,k,Tmattemp,Tvectemp,urate_lag55,1,0,expermeanUtemp,s,rho_hat_urate,0,bLemp,bLunemp);
			lambdabar(k,t,s)    =   lambdervecRhoS(nloc,k,Tmattemp,Tvectemp,urate_lag55,0,1,expermeanUtemp,s,rho_hat_urate,0,bLemp,bLunemp);
			lambda_u_bar(k,t,s) =   lambdervecRhoS(nloc,k,Tmattemp,Tvectemp,urate_lag55,0,0,expermeanUtemp,s,rho_hat_urate,0,bLemp,bLunemp);
		end
	end
end

for s=1:S
	summarize(deltabar(:,4,s));
	summarize(lambda_e_bar(:,4,s));
	summarize(lambdabar(:,4,s));
	summarize(lambda_u_bar(:,4,s));
	summarize(deltabar(:,8,s));
	summarize(lambda_e_bar(:,8,s));
	summarize(lambdabar(:,8,s));
	summarize(lambda_u_bar(:,8,s));
end

lambdabar = prior(1)*lambdabar(:,:,1)+prior(2)*lambdabar(:,:,2);
deltabar  = prior(1)*deltabar(:,:,1) +prior(2)*deltabar(:,:,2);

save LLCresults wHat delta_hat lambda_hat deltabar lambdabar bARw rho_hat_urate sig_hat_urate

%==========================================================================
% Employment summary stats
%==========================================================================
moverS = (YS-nloc*(YS>nloc))~=(LYS-nloc*(LYS>nloc));

% employed stayers
options=struct('Weights',PTypel(flagSl(:) & empFT_lagS(:)==1 & moverS==0));
summarize(pimatS(flagSl(:) & empFT_lagS(:)==1 & moverS==0,:),options);

% employed movers

% non-employed stayers

% non-employed movers


%==========================================================================
% Model Fit --- wages
%==========================================================================
predWage = dat.Xw*bwage;

for yy=2004:2013
    tt=yy-2003;
    wagecomp(tt,1) = wmean(dat.wage(wageflagStemp & calyrS==yy),PTypel(wageflagStemp & calyrS==yy));
    wagecomp(tt,2) = wmean(predWage(wageflagStemp & calyrS==yy),PTypel(wageflagStemp & calyrS==yy));
    wagecomp(tt,3) = std (dat.wage(wageflagStemp & calyrS==yy),PTypel(wageflagStemp & calyrS==yy));
    wagecomp(tt,4) = std (predWage(wageflagStemp & calyrS==yy),PTypel(wageflagStemp & calyrS==yy));
end
% for yy=2004:2013
%     tt=yy-2003;
%     wagecompStationary(tt,1) = wmean(wagevar           (wageflagtemp & calyr==yy));
%     wagecompStationary(tt,2) = wmean(predWageStationary(wageflagtemp & calyr==yy));
%     wagecompStationary(tt,3) = std (wagevar           (wageflagtemp & calyr==yy));
%     wagecompStationary(tt,4) = std (predWageStationary(wageflagtemp & calyr==yy));
% end

% Plot model fit over time for mean wage
% plot(wagecomp(:,1)-wagecomp(:,1),'b-')
% hold on
% plot(wagecomp(:,2)-wagecomp(:,1),'r:')
% plot(wagecompStationary(:,2)-wagecomp(:,1),'r-')
% hold off

% Plot model fit over time for SD of wage
% plot(wagecomp(:,3)-wagecomp(:,3),'b-')
% hold on
% plot(wagecomp(:,4)-wagecomp(:,3),'r:')
% plot(wagecompStationary(:,4)-wagecomp(:,3),'r-')
% hold off

%==========================================================================
% Baseline choice probabilities
%==========================================================================
Pstruc = nan(N*T*S,110);
Pstruc(flagSl,:) = pclogitAdj1(bstruc,YS(flagSl),[],ZS(flagSl,:,:),baseAlt,Adj(flagSl,:),1);

%==========================================================================
% Model Fit --- unemployment rate, LFP rate
%==========================================================================
for yy=2004:2013
    tt=yy-2003;
    uratecomp(tt,1) = wmean(inlfS(flagSl & calyrS(:)==yy) & ~empFTS(flagSl & calyrS(:)==yy),PTypel(flagSl & calyrS(:)==yy));
    uratecomp(tt,2) = wmean(sum(Pstruc(flagSl & calyrS(:)==yy,1:nloc).*(1-pimatS(flagSl & calyrS(:)==yy,:)),2),PTypel(flagSl & calyrS(:)==yy));
end
% Plot model fit over time for urate
% plot(uratecomp(:,1),'b-')
% hold on
% plot(uratecomp(:,2),'b:')
% % plot(uratecompStationary(:,4)-uratecomp(:,3),'r-')
% hold off

for yy=2004:2013
    tt=yy-2003;
    LFPratecomp(tt,1) = wmean(inlfS(flagSl & calyrS(:)==yy),PTypel(flagSl & calyrS(:)==yy));
    LFPratecomp(tt,2) = wmean(sum(Pstruc(flagSl & calyrS(:)==yy,1:nloc),2),PTypel(flagSl & calyrS(:)==yy));
end

% Plot model fit over time for LFP rate
% plot(LFPratecomp(:,1),'b-')
% hold on
% plot(LFPratecomp(:,2),'b:')
% % plot(LFPratecompStationary(:,4),'r-')
% hold off

%==========================================================================
% Model Fit --- migration rate
%==========================================================================
location = YS-nloc*(YS>nloc);
location_lag = LYS-nloc*(LYS>nloc);
distvec = nan(size(LYS));
%--------------------------------------------------------------------------
% Baseline migration rates by employment status
%--------------------------------------------------------------------------
Pmig = nan(N*T*S,1);
for i=1:N*T*S
    if flagSl(i)==1
        Pmig(i) = 1-sum(Pstruc(i,[location_lag(i) location_lag(i)+nloc]),2);
    end
end
% for yy=2004:2013
%     tt=yy-2003;
%     migratecomp(tt,1) = wmean(location(flagSl & calyrS(:)==yy) ~= location_lag(flagSl & calyrS(:)==yy));
%     migratecomp(tt,2) = wmean(Pmig(flagSl & calyrS(:)==yy));
%     migratecomp(tt,3) = wmean(location(flagSl & calyrS(:)==yy & empFT_lagS(:)==1) ~= location_lag(flagSl & calyrS(:)==yy & empFT_lagS(:)==1));
%     migratecomp(tt,4) = wmean(Pmig(flagSl & calyrS(:)==yy & empFT_lagS(:)==1));
%     migratecomp(tt,5) = wmean(location(flagSl & calyrS(:)==yy & empFT_lagS(:)==0) ~= location_lag(flagSl & calyrS(:)==yy & empFT_lagS(:)==0));
%     migratecomp(tt,6) = wmean(Pmig(flagSl & calyrS(:)==yy & empFT_lagS(:)==0));
% end
migratecomp(1,1) = wmean(location(flagSl) ~= location_lag(flagSl),PTypel(flagSl));
migratecomp(1,2) = wmean(Pmig(flagSl),PTypel(flagSl));
migratecomp(2,1) = wmean(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==1));
migratecomp(2,2) = wmean(Pmig(flagSl & empFT_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==1));
migratecomp(3,1) = wmean(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecomp(3,2) = wmean(Pmig(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecomp(4,1) = wmean(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0),PTypel(flagSl & inlf_lagS(:)==0));
migratecomp(4,2) = wmean(Pmig(flagSl & inlf_lagS(:)==0),PTypel(flagSl & inlf_lagS(:)==0));
Nmigratecomp(1,1) = length(location(flagSl) ~= location_lag(flagSl))./S;
Nmigratecomp(1,2) = length(Pmig(flagSl))./S;
Nmigratecomp(2,1) = length(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1))./S;
Nmigratecomp(2,2) = length(Pmig(flagSl & empFT_lagS(:)==1))./S;
Nmigratecomp(3,1) = length(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
Nmigratecomp(3,2) = length(Pmig(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
Nmigratecomp(4,1) = length(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0))./S;
Nmigratecomp(4,2) = length(Pmig(flagSl & inlf_lagS(:)==0))./S;

%--------------------------------------------------------------------------
% Baseline migration rates by employment status and pre/post recession
%--------------------------------------------------------------------------
migratecompBC(1,1) = wmean(location(ismember(calyrS(:),[2004:2008]) & flagSl) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl));
migratecompBC(1,2) = wmean(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl));
migratecompBC(2,1) = wmean(location(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1));
migratecompBC(2,2) = wmean(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1));
migratecompBC(3,1) = wmean(location(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompBC(3,2) = wmean(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompBC(4,1) = wmean(location(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0));
migratecompBC(4,2) = wmean(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0),PTypel(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0));
NmigratecompBC(1,1) = length(location(ismember(calyrS(:),[2004:2008]) & flagSl) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl))./S;
NmigratecompBC(1,2) = length(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl))./S;
NmigratecompBC(2,1) = length(location(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1))./S;
NmigratecompBC(2,2) = length(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==1))./S;
NmigratecompBC(3,1) = length(location(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
NmigratecompBC(3,2) = length(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
NmigratecompBC(4,1) = length(location(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0) ~= location_lag(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0))./S;
NmigratecompBC(4,2) = length(Pmig(ismember(calyrS(:),[2004:2008]) & flagSl & inlf_lagS(:)==0))./S;
migratecompBC(1,3) = wmean(location(ismember(calyrS(:),[2009:2013]) & flagSl) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl));
migratecompBC(1,4) = wmean(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl));
migratecompBC(2,3) = wmean(location(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1));
migratecompBC(2,4) = wmean(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1));
migratecompBC(3,3) = wmean(location(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompBC(3,4) = wmean(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompBC(4,3) = wmean(location(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0));
migratecompBC(4,4) = wmean(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0),PTypel(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0));
NmigratecompBC(1,3) = length(location(ismember(calyrS(:),[2009:2013]) & flagSl) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl))./S;
NmigratecompBC(1,4) = length(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl))./S;
NmigratecompBC(2,3) = length(location(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1))./S;
NmigratecompBC(2,4) = length(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==1))./S;
NmigratecompBC(3,3) = length(location(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
NmigratecompBC(3,4) = length(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
NmigratecompBC(4,3) = length(location(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0) ~= location_lag(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0))./S;
NmigratecompBC(4,4) = length(Pmig(ismember(calyrS(:),[2009:2013]) & flagSl & inlf_lagS(:)==0))./S;

%--------------------------------------------------------------------------
% Migration rates by age and employment status
%--------------------------------------------------------------------------
ageRanger = [18 25;26 35;36 45;46 60];
migratecompAge = nan(size(ageRanger,1),8);
NmigratecompAge = nan(size(ageRanger,1),8);
for j=1:size(ageRanger,1);
    migratecompAge(j,1) = wmean(location(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    migratecompAge(j,2) = wmean(Pmig(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    migratecompAge(j,3) = wmean(location(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    migratecompAge(j,4) = wmean(Pmig(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    migratecompAge(j,5) = wmean(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    migratecompAge(j,6) = wmean(Pmig(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    migratecompAge(j,7) = wmean(location(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    migratecompAge(j,8) = wmean(Pmig(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)),PTypel(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)));
    NmigratecompAge(j,1) = length(location(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
    NmigratecompAge(j,2) = length(Pmig(flagSl & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
    NmigratecompAge(j,3) = length(location(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
    NmigratecompAge(j,4) = length(Pmig(flagSl & empFT_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
    NmigratecompAge(j,5) = length(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
    NmigratecompAge(j,6) = length(Pmig(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
    NmigratecompAge(j,7) = length(location(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)) ~= location_lag(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
    NmigratecompAge(j,8) = length(Pmig(flagSl & inlf_lagS(:)==0 & ageS(:)>=ageRanger(j,1) & ageS(:)<=ageRanger(j,2)))./S;
end

%--------------------------------------------------------------------------
% Migration rates by distance and employment status
%--------------------------------------------------------------------------
distRanger = [1 500;500 1000;1000 1500;1500 2000;2000 5500];
PmigDist= nan(N*T*S,size(distRanger,1));
for i=1:N*T*S
    if flagSl(i)==1
        distvec(i) = dist110(location(i),location_lag(i));
        for j=1:size(distRanger,1)
            PmigDist(i,j) = sum(Pstruc(i,:).*(dist110(location_lag(i),:)>=distRanger(j,1) & dist110(location_lag(i),:)<=distRanger(j,2)),2);
        end
    end
end
migratecompDist = nan(size(distRanger,1),8);
NmigratecompDist = nan(size(distRanger,1),8);
for j=1:size(distRanger,1);
    migratecompDist(j,1) = wmean(location(flagSl) ~= location_lag(flagSl) & distvec(flagSl)>=distRanger(j,1) & distvec(flagSl)<=distRanger(j,2),PTypel(flagSl));
    migratecompDist(j,2) = wmean(PmigDist(flagSl,j),PTypel(flagSl));
    migratecompDist(j,3) = wmean(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & distvec(flagSl & empFT_lagS(:)==1)>=distRanger(j,1) & distvec(flagSl & empFT_lagS(:)==1)<=distRanger(j,2),PTypel(flagSl & empFT_lagS(:)==1));
    migratecompDist(j,4) = wmean(PmigDist(flagSl & empFT_lagS(:)==1,j),PTypel(flagSl & empFT_lagS(:)==1));
    migratecompDist(j,5) = wmean(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) & distvec(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)>=distRanger(j,1) & distvec(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)<=distRanger(j,2),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
    migratecompDist(j,6) = wmean(PmigDist(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1,j),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
    migratecompDist(j,7) = wmean(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0) & distvec(flagSl & inlf_lagS(:)==0)>=distRanger(j,1) & distvec(flagSl & inlf_lagS(:)==0)<=distRanger(j,2),PTypel(flagSl & inlf_lagS(:)==0));
    migratecompDist(j,8) = wmean(PmigDist(flagSl & inlf_lagS(:)==0,j),PTypel(flagSl & inlf_lagS(:)==0));
    NmigratecompDist(j,1) = length(location(flagSl) ~= location_lag(flagSl) & distvec(flagSl)>=distRanger(j,1) & distvec(flagSl)<=distRanger(j,2))./S;
    NmigratecompDist(j,2) = length(PmigDist(flagSl,j))./S;
    NmigratecompDist(j,3) = length(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & distvec(flagSl & empFT_lagS(:)==1)>=distRanger(j,1) & distvec(flagSl & empFT_lagS(:)==1)<=distRanger(j,2))./S;
    NmigratecompDist(j,4) = length(PmigDist(flagSl & empFT_lagS(:)==1,j))./S;
    NmigratecompDist(j,5) = length(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) & distvec(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)>=distRanger(j,1) & distvec(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)<=distRanger(j,2))./S;
    NmigratecompDist(j,6) = length(PmigDist(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1,j))./S;
    NmigratecompDist(j,7) = length(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0) & distvec(flagSl & inlf_lagS(:)==0)>=distRanger(j,1) & distvec(flagSl & inlf_lagS(:)==0)<=distRanger(j,2))./S;
    NmigratecompDist(j,8) = length(PmigDist(flagSl & inlf_lagS(:)==0,j))./S;
end

%--------------------------------------------------------------------------
% Migration rates by employment status: Movement to better-employment areas
%--------------------------------------------------------------------------
PmigPiE=nan(N*T*S,1);
PmigPiN=nan(N*T*S,1);
comp1 = nan(N*T*S,J);
comp2 = nan(N*T*S,J);
comp1data = nan(N*T*S,1);
comp2data = nan(N*T*S,1);
for i=1:N*T*S
    if flagSl(i)==1
        comp1(i,:) = kron(ones(1,2),(lambda_e_hat(:,calyrS(i)-2003)'>lambda_e_hat(location_lag(i),calyrS(i)-2003)));
        comp2(i,:) = kron(ones(1,2),(lambda_hat(:,calyrS(i)-2003)'>lambda_hat(location_lag(i),calyrS(i)-2003)));
        comp1data(i) = lambda_e_hat(location(i),calyrS(i)-2003)'>lambda_e_hat(location_lag(i),calyrS(i)-2003);
        comp2data(i) = lambda_hat(location(i),calyrS(i)-2003)'>lambda_hat(location_lag(i),calyrS(i)-2003);
        PmigPiE(i) = sum(Pstruc(i,:).*comp1(i,:),2);
        PmigPiN(i) = sum(Pstruc(i,:).*comp2(i,:),2);
    end
end
migratecompPi(1,1) = wmean(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==1)==1,PTypel(flagSl & empFT_lagS(:)==1));
migratecompPi(1,2) = wmean(PmigPiE(flagSl & empFT_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==1));
migratecompPi(2,1) = wmean(location(flagSl & empFT_lagS(:)==0) ~= location_lag(flagSl & empFT_lagS(:)==0) & comp2data(flagSl & empFT_lagS(:)==0)==1,PTypel(flagSl & empFT_lagS(:)==0));
migratecompPi(2,2) = wmean(PmigPiN(flagSl & empFT_lagS(:)==0),PTypel(flagSl & empFT_lagS(:)==0));
NmigratecompPi(1,1) = length(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==1)==1)./S;
NmigratecompPi(1,2) = length(PmigPiE(flagSl & empFT_lagS(:)==1))./S;
NmigratecompPi(2,1) = length(location(flagSl & empFT_lagS(:)==0) ~= location_lag(flagSl & empFT_lagS(:)==0) & comp2data(flagSl & empFT_lagS(:)==0)==0)./S;
NmigratecompPi(2,2) = length(PmigPiN(flagSl & empFT_lagS(:)==0))./S;

%--------------------------------------------------------------------------
% Migration rates by employment status: Movement to higher earnings areas
%--------------------------------------------------------------------------
PmigW=nan(N*T*S,1);
comp1 = nan(N*T*S,J);
comp1data = nan(N*T*S,1);
for i=1:N*T*S
    if flagSl(i)==1
        comp1(i,:) = kron(ones(1,2),(wHat(:,calyrS(i)-2003)'>wHat(location_lag(i),calyrS(i)-2003)));
        comp1data(i) = wHat(location(i),calyrS(i)-2003)'>wHat(location_lag(i),calyrS(i)-2003);
        PmigW(i) = sum(Pstruc(i,:).*comp1(i,:),2);
    end
end
migratecompW(1,1) = wmean(location(flagSl) ~= location_lag(flagSl) & comp1data(flagSl)==1,PTypel(flagSl));
migratecompW(1,2) = wmean(PmigW(flagSl),PTypel(flagSl));
migratecompW(2,1) = wmean(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==1)==1,PTypel(flagSl & empFT_lagS(:)==1));
migratecompW(2,2) = wmean(PmigW(flagSl & empFT_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==1));
migratecompW(3,1) = wmean(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)==1,PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompW(3,2) = wmean(PmigW(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompW(4,1) = wmean(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0) & comp1data(flagSl & inlf_lagS(:)==0)==1,PTypel(flagSl & inlf_lagS(:)==0));
migratecompW(4,2) = wmean(PmigW(flagSl & inlf_lagS(:)==0),PTypel(flagSl & inlf_lagS(:)==0));
NmigratecompW(1,1) = length(location(flagSl) ~= location_lag(flagSl) & comp1data(flagSl)==1)./S;
NmigratecompW(1,2) = length(PmigW(flagSl))./S;
NmigratecompW(2,1) = length(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==1)==1)./S;
NmigratecompW(2,2) = length(PmigW(flagSl & empFT_lagS(:)==1))./S;
NmigratecompW(3,1) = length(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)==1)./S;
NmigratecompW(3,2) = length(PmigW(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
NmigratecompW(4,1) = length(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0) & comp1data(flagSl & inlf_lagS(:)==0)==1)./S;
NmigratecompW(4,2) = length(PmigW(flagSl & inlf_lagS(:)==0))./S;

%--------------------------------------------------------------------------
% Migration rates by employment status: Movement to higher amenity areas
%--------------------------------------------------------------------------
%load strucFrict_earn_55locHS bstruc
PmigA=nan(N*T*S,1);
comp1 = nan(N*T*S,J);
comp1data = nan(N*T*S,1);
for i=1:N*T*S
    if flagSl(i)==1
        distvec(i) = dist110(location(i),location_lag(i));
        comp1(i,:) = kron(ones(1,2),(bstruc(1:55)'>bstruc(location_lag(i))));
        comp1data(i) = bstruc(location(i))>bstruc(location_lag(i));
        PmigA(i) = sum(Pstruc(i,:).*comp1(i,:),2);
    end
end
migratecompA(1,1) = wmean(location(flagSl) ~= location_lag(flagSl) & comp1data(flagSl)==1,PTypel(flagSl));
migratecompA(1,2) = wmean(PmigA(flagSl),PTypel(flagSl));
migratecompA(2,1) = wmean(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==1)==1,PTypel(flagSl & empFT_lagS(:)==1));
migratecompA(2,2) = wmean(PmigA(flagSl & empFT_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==1));
migratecompA(3,1) = wmean(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)==1,PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompA(3,2) = wmean(PmigA(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1),PTypel(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1));
migratecompA(4,1) = wmean(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0) & comp1data(flagSl & inlf_lagS(:)==0)==1,PTypel(flagSl & inlf_lagS(:)==0));
migratecompA(4,2) = wmean(PmigA(flagSl & inlf_lagS(:)==0),PTypel(flagSl & inlf_lagS(:)==0));
NmigratecompA(1,1) = length(location(flagSl) ~= location_lag(flagSl) & comp1data(flagSl)==1)./S;
NmigratecompA(1,2) = length(PmigA(flagSl))./S;
NmigratecompA(2,1) = length(location(flagSl & empFT_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==1)==1)./S;
NmigratecompA(2,2) = length(PmigA(flagSl & empFT_lagS(:)==1))./S;
NmigratecompA(3,1) = length(location(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) ~= location_lag(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1) & comp1data(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1)==1)./S;
NmigratecompA(3,2) = length(PmigA(flagSl & empFT_lagS(:)==0 & inlf_lagS(:)==1))./S;
NmigratecompA(4,1) = length(location(flagSl & inlf_lagS(:)==0) ~= location_lag(flagSl & inlf_lagS(:)==0) & comp1data(flagSl & inlf_lagS(:)==0)==1)./S;
NmigratecompA(4,2) = length(PmigA(flagSl & inlf_lagS(:)==0))./S;


wagecomp
uratecomp
LFPratecomp
migratecomp
Nmigratecomp
migratecompBC
NmigratecompBC
migratecompAge
NmigratecompAge
migratecompDist
NmigratecompDist
migratecompPi
NmigratecompPi
migratecompW
NmigratecompW
migratecompA
NmigratecompA





%==========================================================================
% Model Fit --- employment transitions
%==========================================================================
%--------------------------------------------------------------------------
% Baseline migration rates by employment status
%--------------------------------------------------------------------------
Pmig2     = nan(N*T*S,1);
PstayLF   = nan(N*T*S,1);
PstayNILF = nan(N*T*S,1);
PmigLF    = nan(N*T*S,1);
PmigNILF  = nan(N*T*S,1);
piStay    = nan(N*T*S,1);
PmigEmp   = nan(N*T*S,1);
PmigUnemp = nan(N*T*S,1);
tranMat   = nan(3,3,4);
NtranMat  = nan(3,3,4);
for i=1:N*T*S
    if flagSl(i)==1
        PstayLF(i)   = Pstruc(i,location_lag(i));
        PstayNILF(i) = Pstruc(i,location_lag(i)+nloc);
        PmigLF(i)    = sum(Pstruc(i,1:nloc),2)-PstayLF(i);
        PmigNILF(i)  = sum(Pstruc(i,nloc+1:end),2)-PstayNILF(i);
        Pmig2(i)     = sum(Pstruc(i,setdiff(1:2*nloc,[location_lag(i) location_lag(i)+nloc])),2);
        piStay(i)    = pimatS(i,location_lag(i));
        PmigEmp(i)   = sum(Pstruc(i,1:nloc).*pimatS(i,:),2)-piStay(i).*PstayLF(i);
        PmigUnemp(i) = sum(Pstruc(i,1:nloc).*(1-pimatS(i,:)),2)-(1-piStay(i)).*PstayLF(i);
    end
end
Pmig = PmigLF+PmigNILF;
Pstay= PstayLF+PstayNILF;
laggOutcomeStatus = 1*(empFT_lagS(:)==1)+2*(empFT_lagS(:)==0 & inlf_lagS(:)==1)+3*(inlf_lagS(:)==0);
currOutcomeStatus = 1*(empFTS    (:)==1)+2*(empFTS    (:)==0 & inlfS    (:)==1)+3*(inlfS    (:)==0);

save -v7.3 modelFitTest.mat

% Data, conditional on moving
for j=1:3
	for k=1:3
		subby = flagSl & location~=location_lag & laggOutcomeStatus==j;
		tranMat(j,k,1)  = wmean(laggOutcomeStatus(subby)==j & currOutcomeStatus(subby)==k,PTypel(subby));
		NtranMat(j,k,1) = sum((laggOutcomeStatus(subby)==j & currOutcomeStatus(subby)==k).*PTypel(subby));
	end
end
% Data, conditional on staying
for j=1:3
	for k=1:3
		subby = flagSl & location==location_lag & laggOutcomeStatus==j;
		tranMat(j,k,3) = wmean(laggOutcomeStatus(subby)==j & currOutcomeStatus(subby)==k,PTypel(subby));
		NtranMat(j,k,4) = sum((laggOutcomeStatus(subby)==j & currOutcomeStatus(subby)==k).*PTypel(subby));
	end
end

% Model, condition
for j=1:3
	subby = flagSl & laggOutcomeStatus==j;
	tranMat(j,1,2) = wmean(PmigEmp(subby)./Pmig(subby),PTypel(subby));
	tranMat(j,2,2) = wmean(PmigUnemp(subby)./Pmig(subby),PTypel(subby));
	tranMat(j,3,2) = wmean(PmigNILF(subby)./Pmig(subby),PTypel(subby));
end
for j=1:3
	subby = flagSl & laggOutcomeStatus==j;
	tranMat(j,1,4) = wmean((PstayLF(subby).*piStay(subby))./Pstay(subby),PTypel(subby));
	tranMat(j,2,4) = wmean((PstayLF(subby).*(1-piStay(subby)))./Pstay(subby),PTypel(subby));
	tranMat(j,3,4) = wmean(PstayNILF(subby)./Pstay(subby),PTypel(subby));
end

tranMat(:,:,1)
tranMat(:,:,2)
tranMat(:,:,3)
tranMat(:,:,4)
NtranMat(:,:,1)
NtranMat(:,:,4)

