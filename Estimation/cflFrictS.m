if strcmp(time,'trimest')==1
	error('You shouldn''t be using trimesterly data!!');
    Beta = Beta.^(1/3); % convert annual discount factor to trimesterly
end

if Beta==0
    delete([diarystub,num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary ([diarystub,num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
else                               
    delete([diarystub,num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary ([diarystub,num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
end

tic;


load('allBefDebuggingCCPs.mat','draws','distance','dist110','ARcov','ARcov2','birthLocS','birthDivS','empFT_lagS','inlf_lagS','empFTS','inlfS','flagS','flagSl','nloc','N','T','S','J','dat','bwage','wageflagStemp','calyrS','calyr','calmo','LYS','ageS','experS','lambda_e_hat','lambda_hat','lambda_u_hat','delta_hat','PTypel','baseAlt','wHat','bARw','rho_hat_wage','rho_hat_urate','sig_hat_wage','sig_hat_urate','urate_lag55','bLemp','bLunemp');
load([money,'_',num2str(nloc),'loc',sample,'Het.mat'],'bHet');
load([money,'_',num2str(nloc),'loc',sample,'Het.mat'],'bHet');
load(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'.mat'],'bstruc');

% ========================================================================
% Create "fake" cities with the specified feature profiles: set to be same location!!!
% ========================================================================
disp('wage shock quartiles')
summarize(bARw(1:nloc));
prctile(bARw(1:nloc),50)
summarize(bARw(nloc+2:end));
prctile(bARw(nloc+2:end),50)
summarize(rho_hat_wage');
prctile(rho_hat_wage',50)
summarize(rho_hat_urate');
prctile(rho_hat_urate',50)

% shock volatility, drift, and autocorrelation: set to be the same everywhere
% bARw(1:nloc)     = kron(ones(nloc,1),prctile(bARw(1:nloc),50));
% bARw(nloc+2:end) = kron(ones(nloc,1),prctile(bARw(nloc+2:end),50));
% rho_hat_wage     = kron(ones(nloc,1),prctile(rho_hat_wage',50))';
% rho_hat_urate    = kron(ones(nloc,1),prctile(rho_hat_urate',50))';

disp('amenity quartiles')
prctile(bstruc(1:nloc),75)
prctile(bstruc(1:nloc),25)
disp('amenity FL quartiles')
prctile(bHet(1:2:2*(nloc-1)+1),75)
prctile(bHet(1:2:2*(nloc-1)+1),25)
disp('earnings quartiles')
prctile(wHat(:,2:end),75,1)
prctile(wHat(:,2:end),25,1)
disp('urate quartiles')
prctile(urate_lag55,25,1)
prctile(urate_lag55,75,1)
disp('emp logit parm quartiles')
prctile([0;bLemp(2:nloc)],75,1)
prctile([0;bLunemp(2:nloc)],75,1)


if strcmpi(locstring,'hmm')
	LYtemp = 14;
	% amenities (structural)
    bstruc(LYtemp)                 = prctile(bstruc(1:nloc),75);
    % amenities (flexible logit)
    bHet([2*LYtemp-3 2*(nloc-1)+2*LYtemp-1]) = prctile(bHet(1:2:2*(nloc-1)+1),75);
    % earnings
    bwage(nloc+(LYtemp-1)*9+[1:9]) = prctile(wHat(:,2:end),50,1); % earnings level
    bARw(LYtemp)                   = prctile(bARw(1:nloc),50); % earnings drift
    bARw(nloc+1+LYtemp)            = prctile(bARw(nloc+2:end),50); % earnings shock variance
    % employment
    urate_lag55(LYtemp,:)          = prctile(urate_lag55,50,1); % unemployment level
    bLemp(LYtemp)                  = prctile([0;bLemp(2:nloc)],50,1); % unemployment level
    bLunemp(LYtemp)                = prctile([0;bLunemp(2:nloc)],50,1); % unemployment level
    rho_hat_urate(1,LYtemp)        = prctile(rho_hat_urate(1,:),50); % unemployment drift
    rho_hat_urate(2,LYtemp)        = prctile(rho_hat_urate(2,:),50); % unemployment autocorrelation
    sig_hat_urate(1,LYtemp)        = prctile(sig_hat_urate(1,:),50); % unemployment shock variance
elseif strcmpi(locstring,'lmm')
	LYtemp = 14;
	% amenities (structural)
    bstruc(LYtemp)                 = prctile(bstruc(1:nloc),25);
    % amenities (flexible logit)
    bHet([2*LYtemp-3 2*(nloc-1)+2*LYtemp-1]) = prctile(bHet(1:2:2*(nloc-1)+1),25);
    % earnings
    bwage(nloc+(LYtemp-1)*9+[1:9]) = prctile(wHat(:,2:end),50,1); % earnings level
    bARw(LYtemp)                   = prctile(bARw(1:nloc),50); % earnings drift
    bARw(nloc+1+LYtemp)            = prctile(bARw(nloc+2:end),50); % earnings shock variance
    % employment
    urate_lag55(LYtemp,:)          = prctile(urate_lag55,50,1); % unemployment level
    bLemp(LYtemp)                  = prctile([0;bLemp(2:nloc)],50,1); % unemployment level
    bLunemp(LYtemp)                = prctile([0;bLunemp(2:nloc)],50,1); % unemployment level
    rho_hat_urate(1,LYtemp)        = prctile(rho_hat_urate(1,:),50); % unemployment drift
    rho_hat_urate(2,LYtemp)        = prctile(rho_hat_urate(2,:),50); % unemployment autocorrelation
    sig_hat_urate(1,LYtemp)        = prctile(sig_hat_urate(1,:),50); % unemployment shock variance
elseif strcmpi(locstring,'mhm')
	LYtemp = 14;
	% amenities (structural)
    bstruc(LYtemp)                 = prctile(bstruc(1:nloc),50);
    % amenities (flexible logit)
    bHet([2*LYtemp-3 2*(nloc-1)+2*LYtemp-1]) = prctile(bHet(1:2:2*(nloc-1)+1),50);
    % earnings
    bwage(nloc+(LYtemp-1)*9+[1:9]) = prctile(wHat(:,2:end),75,1); % earnings level
    bARw(LYtemp)                   = prctile(bARw(1:nloc),50); % earnings drift
    bARw(nloc+1+LYtemp)            = prctile(bARw(nloc+2:end),50); % earnings shock variance
    % employment
    urate_lag55(LYtemp,:)          = prctile(urate_lag55,50,1); % unemployment level
    bLemp(LYtemp)                  = prctile([0;bLemp(2:nloc)],50,1); % unemployment level
    bLunemp(LYtemp)                = prctile([0;bLunemp(2:nloc)],50,1); % unemployment level
    rho_hat_urate(1,LYtemp)        = prctile(rho_hat_urate(1,:),50); % unemployment drift
    rho_hat_urate(2,LYtemp)        = prctile(rho_hat_urate(2,:),50); % unemployment autocorrelation
    sig_hat_urate(1,LYtemp)        = prctile(sig_hat_urate(1,:),50); % unemployment shock variance
elseif strcmpi(locstring,'mlm')
	LYtemp = 14;
	% amenities (structural)
    bstruc(LYtemp)                 = prctile(bstruc(1:nloc),50);
    % amenities (flexible logit)
    bHet([2*LYtemp-3 2*(nloc-1)+2*LYtemp-1]) = prctile(bHet(1:2:2*(nloc-1)+1),50);
    % earnings
    bwage(nloc+(LYtemp-1)*9+[1:9]) = prctile(wHat(:,2:end),25,1); % earnings level
    bARw(LYtemp)                   = prctile(bARw(1:nloc),50); % earnings drift
    bARw(nloc+1+LYtemp)            = prctile(bARw(nloc+2:end),50); % earnings shock variance
    % employment
    urate_lag55(LYtemp,:)          = prctile(urate_lag55,50,1); % unemployment level
    bLemp(LYtemp)                  = prctile([0;bLemp(2:nloc)],50,1); % unemployment level
    bLunemp(LYtemp)                = prctile([0;bLunemp(2:nloc)],50,1); % unemployment level
    rho_hat_urate(1,LYtemp)        = prctile(rho_hat_urate(1,:),50); % unemployment drift
    rho_hat_urate(2,LYtemp)        = prctile(rho_hat_urate(2,:),50); % unemployment autocorrelation
    sig_hat_urate(1,LYtemp)        = prctile(sig_hat_urate(1,:),50); % unemployment shock variance
elseif strcmpi(locstring,'mmh')
	LYtemp = 14;
	% amenities (structural)
    bstruc(LYtemp)                 = prctile(bstruc(1:nloc),50);
    % amenities (flexible logit)
    bHet([2*LYtemp-3 2*(nloc-1)+2*LYtemp-1]) = prctile(bHet(1:2:2*(nloc-1)+1),50);
    % earnings
    bwage(nloc+(LYtemp-1)*9+[1:9]) = prctile(wHat(:,2:end),50,1); % earnings level
    bARw(LYtemp)                   = prctile(bARw(1:nloc),50); % earnings drift
    bARw(nloc+1+LYtemp)            = prctile(bARw(nloc+2:end),50); % earnings shock variance
    % employment
    urate_lag55(LYtemp,:)          = prctile(urate_lag55,25,1); % unemployment level
    bLemp(LYtemp)                  = prctile([0;bLemp(2:nloc)],75,1); % unemployment level
    bLunemp(LYtemp)                = prctile([0;bLunemp(2:nloc)],75,1); % unemployment level
    rho_hat_urate(1,LYtemp)        = prctile(rho_hat_urate(1,:),50); % unemployment drift
    rho_hat_urate(2,LYtemp)        = prctile(rho_hat_urate(2,:),50); % unemployment autocorrelation
    sig_hat_urate(1,LYtemp)        = prctile(sig_hat_urate(1,:),50); % unemployment shock variance
elseif strcmpi(locstring,'mml')
	LYtemp = 14;
	% amenities (structural)
    bstruc(LYtemp)                 = prctile(bstruc(1:nloc),50);
    % amenities (flexible logit)
    bHet([2*LYtemp-3 2*(nloc-1)+2*LYtemp-1]) = prctile(bHet(1:2:2*(nloc-1)+1),50);
    % earnings
    bwage(nloc+(LYtemp-1)*9+[1:9]) = prctile(wHat(:,2:end),50,1); % earnings level
    bARw(LYtemp)                   = prctile(bARw(1:nloc),50); % earnings drift
    bARw(nloc+1+LYtemp)            = prctile(bARw(nloc+2:end),50); % earnings shock variance
    % employment
    urate_lag55(LYtemp,:)          = prctile(urate_lag55,75,1); % unemployment level
    bLemp(LYtemp)                  = prctile([0;bLemp(2:nloc)],25,1); % unemployment level
    bLunemp(LYtemp)                = prctile([0;bLunemp(2:nloc)],25,1); % unemployment level
    rho_hat_urate(1,LYtemp)        = prctile(rho_hat_urate(1,:),50); % unemployment drift
    rho_hat_urate(2,LYtemp)        = prctile(rho_hat_urate(2,:),50); % unemployment autocorrelation
    sig_hat_urate(1,LYtemp)        = prctile(sig_hat_urate(1,:),50); % unemployment shock variance
end
% LYtemp = 49;
% yr     = 2004;
% birthLoci =   [ 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0];
% birthDivi =   [ 0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0     0     1     0     0     0     0     0     0     1     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     0     0     0     0     0     0];
% ========================================================================
% Evaluate Adj and flow utilities at certain characteristics (employed)
% ========================================================================
l   = LYtemp;
LEi = 1;
typer = 2-typer
wage_shock_amt =  -1.5*bARw(nloc+1+LYtemp); %1 SD decrease in wages
pi_shock_amt   = 2.5*sig_hat_urate(LYtemp); %1 SD increase in UR
MCbonus = .1; % e.g. .1 means a 10% subsidy to fixed cost of moving
experer= wmean(experS(flagS & ageS==ager & empFTS==LEi),PTypel(flagSl & ageS(:)==ager & empFTS(:)==LEi));
tic
%Tmater = zeros(1,10);
%Tmater(yr-2003)=1;
%AdjCompNew = estAdjFun(bHet,LYtemp,baseAlt,1,1,draws,true,true,LEi,ARcov,ARcov2,J,nloc,Tmater,yr-2003,urate_lag55,experer,rho_hat_urate,rho_hat_wage,bLemp,bLunemp,bwage,birthLoci,birthDivi,ager,distance,Beta);
AdjNew = getAdjIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,distance,baseAlt,Beta);
disp(['getAdj function took ',num2str(toc/60),' minutes']);
assert(abs(AdjNew(LYtemp+nloc))<1e-1,['Problem with getAdj, abs(AdjNew(LYtemp+nloc)) = ',num2str(abs(AdjNew(LYtemp+nloc)))]);
AdjMat(5,:) = AdjNew;

%Tmater = zeros(1,10);
%Tmater(yr-2003)=1;
%pvEmp    = LYtemp<nloc & LEi==1;
%pvUnemp  = LYtemp<nloc & LEi==0;

tic
[ZNew      ,pimatNew] = getZintS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,distance,Beta);
[ZNewStatic]          = getZintS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,distance,0);
% [ZNewfun   ,pimaNfun] = constructZ(1,1,J,nloc,draws,true,true,ARcov,ARcov2,yr,Tmater,yr-2003,urate_lag55,rho_hat_urate,rho_hat_wage,bLemp,bLunemp,bwage,experer,ager,distance,birthLoc,birthDiv,LYtemp,LYtemp,pvEmp,pvUnemp,LEi,Beta);
% [ZNewStaticfun]       = constructZ(1,1,J,nloc,draws,true,true,ARcov,ARcov2,yr,Tmater,yr-2003,urate_lag55,rho_hat_urate,rho_hat_wage,bLemp,bLunemp,bwage,experer,ager,distance,birthLoc,birthDiv,LYtemp,LYtemp,pvEmp,pvUnemp,LEi,0);
disp(['getZ function took ',num2str(toc/60),' minutes']);
disp('Adj baseline summary stats')
summarize(AdjNew(:));
disp('pimat baseline summary stats')
summarize(pimatNew(:));
disp('Z baseline summary stats')
summarize(ZNew(:));

%==========================================================================
% Baseline choice probabilities
%==========================================================================
PstrucNew = pclogitAdj1a(bstruc,l,[],ZNew,J,baseAlt,AdjNew,1);

uratebase   = nan(6,2);
LFPratebase = nan(6,2);
migratebase = nan(6,2);
uratecomp   = nan(6,2);
LFPratecomp = nan(6,2);
migratecomp = nan(6,2);

uratebase(1:6,1)   =   kron(ones(6,1),(PstrucNew(LYtemp).*(1-pimatNew(LYtemp))));
LFPratebase(1:6,1) =   kron(ones(6,1),(PstrucNew(LYtemp))./(sum(PstrucNew([LYtemp LYtemp+nloc]),2)));
LFPratebaseA(1:6,1)=   kron(ones(6,1),(PstrucNew(LYtemp)));
LFPratebaseB(1:6,1)=   kron(ones(6,1),(PstrucNew(LYtemp+nloc)));
migratebase(1:6,1) = 1-kron(ones(6,1),sum(PstrucNew([LYtemp LYtemp+nloc]),2));
inmigratebase(1,1) = sum(PstrucNew([susumeLocs susumeLocs+nloc]),2);

%==========================================================================
% Counterfactual Sim #1: 1 s.d. transitory decrease in wage in current loc
%==========================================================================
AdjaNew = getAdjCflWageIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflWageIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
AdjMat(1,:) = AdjaNew;
disp('Adj summary stats 1 emp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 1 emp FLint')
summarize(pimataNew(:));
disp('Z summary stats 1 emp FLint')
summarize(ZaNew(:));

uratecomp(1,1)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(1,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(1,1) = (PstrucTilde(LYtemp));
LFPratecompB(1,1) = (PstrucTilde(LYtemp+nloc));
migratecomp(1,1) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPWshock(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateWshock(:,1) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigWshock(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #2: 1 s.d. transitory decrease in wage in current loc; correlated across all locations
%==========================================================================
AdjaNew = getAdjCflWageAllIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflWageAllIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
AdjMat(2,:) = AdjaNew;
disp('Adj summary stats 2 emp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 2 emp FLint')
summarize(pimataNew(:));
disp('Z summary stats 2 emp FLint')
summarize(ZaNew(:));

uratecomp(2,1)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(2,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(2,1) = (PstrucTilde(LYtemp));
LFPratecompB(2,1) = (PstrucTilde(LYtemp+nloc));
migratecomp(2,1) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPWshockAll(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateWshockAll(:,1) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigWshockAll(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #3: 1 s.d. transitory decrease in pi in current loc
%==========================================================================
AdjaNew = getAdjCflPiIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflPiIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
AdjMat(3,:) = AdjaNew;
disp('Adj summary stats 3 emp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 3 emp FLint')
summarize(pimataNew(:));
disp('Z summary stats 3 emp FLint')
summarize(ZaNew(:));

uratecomp(3,1)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(3,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(3,1) = (PstrucTilde(LYtemp));
LFPratecompB(3,1) = (PstrucTilde(LYtemp+nloc));
migratecomp(3,1) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPPishock(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUratePishock(:,1) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigPishock(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #4: 1 s.d. transitory decrease in pi in current loc; correlated across all locations
%==========================================================================
AdjaNew = getAdjCflPiAllIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflPiAllIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
AdjMat(4,:) = AdjaNew;
disp('Adj summary stats 4 emp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 4 emp FLint')
summarize(pimataNew(:));
disp('Z summary stats 4 emp FLint')
summarize(ZaNew(:));

uratecomp(4,1)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(4,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(4,1) = (PstrucTilde(LYtemp));
LFPratecompB(4,1) = (PstrucTilde(LYtemp+nloc));
migratecomp(4,1) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPPishockAll(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUratePishockAll(:,1) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigPishockAll(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #5: one period of zero search costs
%==========================================================================
bstrucTilde = bstruc;
bstrucTilde(nloc+1)=0;
PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);

uratecomp(5,1)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
LFPratecomp(5,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(5,1) = (PstrucTilde(LYtemp));
LFPratecompB(5,1) = (PstrucTilde(LYtemp+nloc));
migratecomp(5,1) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
HeatMapLFPNoSrchCost(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateNoSrchCost(:,1) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigNoSrchCost(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #6: one period of lower moving costs ($20,000 discount)
%==========================================================================
bstrucTilde = bstruc;
bstrucTilde(end-7)=bstrucTilde(end-7).*(1-MCbonus);
PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);

uratecomp(6,1)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
LFPratecomp(6,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(6,1) = (PstrucTilde(LYtemp));
LFPratecompB(6,1) = (PstrucTilde(LYtemp+nloc));
migratecomp(6,1) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
HeatMapLFPMCsubsidy(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateMCsubsidy(:,1) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigMCsubsidy(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #7: one period of lower home production benefits (13%
% reduction)
%==========================================================================
% bstrucTilde = bstruc;
% bstrucTilde(nloc+2)=bstrucTilde(nloc+2).*(1-0.1);
% PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);
% 
% uratecomp(7,1)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
% LFPratecomp(7,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
% LFPratecompA(7,1) = (PstrucTilde(LYtemp));
% LFPratecompB(7,1) = (PstrucTilde(LYtemp+nloc));
% migratecomp(7,1) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
% HeatMapLFPHPsubsidy(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
% HeatMapUrateHPsubsidy(:,1) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
% HeatMapMigHPsubsidy(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #8: one period of lower moving costs conditional on
% moving to Denver, Cincinnati, or San Diego [really good labor markets]
%==========================================================================
% bstrucTilde = bstruc;
% [ZcNewStatic] = getZcflMC(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,distance,MCbonus,susumeLocs,0);
% PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZcNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);
% 
% uratecomp(8,1)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
% LFPratecomp(8,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
% LFPratecompA(8,1) = (PstrucTilde(LYtemp));
% LFPratecompB(8,1) = (PstrucTilde(LYtemp+nloc));
% migratecomp(8,1) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
% inmigratecomp(1,1) = sum(PstrucTilde([susumeLocs susumeLocs+nloc]),2);
% HeatMapLFPTgtMCsubsidy(:,1) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
% HeatMapUrateTgtMCsubsidy(:,1) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
% HeatMapMigTgtMCsubsidy(:,1) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #9: one period of lower moving costs (varying discounts)
%==========================================================================
MCvec  = [-1;-0.93833684;-0.78194736;-0.62555789;-0.46916842;-0.31277895;-0.18766737;-0.15638947;-0.12511158;-0.09383368;-0.06255579;-0.03127789;-0.01563895;-0.01251116;-0.00625558;];
MCamnt = [1.00e6;1.50e05;1.25e05;1.00e05;7.50e04;5.00e04;3.00e04;2.50e04;2.00e04;1.50e04;1.00e04;5.00e03;2.50e03;2.00e03;1.00e03];
urateMCcomp    = nan(length(MCvec),2);
LFPrateMCcomp  = nan(length(MCvec),2);
LFPrateMCcompA = nan(length(MCvec),2);
LFPrateMCcompB = nan(length(MCvec),2);
migrateMCcomp  = nan(length(MCvec),2);
for cc=1:length(MCvec);
    bstrucTilde = bstruc;
    bstrucTilde(end-7)=bstrucTilde(end-7).*(1+MCvec(cc));
    PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);
    
    urateMCcomp(cc,1)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
    LFPrateMCcomp(cc,1) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
    LFPrateMCcompA(cc,1) = (PstrucTilde(LYtemp));
    LFPrateMCcompB(cc,1) = (PstrucTilde(LYtemp+nloc));
    migrateMCcomp(cc,1) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
end

% ========================================================================
% Evaluate Adj and flow utilities at certain characteristics (unemployed)
% ========================================================================
LEi = 0;
experer= wmean(experS(flagS & ageS==ager & empFTS==LEi),PTypel(flagSl & ageS(:)==ager & empFTS(:)==LEi));
tic
AdjNew = getAdjIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,distance,baseAlt,Beta);
disp(['getAdj function took ',num2str(toc/60),' minutes']);

tic
[ZNew      ,pimatNew] = getZintS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,distance,Beta);
disp(['getZ function took ',num2str(toc/60),' minutes']);

%==========================================================================
% Baseline choice probabilities
%==========================================================================
PstrucNew = pclogitAdj1a(bstruc,l,[],ZNew,J,baseAlt,AdjNew,1);

uratebase(1:6,2)   =   kron(ones(6,1),(PstrucNew(LYtemp).*(1-pimatNew(LYtemp))));
LFPratebase(1:6,2) =   kron(ones(6,1),(PstrucNew(LYtemp))./(sum(PstrucNew([LYtemp LYtemp+nloc]),2)));
LFPratebaseA(1:6,2)=   kron(ones(6,1),(PstrucNew(LYtemp)));
LFPratebaseB(1:6,2)=   kron(ones(6,1),(PstrucNew(LYtemp+nloc)));
migratebase(1:6,2) = 1-kron(ones(6,1),sum(PstrucNew([LYtemp LYtemp+nloc]),2));
inmigratebase(1,2) = sum(PstrucNew([susumeLocs susumeLocs+nloc]),2);

%==========================================================================
% Counterfactual Sim #1: 1 s.d. transitory decrease in wage in current loc
%==========================================================================
AdjaNew = getAdjCflWageIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflWageIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
disp('Adj summary stats 1 unemp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 1 unemp FLint')
summarize(pimataNew(:));
disp('Z summary stats 1 unemp FLint')
summarize(ZaNew(:));

uratecomp(1,2)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(1,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(1,2) = (PstrucTilde(LYtemp));
LFPratecompB(1,2) = (PstrucTilde(LYtemp+nloc));
migratecomp(1,2) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPWshock(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateWshock(:,2) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigWshock(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #2: 1 s.d. transitory decrease in wage in current loc; correlated across all locations
%==========================================================================
AdjaNew = getAdjCflWageAllIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflWageAllIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,wage_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
disp('Adj summary stats 2 unemp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 2 unemp FLint')
summarize(pimataNew(:));
disp('Z summary stats 2 unemp FLint')
summarize(ZaNew(:));

uratecomp(2,2)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(2,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(2,2) = (PstrucTilde(LYtemp));
LFPratecompB(2,2) = (PstrucTilde(LYtemp+nloc));
migratecomp(2,2) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPWshockAll(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateWshockAll(:,2) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigWshockAll(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #3: 1 s.d. transitory decrease in pi in current loc
%==========================================================================
AdjaNew = getAdjCflPiIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflPiIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
disp('Adj summary stats 3 unemp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 3 unemp FLint')
summarize(pimataNew(:));
disp('Z summary stats 3 unemp FLint')
summarize(ZaNew(:));

uratecomp(3,2)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(3,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(3,2) = (PstrucTilde(LYtemp));
LFPratecompB(3,2) = (PstrucTilde(LYtemp+nloc));
migratecomp(3,2) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPPishock(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUratePishock(:,2) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigPishock(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #4: 1 s.d. transitory decrease in pi in current loc; correlated across all locations
%==========================================================================
AdjaNew = getAdjCflPiAllIntS(bHet,draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,baseAlt,Beta);
[ZaNew      ,pimataNew] = getZcflPiAllIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,pi_shock_amt,distance,Beta);
PstrucTilde = pclogitAdj1a(bstruc,l,[],ZaNew,J,baseAlt,AdjaNew,1);
disp('Adj summary stats 4 unemp FLint')
summarize(AdjaNew(:));
disp('pimat summary stats 4 unemp FLint')
summarize(pimataNew(:));
disp('Z summary stats 4 unemp FLint')
summarize(ZaNew(:));

uratecomp(4,2)   = sum(PstrucTilde(LYtemp).*(1-pimataNew(LYtemp)));
LFPratecomp(4,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(4,2) = (PstrucTilde(LYtemp));
LFPratecompB(4,2) = (PstrucTilde(LYtemp+nloc));
migratecomp(4,2) = 1-sum(PstrucTilde([LYtemp LYtemp+nloc]),2);
HeatMapLFPPishockAll(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUratePishockAll(:,2) = PstrucTilde(1:nloc).*(1-pimataNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigPishockAll(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #5: one period of zero search costs
%==========================================================================
bstrucTilde = bstruc;
bstrucTilde(nloc+1)=0;
PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);

uratecomp(5,2)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
LFPratecomp(5,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(5,2) = (PstrucTilde(LYtemp));
LFPratecompB(5,2) = (PstrucTilde(LYtemp+nloc));
migratecomp(5,2) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
HeatMapLFPNoSrchCost(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateNoSrchCost(:,2) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigNoSrchCost(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #6: one period of lower moving costs ($20,000 discount)
%==========================================================================
bstrucTilde = bstruc;
bstrucTilde(end-7)=bstrucTilde(end-7).*(1-MCbonus);
PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);

uratecomp(6,2)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
LFPratecomp(6,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
LFPratecompA(6,2) = (PstrucTilde(LYtemp));
LFPratecompB(6,2) = (PstrucTilde(LYtemp+nloc));
migratecomp(6,2) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
HeatMapLFPMCsubsidy(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
HeatMapUrateMCsubsidy(:,2) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
HeatMapMigMCsubsidy(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #7: one period of lower home production benefits (13%
% reduction)
%==========================================================================
% bstrucTilde = bstruc;
% bstrucTilde(nloc+2)=bstrucTilde(nloc+2).*(1-0.13);
% PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);
% 
% uratecomp(7,2)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
% LFPratecomp(7,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
% LFPratecompA(7,2) = (PstrucTilde(LYtemp));
% LFPratecompB(7,2) = (PstrucTilde(LYtemp+nloc));
% migratecomp(7,2) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
% HeatMapLFPHPsubsidy(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
% HeatMapUrateHPsubsidy(:,2) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
% HeatMapMigHPsubsidy(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #8: one period of lower moving costs conditional on
% moving to Denver, Cincinnati, or San Diego [really good labor markets]
%==========================================================================
% bstrucTilde = bstruc;
% [ZcNewStatic] = getZcflMC(draws,ARcov,ARcov2,yr,ager,LEi,LYtemp,experer,typer,birthLoci,birthDivi,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,distance,MCbonus,susumeLocs,0);
% PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZcNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);
% 
% uratecomp(8,2)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
% LFPratecomp(8,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
% LFPratecompA(8,2) = (PstrucTilde(LYtemp));
% LFPratecompB(8,2) = (PstrucTilde(LYtemp+nloc));
% migratecomp(8,2) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
% inmigratecomp(1,2) = sum(PstrucTilde([susumeLocs susumeLocs+nloc]),2);
% HeatMapLFPTgtMCsubsidy(:,2) = PstrucTilde(1:nloc)'./sum(reshape(PstrucTilde,nloc,2),2) - PstrucNew(1:nloc)'./sum(reshape(PstrucNew,nloc,2),2);
% HeatMapUrateTgtMCsubsidy(:,2) = PstrucTilde(1:nloc).*(1-pimatNew)-PstrucNew(1:nloc).*(1-pimatNew);
% HeatMapMigTgtMCsubsidy(:,2) = sum(reshape(PstrucTilde,nloc,2),2)-sum(reshape(PstrucNew,nloc,2),2);

%==========================================================================
% Counterfactual Sim #9: one period of lower moving costs (varying discounts)
%==========================================================================
for cc=1:length(MCvec);
    bstrucTilde = bstruc;
    bstrucTilde(end-7)=bstrucTilde(end-7).*(1+MCvec(cc));
    PstrucTilde = pclogitAdj2a(bstrucTilde,bstruc,l,[],[],ZNewStatic,ZNew-ZNewStatic,J,baseAlt,AdjNew,1);
    
    urateMCcomp(cc,2)   = sum(PstrucTilde(LYtemp).*(1-pimatNew(LYtemp)));
    LFPrateMCcomp(cc,2) = (PstrucTilde(LYtemp))./(sum(PstrucTilde([LYtemp LYtemp+nloc]),2));
    LFPrateMCcompA(cc,2) = (PstrucTilde(LYtemp));
    LFPrateMCcompB(cc,2) = (PstrucTilde(LYtemp+nloc));
    migrateMCcomp(cc,2) = sum(PstrucTilde(setdiff(1:J,[LYtemp LYtemp+nloc])),2);
end

% uratecomp(8,2)   = 0;
% LFPratecomp(8,2) = 0;
% LFPratecompA(8,2) = 0;
% LFPratecompB(8,2) = 0;
% migratecomp(8,2) = 0;
% inmigratecomp(1,2) = 0;

[uratecomp     uratebase    ]
[LFPratecomp   LFPratebase  ]
[LFPratecompA  LFPratebaseA ]
[LFPratecompB  LFPratebaseB ]
[migratecomp   migratebase  ]
% [inmigratecomp inmigratebase]

[uratecomp-uratebase   ]
[LFPratecomp-LFPratebase ]
[LFPratecompA-LFPratebaseA]
[LFPratecompB-LFPratebaseB]
[migratecomp-migratebase ]
% inmigratecomp-inmigratebase

dlmwrite(['UrateCfl',   num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[uratecomp-uratebase]);
dlmwrite(['LFPrateCfl', num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[LFPratecomp-LFPratebase]);
dlmwrite(['LFPrateACfl',num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[LFPratecompA-LFPratebaseA]);
dlmwrite(['LFPrateBCfl',num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[LFPratecompB-LFPratebaseB]);
dlmwrite(['MigRateCfl', num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[migratecomp-migratebase]);

dlmwrite(['UrateCfl',   num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[uratebase(1,:)],'-append');
dlmwrite(['LFPrateCfl', num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[LFPratebase(1,:)],'-append');
dlmwrite(['LFPrateACfl',num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[LFPratebaseA(1,:)],'-append');
dlmwrite(['LFPrateBCfl',num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[LFPratebaseB(1,:)],'-append');
dlmwrite(['MigRateCfl', num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[migratebase(1,:)],'-append');

dlmwrite(['MigRateCflMCvec',num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer),'.csv'],[MCamnt migrateMCcomp-kron(ones(length(MCvec),1),migratebase(3,:))]);

CovarE = [ones(nloc,1) (bstruc(1:nloc)-bstruc(LYtemp))./std(bstruc(1:nloc)) (wHat(:,yr-2003)-wHat(LYtemp,yr-2003))./std(wHat(:,yr-2003)) (lambda_e_hat(:,yr-2003)-lambda_e_hat(LYtemp,yr-2003))./std(lambda_e_hat(:,yr-2003)) bARw(1:nloc)./std(bARw(1:nloc)) bARw(nloc+2:end)./std(bARw(nloc+2:end)) rho_hat_urate(1,:)'./std(rho_hat_urate(1,:)) rho_hat_urate(2,:)'./std(rho_hat_urate(2,:)) sig_hat_urate(1,:)'./std(sig_hat_urate(1,:)) log(dist110(1:nloc,LYtemp)) birthLoci(1:nloc)' birthDivi(1:nloc)'];
CovarU = [ones(nloc,1) (bstruc(1:nloc)-bstruc(LYtemp))./std(bstruc(1:nloc)) (wHat(:,yr-2003)-wHat(LYtemp,yr-2003))./std(wHat(:,yr-2003)) (lambda_u_hat(:,yr-2003)-lambda_u_hat(LYtemp,yr-2003))./std(lambda_u_hat(:,yr-2003)) bARw(1:nloc)./std(bARw(1:nloc)) bARw(nloc+2:end)./std(bARw(nloc+2:end)) rho_hat_urate(1,:)'./std(rho_hat_urate(1,:)) rho_hat_urate(2,:)'./std(rho_hat_urate(2,:)) sig_hat_urate(1,:)'./std(sig_hat_urate(1,:)) log(dist110(1:nloc,LYtemp)) birthLoci(1:nloc)' birthDivi(1:nloc)'];

save(['Results',num2str(ager),locstring,num2str(yr),'Born',bstatestring,'t',num2str(typer)],'AdjMat','migratecomp','migratebase','LFPratecomp*','LFPratebase*','uratebase','uratecomp','HeatMap*','Covar*','bstruc','wHat','lambda_e_hat','lambda_u_hat','dist110','birthLoci','birthDivi','nloc','LYtemp','yr');

disp('Urate Heat Map Regression of Cfl1 (wage shock)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateWshock(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateWshock(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Urate Heat Map Regression of Cfl1 (wage shock everywhere)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateWshockAll(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateWshockAll(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Urate Heat Map Regression of Cfl2 (pi shock)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUratePishock(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUratePishock(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Urate Heat Map Regression of Cfl2 (pi shock everywhere)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUratePishockAll(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUratePishockAll(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Urate Heat Map Regression of Cfl3 (0 search cost)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateNoSrchCost(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateNoSrchCost(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Urate Heat Map Regression of Cfl4 (MC subsidy)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateMCsubsidy(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateMCsubsidy(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
% disp('Urate Heat Map Regression of Cfl5 (HP benefit decrease)');
% [bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateHPsubsidy(setdiff(1:nloc,LYtemp),1));
% [bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateHPsubsidy(setdiff(1:nloc,LYtemp),2));
% [bHMe seHMe bHMe./seHMe]
% [bHMu seHMu bHMu./seHMu]
% disp('Urate Heat Map Regression of Cfl6 (Targeted move subsidy)');
% [bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateTgtMCsubsidy(setdiff(1:nloc,LYtemp),1));
% [bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapUrateTgtMCsubsidy(setdiff(1:nloc,LYtemp),2));
% [bHMe seHMe bHMe./seHMe]
% [bHMu seHMu bHMu./seHMu]


disp('LFP Heat Map Regression of Cfl1 (wage shock)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPWshock(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPWshock(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('LFP Heat Map Regression of Cfl1 (wage shock everywhere)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPWshockAll(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPWshockAll(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('LFP Heat Map Regression of Cfl2 (pi shock)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPPishock(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPPishock(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('LFP Heat Map Regression of Cfl2 (pi shock everywhere)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPPishockAll(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPPishockAll(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('LFP Heat Map Regression of Cfl3 (0 search cost)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPNoSrchCost(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPNoSrchCost(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('LFP Heat Map Regression of Cfl4 (MC subsidy)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPMCsubsidy(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPMCsubsidy(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
% disp('LFP Heat Map Regression of Cfl5 (HP benefit decrease)');
% [bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPHPsubsidy(setdiff(1:nloc,LYtemp),1));
% [bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPHPsubsidy(setdiff(1:nloc,LYtemp),2));
% [bHMe seHMe bHMe./seHMe]
% [bHMu seHMu bHMu./seHMu]
% disp('LFP Heat Map Regression of Cfl6 (Targeted move subsidy)');
% [bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPTgtMCsubsidy(setdiff(1:nloc,LYtemp),1));
% [bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapLFPTgtMCsubsidy(setdiff(1:nloc,LYtemp),2));
% [bHMe seHMe bHMe./seHMe]
% [bHMu seHMu bHMu./seHMu]


disp('Migration Heat Map Regression of Cfl1 (wage shock)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigWshock(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigWshock(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Migration Heat Map Regression of Cfl1 (wage shock everywhere)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigWshockAll(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigWshockAll(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Migration Heat Map Regression of Cfl2 (pi shock)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigPishock(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigPishock(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Migration Heat Map Regression of Cfl2 (pi shock everywhere)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigPishockAll(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigPishockAll(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Migration Heat Map Regression of Cfl3 (0 search cost)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigNoSrchCost(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigNoSrchCost(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
disp('Migration Heat Map Regression of Cfl4 (MC subsidy)');
[bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),1));
[bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigMCsubsidy(setdiff(1:nloc,LYtemp),2));
[bHMe seHMe bHMe./seHMe]
[bHMu seHMu bHMu./seHMu]
% disp('Migration Heat Map Regression of Cfl5 (HP benefit decrease)');
% [bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigHPsubsidy(setdiff(1:nloc,LYtemp),1));
% [bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigHPsubsidy(setdiff(1:nloc,LYtemp),2));
% [bHMe seHMe bHMe./seHMe]
% [bHMu seHMu bHMu./seHMu]
% disp('Migration Heat Map Regression of Cfl6 (Targeted move subsidy)');
% [bHMe,seHMe] = lscov(CovarE(setdiff(1:nloc,LYtemp),:),100*HeatMapMigTgtMCsubsidy(setdiff(1:nloc,LYtemp),1));
% [bHMu,seHMu] = lscov(CovarU(setdiff(1:nloc,LYtemp),:),100*HeatMapMigTgtMCsubsidy(setdiff(1:nloc,LYtemp),2));
% [bHMe seHMe bHMe./seHMe]
% [bHMu seHMu bHMu./seHMu]
