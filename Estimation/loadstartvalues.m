%==========================================================================
% Estimate friction parameters
%==========================================================================
loc_dummies = [];
for j=setdiff(1:nloc,baseLoc);
	loc_dummies = cat(2,loc_dummies,Choice(:)==j);
end

% logit of empFT conditional on inlf and empFT_lag
conditioner = inlf(:)==1 & empFT_lag(:)==1 & flag==1;
conditionerS = inlfS(:)==1 & empFT_lagS(:)==1 & flagS(:)==1;
diffloc = Choicelag(:)~=Choice(:) & Choicelag(:)~=Choice(:)-nloc & Choicelag(:)-55~=Choice(:);
if nloc<=53
	Xtild = [loc_dummies UrateLag(:) diffloc exper(:) exper(:).^2./100];
elseif nloc==55 && strcmp(sample,'HS')==1
	Xtild = [loc_dummies UrateLag(:) diffloc exper(:) exper(:).^2./100];
elseif nloc==55 && strcmp(sample,'HS')==0
	% Lump sparse locations together for college grads
	Xtild = [cat(2,loc_dummies(:,1:end-2),Choice(:)==54 | Choice(:)==55) UrateLag(:) diffloc exper(:) exper(:).^2./100];
end
% Xtild = [loc_dummies time_dummies loc_time_dummies diffloc];
Ytild = empFT(:);
bLemp = glmfit(Xtild(conditioner,:),Ytild(conditioner,:), 'binomial', 'link', 'logit');
bLempTemp = cat(1,bLemp,.2*rand(S-1,1));
dat.Pemp = glmval(bLempTemp,cat(2,kron(ones(S,1),Xtild),kron(eye(S,S-1),ones(size(Xtild,1),1))),'logit');
dat.Pemp(~conditionerS) = 0;
dat.Pemp = reshape(dat.Pemp,[N T S]);
bLempMat = bLempTemp;
if nloc==55 && strcmp(sample,'HS')==0
   bLemp = [bLemp(1:end-4);bLemp(end-4);bLemp(end-3:end)]; 
end
bLempSEmat = zeros(size(bLempMat));

% logit of empFT conditional on inlf and ~empFT_lag
conditioner1 = inlf(:)==1 & empFT_lag(:)==0 & flag==1;
conditioner1S = inlfS(:)==1 & empFT_lagS(:)==0 & flagS(:)==1;
if nloc<=53
	Xtild = [loc_dummies UrateLag(:) diffloc exper(:) exper(:).^2./100];
elseif nloc==55 && strcmp(sample,'HS')==1
	Xtild = [cat(2,loc_dummies(:,1:end-2),Choice(:)==54 | Choice(:)==55) UrateLag(:) diffloc exper(:) exper(:).^2./100];
elseif nloc==55 && strcmp(sample,'HS')==0
	% Lump sparse locations together for college grads
	Xtild = [cat(2,loc_dummies(:,[1:end-4 end-2]),(Choice(:)==52 | Choice(:)==54 | Choice(:)==55)) UrateLag(:) diffloc exper(:) exper(:).^2./100];
end
Ytild = empFT(:);
bLunemp = glmfit(Xtild(conditioner1,:),Ytild(conditioner1,:), 'binomial', 'link', 'logit');
bLunempTemp = cat(1,bLunemp,.2*rand(S-1,1));
dat.Punemp = glmval(bLunempTemp,cat(2,kron(ones(S,1),Xtild),kron(eye(S,S-1),ones(size(Xtild,1),1))),'logit');
dat.Punemp(~conditioner1S) = 0;
dat.Punemp = reshape(dat.Punemp,[N T S]);
if nloc==55 && strcmp(sample,'HS')==1
	bLunemp = [bLunemp(1:end-4);bLunemp(end-4);bLunemp(end-3:end)]; 
elseif nloc==55 && strcmp(sample,'HS')==0
	bLunemp = [bLunemp(1:end-6);bLunemp(end-4);bLunemp(end-5);bLunemp(end-4);bLunemp(end-4);bLunemp(end-3:end)]; 
end
bLunempMat = cat(1,bLunemp,bLunempTemp(end));
bLunempSEmat = zeros(size(bLunempMat));

whos bLemp bLunemp

if iteration==1;
	% Non-stationary case
	[bLemp bLunemp]
	[sum(conditioner) sum(conditioner1)]
	delta_hat = nan(nloc,size(urate55,2));
	lambda_hat = nan(nloc,size(urate55,2));
	lambda_e_hat = nan(nloc,size(urate55,2));
	lambda_u_hat = nan(nloc,size(urate55,2));
	for j=1:nloc;
		if j<baseLoc
			locer = [zeros(j,1); 1; zeros(nloc-j-1,1);];
		elseif j==baseLoc
			locer = zeros(nloc-1,1);
		elseif j>=baseLoc
			locer = [zeros(j-2,1); 1; zeros(nloc-j,1)];
		end
		for t=1:size(urate55,2)
			indicter = [1; locer; 100*urate_lag55(j,t); 0; 0; 0]';
			delta_hat(j,t)   = 1./(1+exp(indicter*bLemp));
			lambda_hat(j,t)  = exp(indicter*bLunemp)./(1+exp(indicter*bLunemp));
			indicter = [1;locer; 100*urate_lag55(j,t); 1; 0; 0]';
			lambda_e_hat(j,t)= exp(indicter*bLemp)./(1+exp(indicter*bLemp));
			lambda_u_hat(j,t)= exp(indicter*bLunemp)./(1+exp(indicter*bLunemp));
		end
	end

	disp('Delta hat in 2004, 2009, 2013')
	[delta_hat(:,1) delta_hat(:,5) delta_hat(:,10)]
	disp('Lambda e hat in 2004, 2009, 2013')
	[lambda_e_hat(:,1) lambda_e_hat(:,5) lambda_e_hat(:,10)]
	disp('Lambda hat in 2004, 2009, 2013')
	[lambda_hat(:,1) lambda_hat(:,5) lambda_hat(:,10)]
	disp('Lambda u hat in 2004, 2009, 2013')
	[lambda_u_hat(:,1) lambda_u_hat(:,5) lambda_u_hat(:,10)]
end

%==========================================================================
% Estimate unemployment autocorrelation
%==========================================================================
rho_hat_urate    = zeros(nloc,2)';
sig_hat_urate    = zeros(nloc,1)';
resid_urate      = zeros(nloc,9)';

for j=1:nloc
   [rho_hat_urate(:,j)   ,~,resid_urate(:,j)   ,~,stemp] = regress(urate_lag55(j,2:10)'   ,[ones(9,1) urate_lag55(j,1:9)'   ]);
   sig_hat_urate(j) = sqrt(stemp(end));
end
disp('Summary statistics on urate autocorrelations');
summarize([rho_hat_urate([2 1],:)' sig_hat_urate']);

pimatO     = nan(N*T,nloc);
pimat      = nan(N*T,nloc);
pimatVec   = nan(N*T,nloc);
pimat1Vec  = nan(N*T,1);
pimat1Int  = nan(N*T,nloc);
pimat1O    = nan(N*T,nloc);
pimat1     = nan(N*T,nloc);
omegaO     = nan(N*T,nloc);
omega      = nan(N*T,nloc);
omegaVec   = nan(N*T,nloc);
omegaInt  = nan(N*T,nloc);
pimat1w   = nan(N*T,nloc,draws);
omegaw    = nan(N*T,nloc,draws);

Tmat    = zeros(N*T,size(urate55,2));
Tvec    = calyr(:)-2003;
for i=1:N*T
	if flag(i)==1
	    Tmat(i,Tvec(i)) = 1;
	else
	    Tmat(i,1)       = 1;
	end
end

tic
for k=1:nloc
	ftemp = flag & empFT_lag(:)==0 & (k==LY(:)-nloc*(LY(:)>nloc)); 
	pimat(ftemp,k) =   lambdervecRho(nloc,k,Tmat(ftemp,:),Tvec(ftemp),urate_lag55,empFT_lag(ftemp)==1,1,exper(ftemp),rho_hat_urate  ,0,bLemp,bLunemp);

	ftemp = flag & empFT_lag(:)==1 & (k==LY(:)-nloc*(LY(:)>nloc));
	pimat(ftemp,k) = 1-lambdervecRho(nloc,k,Tmat(ftemp,:),Tvec(ftemp),urate_lag55,empFT_lag(ftemp)==1,1,exper(ftemp),rho_hat_urate   ,0,bLemp,bLunemp);

	ftemp = flag & empFT_lag(:)==0 & (k~=LY(:)-nloc*(LY(:)>nloc));
	pimat(ftemp,k) =   lambdervecRho(nloc,k,Tmat(ftemp,:),Tvec(ftemp),urate_lag55,empFT_lag(ftemp)==1,0,exper(ftemp),rho_hat_urate,0,bLemp,bLunemp);

	ftemp = flag & empFT_lag(:)==1 & (k~=LY(:)-nloc*(LY(:)>nloc));
	pimat(ftemp,k) =   lambdervecRho(nloc,k,Tmat(ftemp,:),Tvec(ftemp),urate_lag55,empFT_lag(ftemp)==1,0,exper(ftemp),rho_hat_urate,0,bLemp,bLunemp);
	
	ftemp = flag & (k==LY(:)-nloc*(LY(:)>nloc));
	pimat1(ftemp,1) =   lambdervecRho(nloc,k,Tmat(ftemp,:),Tvec(ftemp),urate_lag55,0,1,exper(ftemp),rho_hat_urate,1,bLemp,bLunemp);
end

for k=1:nloc
	omega(:,k) = pimat(:,k)./pimat1(:,1);
end
if any(any(any(isnan(pimat(flag,:)))))
	error('pi''s not set up correctly')
end
if any(any(any(isnan(pimat1(flag,1)))))
	error('pi_t+1''s not set up correctly')
end
disp(['Pimat construction took ',num2str(toc/60),' minutes']); %takes 1.5 minutes without quadrature, 11 minutes with quadrature

%==========================================================================
% Set up wage regression matrices
%==========================================================================
loc_dummies = zeros(N*T,nloc-1);
for j=setdiff(1:nloc,baseLoc);
	loc_dummies(:,j-1) = Choice(:)==j;
end

time_dummies = zeros(N*T,9);
for k=2005:2013
   time_dummies(:,k-2004) = calyr(:)==k;
end

loc_time_dummies = zeros(N*T,(nloc-1)*9);
for j=setdiff(1:nloc,baseLoc);
	for k=2005:2013
	    loc_time_dummies(:,(j-2)*9+k-2004) = Choice(:)==j & calyr(:)==k;
	end
end

if strcmp(money,'wage')==1
	wagevar = lnWageHr;
elseif strcmp(money,'earn')==1
	wagevar = lnearnfinalJb;
end

%==========================================================================
% Estimate wage parameters
%==========================================================================
% construct flag for wages that will actually be used
wageflagtemp = wageflag==1 & Choice<=nloc & flagw & empFT==1; % Note the addition that employmentFT must hold!!
size(wageflagtemp)
% [bwage  ,sebhatw  ] = lscov([ones(sum(wageflagtemp(:)),1) loc_dummies(wageflagtemp,:) exper(wageflagtemp) exper(wageflagtemp).^2./100],wagevar(wageflagtemp));
% [bwage   sebhatw  ]

if nloc==55
	[bwager  ,sebhatwr  , mser] = lscov([ones(sum(wageflagtemp(:)),1) loc_dummies(wageflagtemp,:) time_dummies(wageflagtemp,:) loc_time_dummies(wageflagtemp,1:end-18) exper(wageflagtemp) exper(wageflagtemp).^2./100],wagevar(wageflagtemp));
	[bwager   sebhatwr  ];
else
	[bwager  ,sebhatwr  ] = lscov([ones(sum(wageflagtemp(:)),1) loc_dummies(wageflagtemp,:) time_dummies(wageflagtemp,:) loc_time_dummies(wageflagtemp,:) exper(wageflagtemp) exper(wageflagtemp).^2./100],wagevar(wageflagtemp));
	[bwager   sebhatwr  ];
end

optionswage = optimset('Disp','iter','LargeScale','on','maxiter',1e8,'maxfuneval',1e8,'TolX',1e-6,'Tolfun',1e-6,'DerivativeCheck','off','GradObj','on','FinDiffType','central');
[parmwage,~,~,~,~,hwage] = fminunc('normalMLEw',[bwager;sqrt(mser)],optionswage,[],wagevar(wageflagtemp),[ones(sum(wageflagtemp(:)),1) loc_dummies(wageflagtemp,:) time_dummies(wageflagtemp,:) loc_time_dummies(wageflagtemp,1:end-18) exper(wageflagtemp) exper(wageflagtemp).^2./100],[]);

bwage = parmwage(1:end-1);
sebhatw = sqrt(diag(eye(size(hwage))/hwage));
sebhatw = sebhatw(1:end-1);
parms.wageBeta = cat(1,bwage,.2*rand(S-1,1));
parms.wageSig  = parmwage(end);
bwageMat = parms.wageBeta;

NTwage = sum(wageflagtemp(:))
Nwage  = numel(unique(ID(wageflagtemp)))
[bwage sebhatw];
wHat      = zeros(nloc,10);
wHatOther = zeros(nloc,10);
for j=1:nloc
	for t=1:10
	    wHat(j,t) = ewage(nloc,j,t,0,bwage)-bwage(1);
	end
end

%==========================================================================
% Estimate wage autocorrelations
%==========================================================================
rho_hat_wageOld      = zeros(nloc,2)';
resid_wageOld        = zeros(nloc,9)';
sig_hat_wageOld      = zeros(nloc,1)';
for j=1:nloc
   [rho_hat_wageOld(:,j),~,resid_wageOld(:,j),~,stemp] = regress(wHat(j,2:10)',[ones(9,1) wHat(j,1:9)']);
   sig_hat_wageOld(j)   = sqrt(stemp(end));
end
disp('Summary statistics on wage autocorrelations');
summarize([rho_hat_wageOld([2 1],:)' sig_hat_wageOld']);

tm=sum(abs(rho_hat_wageOld(2,:)>1));
disp(['Number of locations with autocorrelation outside the unit circle: ',num2str(tm)]);

%==========================================================================
% Estimate wage autocorrelations --- pooled across locations, but
% heteroskedastic errors and heterogeneous drift
%==========================================================================
rho_hat_wage    = zeros(nloc,2)';
resid_wage      = zeros(nloc,9)';
sig_hat_wage    = zeros(nloc,1)';
YARw = reshape(wHat(:,2:10)',nloc*9,1);
XARw = zeros(nloc*9,nloc+1);
for jj=1:nloc
   XARw((jj-1)*9+1:jj*9,jj)=ones(9,1);
end
XARw(:,end) = reshape(wHat(:,1:9)',nloc*9,1);
dw   = reshape(([1:nloc]'*ones(1,9))',nloc*9,1);
if exist(['bARw',num2str(nloc),'loc',sample,'.mat'],'file')==2
	load(['bARw',num2str(nloc),'loc',sample,'.mat'],'bARw');
	svalARw = bARw;
else
	svalARw = [rho_hat_wageOld(1,:)';.70;sig_hat_wageOld'];
end
optionswage = optimset('Disp','iter','LargeScale','on','maxiter',1e8,'maxfuneval',1e8,'TolX',1e-6,'Tolfun',1e-6,'DerivativeCheck','off','GradObj','on','FinDiffType','central');
[bARw,~,~,~,~,hARw] = fminunc('normalMLE',svalARw,optionswage,[],YARw,XARw,dw);
save(['bARw',num2str(nloc),'loc',sample,'.mat'],'bARw');
[bARw sqrt(diag(inv(full(hARw))))]
residCoolVec = YARw-XARw*bARw(1:nloc+1);
for jj=1:nloc
	resid_wage(1:9,jj) = residCoolVec((jj-1)*9+1:jj*9,1);
end
rho_hat_wage(1,:)   = bARw(1:nloc)';
rho_hat_wage(2,:)   = bARw(nloc+1)*ones(1,nloc);
sig_hat_wage(1:end) = bARw(end-(nloc-1):end)';

disp('Summary statistics on cool wage autocorrelations');
summarize([rho_hat_wage([2 1],:)' sig_hat_wage']);

%==========================================================================
% Estimate wage autocorrelations --- pooled across locations
%==========================================================================
rho_hat_wageUncool      = zeros(1,2)';
resid_wageUnCool        = zeros(nloc,9)';
sig_hat_wageUncool      = zeros(1,1)';
[rho_hat_wageUncool,~,resid_wageUncool,~,stempUncool] = regress(reshape(wHat(:,2:10)',nloc*9,1),[ones(nloc*9,1) reshape(wHat(:,1:9)',nloc*9,1)]);
% [rho_hat_wage2,~,~,~,stemp2] = regress(reshape(wHat(:,3:10)',nloc*8,1),[ones(nloc*8,1) reshape(wHat(:,2:9)',nloc*8,1) reshape(wHat(:,1:8)',nloc*8,1)]);
[temp,setemp] = lscov([ones(nloc*9,1) reshape(wHat(:,1:9)',nloc*9,1)],reshape(wHat(:,2:10)',nloc*9,1));
% [temp2,setemp2] = lscov([ones(nloc*8,1) reshape(wHat(:,2:9)',nloc*8,1) reshape(wHat(:,1:8)',nloc*8,1)],reshape(wHat(:,3:10)',nloc*8,1));
% [temp3,setemp3] = lscov([ones(nloc*7,1) reshape(wHat(:,3:9)',nloc*7,1) reshape(wHat(:,2:8)',nloc*7,1) reshape(wHat(:,1:7)',nloc*7,1)],reshape(wHat(:,4:10)',nloc*7,1));
sig_hat_wageUncool      = sqrt(stemp(end));

disp('Summary statistics on uncool wage autocorrelations');
summarize([rho_hat_wageUncool([2 1],:)' sig_hat_wageUncool']);

%==========================================================================
% Correlations between wage and employment shocks
%==========================================================================
[cor,Pcor]=corr([resid_wageOld(:) resid_urate(:)])
[cor,Pcor]=corr([resid_wage(:) resid_urate(:)])
[cor,Pcor]=corr([resid_wageUncool(:) resid_urate(:)])

tic;
shocker = nan(N*T,2*nloc,draws);
empshocker = nan(size(resid_wage,1),nloc);
mastershocker = zeros(size(resid_wage,1),nloc*2);
rhoshocker = zeros(nloc*2,1);
for jj=1:nloc
	empshocker(:,jj) = resid_urate(:,jj);
	mastershocker(:,(jj-1)*2+1:jj*2) = [resid_wage(:,jj) resid_urate(:,jj)];
	rhoshocker((jj-1)*2+1:jj*2,1) = [bARw(nloc+1) rho_hat_urate(2,j)];
end

[ARcor,pARcor]=corr(mastershocker)
ARcov=cov(mastershocker)
ARcov2=ARcov+ARcov.*(rhoshocker*rhoshocker');
disp(['Time spent drawing shocks: ',num2str(toc/60),' minutes']);

%==========================================================================
% Correlation matrices of interest
%==========================================================================
for t=1:10
%     disp(['Year ',num2str(2003+t),':']);
	[CorrLLM(:,:,t)    ,PLLM(:,:,t)    ]=corr([wHat(:,t) delta_hat(:,t) lambda_hat(:,t)]);
	[RankCorrLLM(:,:,t),RankPLLM(:,:,t)]=corr([wHat(:,t) delta_hat(:,t) lambda_hat(:,t)],'type','Spearman');
end
save LLCresults wHat delta_hat lambda_hat

%========================================================================
% Helpful matrices
%==========================================================================
calyr(isnan(calyr))=2004;
LYi     = ones(N*T,1);
homeLoc = ones(N*T,1);
Lh      = zeros(N*T,J);
Lw      = zeros(N*T,J);
Tmat    = zeros(N*T,size(urate55,2));
Tvec    = calyr(:)-2003;
tic;
for i=1:N*T
	if flag(i)==1
	    LYi             = LY(i);
	    l               = LYi-nloc*(LYi>nloc);
	    homeLoc(i)      = l;
	    Lh(i,nloc+l)    = 1;
	    Lw(i,     l)    = 1;
	    Tmat(i,Tvec(i)) = 1;
	else
	    Tmat(i,1)       = 1;
	end
end
LhS = Lh;
LwS = Lw;
Lh = repmat(Lh,[1 1 draws]);
Lw = repmat(Lw,[1 1 draws]);
homeLocW= repmat(homeLoc,[1 1 draws]);
pvEmpL   = inlf_lag(:)    ==1 & empFT_lag(:)    ==1;
pvUnempL = inlf_lag(:)    ==1 & empFT_lag(:)    ==0;
pvEmp    = inlf_lag(flagw)==1 & empFT_lag(flagw)==1;
pvUnemp  = inlf_lag(flagw)==1 & empFT_lag(flagw)==0;
disp(['Helpful matrices loop took ',num2str(toc),' seconds']);

%========================================================================
% Flexible logit starting values
%==========================================================================
Adj = zeros(N*T,J);
if Beta>0
	%--------------------------------------------------------------------------
	% prepare data
	%--------------------------------------------------------------------------
	X = [ones(N*T,2)];
	Z = zeros(N*T,15,J);
	tic
	for j=1:nloc
	    k=j+nloc;
	    Z(:,:,j) = [(pvEmpL==1).*pimat(:,j) (pvUnempL==1).*pimat(:,j) ewagevecRho(nloc,j,Tmat,Tvec,rho_hat_wage,0,exper(:),bwage)  birthLoc(:,j) birthDiv(:,j) switchervec(nloc,LY(:),j,age(:)) movervec(nloc,LY(:),j,distance,age(:),pvEmpL,pvUnempL)];
	    Z(:,:,k) = [zeros(N*T,3)                                                                                                   birthLoc(:,k) birthDiv(:,k) switchervec(nloc,LY(:),k,age(:)) movervec(nloc,LY(:),k,distance,age(:),pvEmpL,pvUnempL)];
	end
	disp(['Flexible logit data matrices loop took ',num2str(toc/60),' minutes']);
	Y = Choice(:);
	summarize(Z(flag,:,1));

	%--------------------------------------------------------------------------
	% Load no-het estimates
	%--------------------------------------------------------------------------
	load([money,'_',num2str(nloc),'loc',sample,'.mat'],'b');
	% Predict probabilities
	P = pclogit1(b,Y(flag),X(flag,:),Z(flag,:,:),size(X,2),size(Z,2),J,baseAlt);
	summarize(P);
	P([1 10 11 16 17],:)
	NflexLogit       = sum(flag)
	NflexLogitUnique = numel(unique(ID(flagw)))
	% load flexLogitStartVal b
end





%==========================================================================
% New versions of these matrices (to generate type-specific probabilities)
%==========================================================================
pimatSO     = nan(N*T*S,nloc);
pimatS      = nan(N*T*S,nloc);
pimatSVec   = nan(N*T*S,nloc);
pimatS1Vec  = nan(N*T*S,1);
pimatS1Int  = nan(N*T*S,nloc);
pimatS1O    = nan(N*T*S,nloc);
pimatS1     = nan(N*T*S,nloc);
omegaSO     = nan(N*T*S,nloc);
omegaS      = nan(N*T*S,nloc);
omegaSVec   = nan(N*T*S,nloc);
omegaSInt  = nan(N*T*S,nloc);
pimatS1w   = nan(N*T*S,nloc,draws);
omegaSw    = nan(N*T*S,nloc,draws);

TmatS    = zeros(N*T*S,size(urate55,2));
TvecS    = calyrS(:)-2003;
for i=1:N*T*S
	if flagSl(i)==1
	    TmatS(i,TvecS(i)) = 1;
	else
	    TmatS(i,1)       = 1;
	end
end

tic
for k=1:nloc
	ftempS = flagS(:) & empFT_lagS(:)==0 & (k==LYS(:)-nloc*(LYS(:)>nloc)); 
	pimatS(ftempS,k) =   lambdervecRho(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,1,experS(ftempS),rho_hat_urate  ,0,bLemp,bLunemp);

	ftempS = flagS(:) & empFT_lagS(:)==1 & (k==LYS(:)-nloc*(LYS(:)>nloc));
	pimatS(ftempS,k) = 1-lambdervecRho(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,1,experS(ftempS),rho_hat_urate   ,0,bLemp,bLunemp);

	ftempS = flagS(:) & empFT_lagS(:)==0 & (k~=LYS(:)-nloc*(LYS(:)>nloc));
	pimatS(ftempS,k) =   lambdervecRho(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,0,experS(ftempS),rho_hat_urate,0,bLemp,bLunemp);

	ftempS = flagS(:) & empFT_lagS(:)==1 & (k~=LYS(:)-nloc*(LYS(:)>nloc));
	pimatS(ftempS,k) =   lambdervecRho(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,0,experS(ftempS),rho_hat_urate,0,bLemp,bLunemp);
	
	ftempS = flagS(:) & (k==LYS(:)-nloc*(LYS(:)>nloc));
	pimatS1(ftempS,1) =   lambdervecRho(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,0,1,experS(ftempS),rho_hat_urate,1,bLemp,bLunemp);
end

for k=1:nloc
	omegaS(:,k) = pimatS(:,k)./pimatS1(:,1);
end
if any(any(any(isnan(pimatS(flagSl,:)))))
	error('pi''s not set up correctly')
end
if any(any(any(isnan(pimatS1(flagSl,1)))))
	error('pi_t+1''s not set up correctly')
end
disp(['PimatS construction took ',num2str(toc/60),' minutes']); %takes 1.5 minutes without quadrature, 11 minutes with quadrature

%========================================================================
% Helpful matrices
%==========================================================================
calyr(isnan(calyr))=2004;
LYSi     = ones(N*T*S,1);
homeLocS = ones(N*T*S,1);
LhS      = zeros(N*T*S,J);
LwS      = zeros(N*T*S,J);
TmatS    = zeros(N*T*S,size(urate55,2));
TvecS    = calyrS(:)-2003;
tic;
for i=1:N*T*S
	if flagSl(i)==1
	    LYSi              = LYS(i);
	    l                 = LYSi-nloc*(LYSi>nloc);
	    homeLocS(i)       = l;
	    LhS(i,nloc+l)     = 1;
	    LwS(i,     l)     = 1;
	    TmatS(i,TvecS(i)) = 1;
	else
	    TmatS(i,1)        = 1;
	end
end
LhSS = LhS;
LwSS = LwS;
LhS = repmat(LhS,[1 1 draws]);
LwS = repmat(LwS,[1 1 draws]);
homeLocWS= repmat(homeLocS,[1 1 draws]);
pvEmpLS   = inlf_lagS(:)    ==1 & empFT_lagS(:)    ==1;
pvUnempLS = inlf_lagS(:)    ==1 & empFT_lagS(:)    ==0;
pvEmpS    = inlf_lagS(flagS)==1 & empFT_lagS(flagS)==1;
pvUnempS  = inlf_lagS(flagS)==1 & empFT_lagS(flagS)==0;
disp(['Helpful matrices loop took ',num2str(toc),' seconds']);

%========================================================================
% Flexible logit starting values
%==========================================================================
Adj = zeros(N*T,J);
if Beta>0
	%--------------------------------------------------------------------------
	% prepare data
	%--------------------------------------------------------------------------
	XS = [ones(N*T*S,2)];
	ZS = zeros(N*T*S,17,J);
	tic
	for j=1:nloc
	    k=j+nloc;
	    ZS(:,:,j) = [(pvEmpLS==1).*pimatS(:,j) (pvUnempLS==1).*pimatS(:,j) ewagevecRho(nloc,j,TmatS,TvecS,rho_hat_wage,0,experS(:),bwage)  birthLocS(:,j) birthDivS(:,j) switchervecS(nloc,LYS(:),j,ageS(:),typeS(:)) movervecS(nloc,LYS(:),j,distance,ageS(:),pvEmpLS,pvUnempLS,typeS(:))];
	    ZS(:,:,k) = [zeros(N*T*S,3)                                                                                                        birthLocS(:,k) birthDivS(:,k) switchervecS(nloc,LYS(:),k,ageS(:),typeS(:)) movervecS(nloc,LYS(:),k,distance,ageS(:),pvEmpLS,pvUnempLS,typeS(:))];
	end
	disp(['Flexible logit data matrices loop took ',num2str(toc/60),' minutes']);
	summarize(ZS(flagSl,:,1));

	%--------------------------------------------------------------------------
	% Load no-het estimates
	%--------------------------------------------------------------------------
	load([money,'_',num2str(nloc),'loc',sample,'.mat'],'b');
	bHet = [b(1:end-7);.2*rand(S-1,1);b(end-6:end);.2*rand(S-1,1)];
    bHetMat = bHet;
	save([money,'_',num2str(nloc),'loc',sample,'Het.mat'],'bHet');
	% Predict probabilities
	dat.Pc = pclogit1(bHet,YS,XS,ZS,size(XS,2),size(ZS,2),J,baseAlt);
	summarize(dat.Pc);
	dat.Pc([1 10 11 16 17],:)
	summarize(ZS(flagSl,:,1));
	NflexLogit       = sum(flag)
	NflexLogitUnique = numel(unique(ID(flagw)))
	% load flexLogitStartVal b
end

