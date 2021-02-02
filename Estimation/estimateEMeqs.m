%==========================================================================
% Estimate friction parameters
%==========================================================================
loc_dummies = [];
for j=setdiff(1:nloc,baseLoc);
	loc_dummies = cat(2,loc_dummies,ChoiceS(:)==j);
end

% logit of empFT conditional on inlf and empFT_lag
conditioner = inlfS(:)==1 & empFT_lagS(:)==1 & flagS(:)==1;
diffloc = ChoicelagS(:)~=ChoiceS(:) & ChoicelagS(:)~=ChoiceS(:)-nloc & ChoicelagS(:)-55~=ChoiceS(:);
Xtild = [loc_dummies UrateLagS(:) diffloc experS(:) experS(:).^2./100 typeS(:)];
Ytild = empFTS(:);
[bLemp,~,stats] = glmfit(Xtild(conditioner,:),Ytild(conditioner,:),'binomial', 'weights', PTypel(conditioner,:), 'link', 'logit');
dat.Pemp = glmval(bLemp,Xtild,'logit');
dat.Pemp(~conditioner) = 0;
dat.Pemp = reshape(dat.Pemp,[N T S]);
bLempMat = cat(2,bLempMat,bLemp);
bLempSEmat = cat(2,bLempSEmat,stats.se);

% logit of empFT conditional on inlf and ~empFT_lag
conditioner1 = inlfS(:)==1 & empFT_lagS(:)==0 & flagS(:)==1;
Xtild = [cat(2,loc_dummies(:,1:end-2),ChoiceS(:)==54 | ChoiceS(:)==55) UrateLagS(:) diffloc experS(:) experS(:).^2./100 typeS(:)];
Ytild = empFTS(:);
[bLunemp,~,stats] = glmfit(Xtild(conditioner1,:),Ytild(conditioner1,:),'binomial', 'weights', PTypel(conditioner1,:), 'link', 'logit');
dat.Punemp = glmval(bLunemp,Xtild,'logit');
dat.Punemp(~conditioner1) = 0;
dat.Punemp = reshape(dat.Punemp,[N T S]);
bLunemp = [bLunemp(1:end-(S-1)-4);bLunemp(end-(S-1)-4);bLunemp(end-(S-1)-3:end)]; 
bLunempMat = cat(2,bLunempMat,bLunemp);
bLempSEmat = cat(2,bLempSEmat,stats.se);

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
			indicter = [1; locer; 100*urate_lag55(j,t); 0; 0; 0; 0]';
			delta_hat(j,t)   = 1./(1+exp(indicter*bLemp));
			lambda_hat(j,t)  = exp(indicter*bLunemp)./(1+exp(indicter*bLunemp));
			indicter = [1;locer; 100*urate_lag55(j,t); 1; 0; 0; 0]';
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
	pimatS(ftempS,k) =   lambdervecRhoS(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,1,experS(ftempS),typeS(ftempS),rho_hat_urate  ,0,bLemp,bLunemp);

	ftempS = flagS(:) & empFT_lagS(:)==1 & (k==LYS(:)-nloc*(LYS(:)>nloc));
	pimatS(ftempS,k) = 1-lambdervecRhoS(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,1,experS(ftempS),typeS(ftempS),rho_hat_urate   ,0,bLemp,bLunemp);

	ftempS = flagS(:) & empFT_lagS(:)==0 & (k~=LYS(:)-nloc*(LYS(:)>nloc));
	pimatS(ftempS,k) =   lambdervecRhoS(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,0,experS(ftempS),typeS(ftempS),rho_hat_urate,0,bLemp,bLunemp);

	ftempS = flagS(:) & empFT_lagS(:)==1 & (k~=LYS(:)-nloc*(LYS(:)>nloc));
	pimatS(ftempS,k) =   lambdervecRhoS(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,empFT_lagS(ftempS)==1,0,experS(ftempS),typeS(ftempS),rho_hat_urate,0,bLemp,bLunemp);

	ftempS = flagS(:) & (k==LYS(:)-nloc*(LYS(:)>nloc));
	pimatS1(ftempS,1) =   lambdervecRhoS(nloc,k,TmatS(ftempS,:),TvecS(ftempS),urate_lag55,0,1,experS(ftempS),typeS(ftempS),rho_hat_urate,1,bLemp,bLunemp);
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

%==========================================================================
% Set up wage regression matrices
%==========================================================================
loc_dummiesS = zeros(N*T*S,nloc-1);
for j=setdiff(1:nloc,baseLoc);
	loc_dummiesS(:,j-1) = ChoiceS(:)==j;
end

time_dummiesS = zeros(N*T*S,9);
for k=2005:2013
   time_dummiesS(:,k-2004) = calyrS(:)==k;
end

loc_time_dummiesS = zeros(N*T*S,(nloc-1)*9);
for j=setdiff(1:nloc,baseLoc);
	for k=2005:2013
	    loc_time_dummiesS(:,(j-2)*9+k-2004) = ChoiceS(:)==j & calyrS(:)==k;
	end
end

%==========================================================================
% Estimate wage parameters
%==========================================================================
% construct flag for wages that will actually be used
wageflagStempl = wageflagS(:)==1 & ChoiceS(:)<=nloc & flagS(:) & empFTS(:)==1; % Note the addition that employmentFT must hold!!
size(wageflagStemp)
% [bwage  ,sebhatw  ] = lscov([ones(sum(wageflagtemp(:)),1) loc_dummies(wageflagtemp,:) exper(wageflagtemp) exper(wageflagtemp).^2./100],wagevar(wageflagtemp));
% [bwage   sebhatw  ]

[bwager  ,sebhatwr  , mser] = lscov([ones(sum(wageflagStemp(:)),1) loc_dummiesS(wageflagStempl,:) time_dummiesS(wageflagStempl,:) loc_time_dummiesS(wageflagStempl,1:end-18) experS(wageflagStemp) experS(wageflagStemp).^2./100 typeS(wageflagStemp)==1],dat.wage(wageflagStemp),PTypeS(wageflagStemp));
[bwager   sebhatwr  ];

[bwagertest  ,sebhatwr  , msertest] = lscov(dat.Xw(wageflagStempl,:),dat.wage(wageflagStemp),PTypeS(wageflagStemp));
[bwagertest   sebhatwr  ];

optionswage = optimset('Disp','iter','LargeScale','on','maxiter',1e8,'maxfuneval',1e8,'TolX',1e-6,'Tolfun',1e-6,'DerivativeCheck','off','GradObj','on','FinDiffType','central');
[parmwage,~,~,~,~,hwage] = fminunc('normalMLEw',[bwager;sqrt(mser)],optionswage,[],dat.wage(wageflagStemp),dat.Xw(wageflagStempl,:),PTypeS(wageflagStemp));

bwage = parmwage(1:end-1);
sebhatw = sqrt(diag(eye(size(hwage))/hwage));
sebhatw = sebhatw(1:end-1);
parms.wageBeta = bwage;
parms.wageSig  = parmwage(end);
bwageMat = cat(2,bwageMat,bwage);

NTSwage = sum(wageflagStemp(:))
Nwage  = numel(unique(IDS(wageflagStemp)))
[bwage sebhatw];
wHat      = zeros(nloc,10);
wHatOther = zeros(nloc,10);
for j=1:nloc
	for t=1:10
	    wHat(j,t) = ewageS(nloc,j,t,0,0,bwage)-bwage(1);
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
[bARw sqrt(diag(inv(full(hARw))))];
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
[cor,Pcor]=corr([resid_wageOld(:) resid_urate(:)]);
[cor,Pcor]=corr([resid_wage(:) resid_urate(:)]);
[cor,Pcor]=corr([resid_wageUncool(:) resid_urate(:)]);

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

[ARcor,pARcor]=corr(mastershocker);
ARcov=cov(mastershocker);
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
disp(['Helpful matrices loop took ',num2str(toc/60),' minutes']);

%========================================================================
% Estimate flexible logit
%==========================================================================
Adj = zeros(N*T*S,J);
if Beta>0
	%--------------------------------------------------------------------------
	% prepare data
	%--------------------------------------------------------------------------
	XS = [ones(N*T*S,2)];
	ZS = zeros(N*T*S,17,J);
	tic
	for j=1:nloc
	    k=j+nloc;
	    ZS(:,:,j) = [(pvEmpLS==1).*pimatS(:,j) (pvUnempLS==1).*pimatS(:,j) ewagevecRhoS(nloc,j,TmatS,TvecS,rho_hat_wage,0,experS(:),typeS(:),bwage)  birthLocS(:,j) birthDivS(:,j) switchervecS(nloc,LYS(:),j,ageS(:),typeS(:)) movervecS(nloc,LYS(:),j,distance,ageS(:),pvEmpLS,pvUnempLS,typeS(:))];
	    ZS(:,:,k) = [zeros(N*T*S,3)                                                                                                                  birthLocS(:,k) birthDivS(:,k) switchervecS(nloc,LYS(:),k,ageS(:),typeS(:)) movervecS(nloc,LYS(:),k,distance,ageS(:),pvEmpLS,pvUnempLS,typeS(:))];
	end
	disp(['Flexible logit data matrices loop took ',num2str(toc/60),' minutes']);
	summarize(ZS(flagSl,:,1));

	%--------------------------------------------------------------------------
	% estimate mlogit
	%--------------------------------------------------------------------------
	% restriction matrix
	% restrMat1(1 ,:) = [1 7  1 1 0]; % amenities loc 2 (work)
	% restrMat1(2 ,:) = [2 0  0 0 0]; % unemp benefits loc 2 (work)
	% restrMat1(3 ,:) = [3 9  1 1 0]; % amenities loc 3 (work)
	% restrMat1(4 ,:) = [4 0  0 0 0]; % unemp benefits loc 3 (work)
	% restrMat1(5 ,:) = [5 0  0 0 0]; % amenities loc 1 (home)
	% restrMat1(6 ,:) = [6 10 1 1 0]; % unemp benefits loc 1 (home)
	% restrMat1(7 ,:) = [8 10 1 1 0]; % unemp benefits loc 2 (home)
	% restrMat1

	nbstruc = 2;
	restrMat1 = nan(1,5);
	k = 0;
	for j=1:(2*nloc-1)*nbstruc
	    % Amenity is equal in all locations across work/home decisions
	    if (mod(j,nbstruc)==1) && j<(nloc-1)*nbstruc
	        k = k+1;
	        if ceil(j/nbstruc)>=baseLoc
	        restrMat1(k,:) = [j j+nloc*nbstruc 1 1 0];
	        elseif ceil(j/nbstruc)<baseLoc
	        restrMat1(k,:) = [j j+(nloc-1)*nbstruc 1 1 0];
	        end
	    % Amenity is equal to zero in base location
	    elseif (mod(j,nbstruc)==1) && j==(nloc-1)*nbstruc+(baseLoc-1)*nbstruc+1
	        k = k+1;
	    restrMat1(k,:) = [j 0 0 0 0];
	    % "Unemp benefits" (a.k.a. "work costs") equal 0 in all locations within work decisions
	    elseif (mod(j,nbstruc)==nbstruc-2) && j<=(nloc-1)*nbstruc
	        k = k+1;
	    restrMat1(k,:) = [j 0 0 0 0];
	    % "Unemp benefits" (a.k.a. "work costs") equal across all locations within home decisions
	    elseif (mod(j,nbstruc)==nbstruc-2) && j>(nloc-1)*nbstruc+2
	        k = k+1;
	    restrMat1(k,:) = [j (nloc-1)*nbstruc+2 1 1 0];
	    end
	end
	%restrMat1
	
	options = optimset('Disp','iter','LargeScale','on','maxiter',1e8,'maxfuneval',1e8,'TolX',1e-1,'Tolfun',1e-1,'DerivativeCheck','off','GradObj','on','FinDiffType','central');
	optionstest = optimset('Disp','iter','LargeScale','on','maxiter',1,'maxfuneval',1,'GradObj','off','DerivativeCheck','off');
	load([money,'_',num2str(nloc),'loc',sample,'Het.mat'],'bHet');
	startval = bHet;
	%restrMat1
	[bHet,l,e,o,g,h] =  fminunc('clogitW',startval,options,restrMat1,YS(flagSl),XS(flagSl,:),ZS(flagSl,:,:),baseAlt,PTypel(flagSl));
	[bHet,invH] = applyRestr(restrMat1,bHet,h);
	[bHet sqrt(diag(invH))];
    bHetMat = cat(2,bHetMat,bHet);
	save([money,'_',num2str(nloc),'loc',sample,'Het.mat'],'bHet');
	dat.Pc = pclogit1(bHet,YS,XS,ZS,size(XS,2),size(ZS,2),J,baseAlt);
	dat.Pc(~flagSl,:) = 0;
	summarize(dat.Pc);
	dat.Pc([1 10 11 16 17],:)
	NflexLogit       = sum(flagSl)
	NflexLogitUnique = numel(unique(IDS(flagS)))
	
	save intermediateResults bwageMat bLunempMat bLempMat bHetMat qMat likevec
end	

