if strcmp(time,'trimest')==1
	error('You shouldn''t be using trimesterly data!!');
    Beta = Beta.^(1/3); % convert annual discount factor to trimesterly
end

if Beta==0
    delete(['estimationFrict_',money,'_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['estimationFrict_',money,'_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
else                               
    delete(['estimationFrict_',money,'_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
    diary (['estimationFrict_',money,'_',num2str(nloc),'loc',sample,'nonstat_beta',num2str(10*Beta),'.diary']);
end

tic;

%==========================================================================
% Read in the data
%==========================================================================
load [REDACTED]sippCombinedNHWmaleAnnual.mat

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

if nloc==3
    [N,T] = size(choice6);
    Choice = choicefrict6;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict6_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist6;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag3(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate3    (j,t);
        end
    end
elseif nloc==27
    [N,T] = size(choice54);
    Choice = choicefrict54;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict54_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist54;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag27(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate27    (j,t);
        end
    end
elseif nloc==28
    [N,T] = size(choice56);
    Choice = choicefrict56;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict56_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist56;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag28(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate28    (j,t);
        end
    end
elseif nloc==29
    [N,T] = size(choice58);
    Choice = choicefrict58;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict58_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist58;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag29(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate29    (j,t);
        end
    end
elseif nloc==48
    [N,T] = size(choice96);
    Choice = choicefrict96;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict96_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist96;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag48(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate48    (j,t);
        end
    end
elseif nloc==49
    [N,T] = size(choice98);
    Choice = choicefrict98;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict98_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist98;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag49(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate49    (j,t);
        end
    end
elseif nloc==50
    [N,T] = size(choice100);
    Choice = choicefrict100;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict100_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist100;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag50(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate50    (j,t);
        end
    end
elseif nloc==53
    [N,T] = size(choice106);
    Choice = choicefrict106;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict106_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist106;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag53(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate53    (j,t);
        end
    end
elseif nloc==54
    [N,T] = size(choice108);
    Choice = choicefrict108;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict108_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist108;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag54(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate54    (j,t);
        end
    end
elseif nloc==55
    [N,T] = size(choice110);
    Choice = choicefrict110;
    Choice(choiceflag==0)=0;
    Choicelag = choicefrict110_lag;
    Choicelag(choiceflag==0)=0;
    distance = dist110;
    UrateLag = nan(N,T);
    Urate    = nan(N,T);
    for j=1:nloc
        for t=1:10
            UrateLag((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate_lag55(j,t);
            Urate   ((Choice==j | Choice==j+nloc) & calyr==2003+t) = 100*urate55    (j,t);
        end
    end
end
Y  = Choice(:);
LY = Choicelag(:);

if strcmp(sample,'HS')==1
    flag  = Choice(:)>0 & Choicelag(:)>0 & collgrad(:)==0;
    flagw = Choice>0 & Choicelag>0 & collgrad==0;
    disp('tabulation of Y');
	tabulate(Y(flag));
    disp('tabulation of LY');
	tabulate(LY(flag));
elseif strcmp(sample,'BA')==1
    flag  = Choice(:)>0 & Choicelag(:)>0 & collgrad(:)==1;
    flagw = Choice>0 & Choicelag>0 & collgrad==1;
    disp('tabulation of Y');
	tabulate(Y(flag));
    disp('tabulation of LY');
	tabulate(LY(flag));
end

draws=20;
J = 2*nloc;
ID = [1:N]'*ones(1,T);
baseAlt = 1;
if baseAlt>nloc
    baseLoc = baseAlt-nloc;
else
    baseLoc = baseAlt;
end

% assert experience variable
exper = exper_IRS;
exper_increment = 1*annual + (1/3)*(1-annual)

%==========================================================================
% Initialize permanent unobserved heterogeneity parameters
%==========================================================================
S         = 2; % number of types
alph      = 0.5; % threshold for starting value noise
PType     = rand(N,S);
PType     = PType./repmat(sum(PType,2),[1 S]);
PTypeS    = zeros(N,T,S);
for s=1:S
	PTypeS(:,:,s) = repmat(PType(:,s),[1 T]);
end
PTypel    = PTypeS(:);
oPType    = zeros(size(PType));
criter    = norm(PType(:)-oPType(:), Inf);
iteration = 1;
prior     = [.6 .4];
likevec   = [];
qMat      = [];

% replicate data (N x S x T) for EM algorithm:
IDS          = repmat(ID,[1 1 S]);            %kron(ones(S,1),ID);
inlfS        = repmat(inlf,[1 1 S]);          %kron(ones(S,1),inlf);
inlf_lagS    = repmat(inlf_lag,[1 1 S]);      %kron(ones(S,1),inlf_lag);
empFTS       = repmat(empFT,[1 1 S]);         %kron(ones(S,1),empFT);
empFT_lagS   = repmat(empFT_lag,[1 1 S]);     %kron(ones(S,1),empFT_lag);
flagS        = repmat(flagw,[1 1 S]);         %kron(ones(S,1),flagw);
wageflagS    = repmat(wageflag,[1 1 S]);      %kron(ones(S,1),wageflag);
ChoicelagS   = repmat(Choicelag,[1 1 S]);     %kron(ones(S,1),Choicelag);
ChoiceS      = repmat(Choice,[1 1 S]);        %kron(ones(S,1),Choice);
UrateLagS    = repmat(UrateLag,[1 1 S]);      %kron(ones(S,1),UrateLag);
ageS         = repmat(age,[1 1 S]);           %kron(ones(S,1),age);
experS       = repmat(exper,[1 1 S]);         %kron(ones(S,1),exper);
calyrS       = repmat(calyr,[1 1 S]);         %kron(ones(S,1),calyr);
birthLoc110S = repmat(birthLoc110,[1 1 1 S]); %kron(ones(S,1),calyr);
birthDiv110S = repmat(birthDiv110,[1 1 1 S]); %kron(ones(S,1),calyr);
birthLocS    = reshape(permute(birthLoc110S,[1 3 4 2]),N*T*S,110);
birthDivS    = reshape(permute(birthDiv110S,[1 3 4 2]),N*T*S,110);
YS         = ChoiceS(:);
LYS        = ChoicelagS(:);
flagSl     = flagS(:);

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

if strcmp(money,'wage')==1
	dat.wage = repmat(lnWageHr,[1 1 S]);
elseif strcmp(money,'earn')==1
	dat.wage = repmat(lnearnfinalJb,[1 1 S]);
end

wageflagStemp = wageflagS==1 & ChoiceS<=nloc & flagS & empFTS==1;
typeS     = cat(3,ones(N,T),zeros(N,T));

% data in N x T x S format
dat.inlfS         = inlfS;
dat.empFTS        = empFTS;
dat.empFT_lagS    = empFT_lagS;
dat.ChoicelagS    = ChoicelagS;
dat.ChoiceS       = ChoiceS;
dat.flagS         = flagS;
dat.typeS         = typeS;
dat.N             = N;
dat.T             = T;
dat.S             = S;
dat.J             = J;
dat.wageflagStemp = wageflagStemp;
%dat.Xw = [ones(sum(wageflagStemp(:)),1) loc_dummiesS(wageflagStemp,:) time_dummiesS(wageflagStemp,:) loc_time_dummiesS(wageflagStemp,1:end-18) experS(wageflagStemp) experS(wageflagStemp).^2./100 typeS(wageflagStemp)];
dat.Xw = [ones(length(wageflagStemp(:)),1) loc_dummiesS time_dummiesS loc_time_dummiesS(:,1:end-18) experS(:) experS(:).^2./100 typeS(:)];

skip_est = true;
first_time_est = false;
if skip_est == false	
	if first_time_est == true
		%while iteration<3
		while criter>1e-5
			oPType = PType;
	
			% load starting values
			if iteration==1
				tstartval = tic;
				loadstartvalues
				telapsed = toc(tstartval);
				disp(['Time spent loading starting values: ',num2str(telapsed/60),' minutes']);
			end
			
			% compute full likelihood
			tic
			full_like = likecalc(dat,parms);
			disp(['Time spent calculating full likelihood: ',num2str(toc/60),' minutes']);
	
			% Update q's
			[prior,PType,jointlike] = typeprob(prior,full_like);
			disp(['Full likelihood after updating q''s is ',num2str(jointlike)]);
			disp(['Pr(type==1) is ',num2str(prior(1))]);
			likevec = [likevec; jointlike];
			qMat = cat(3,qMat,PType);

			% update relevant q matrices
			for s=1:S
				PTypeS(:,:,s) = repmat(PType(:,s),[1 T 1]);
			end
			PTypel = reshape(PTypeS,[N*T*S 1]);

			% estimation
			estimateEMeqs

			criter = norm(PType(:)-oPType(:),Inf);
			disp(['Iteration is ',num2str(iteration)]);
			disp(['EM criterion is ',num2str(criter)]);
			iteration = iteration+1;
		end
	else
		while criter>1e-5
			oPType = PType;
	
			% load starting values
			if iteration==1
				tstartval = tic;
				loadstartvalues
				telapsed = toc(tstartval);
				disp(['Time spent loading starting values: ',num2str(telapsed/60),' minutes']);
				
				load intermediateResultsFinished88iters;
				prior =[.5 .5];
				PType = qMat(:,:,end);
				for s=1:S
					PTypeS(:,:,s) = repmat(PType(:,s),[1 T 1]);
				end
				PTypel = reshape(PTypeS,[N*T*S 1]);
				bHet = bHetMat(:,end);
				bLemp = bLempMat(:,end);
				bLunemp = bLunempMat(:,end);
				bwage = bwageMat(:,end);
				iteration = size(qMat,3);
			end
		
			% compute full likelihood
			tic
			full_like = likecalc(dat,parms);
			disp(['Time spent calculating full likelihood: ',num2str(toc/60),' minutes']);
	
			% Update q's
			[prior,PType,jointlike] = typeprob(prior,full_like);
			disp(['Full likelihood after updating q''s is ',num2str(jointlike)]);
			disp(['Pr(type==1) is ',num2str(prior(1))]);
			likevec = [likevec; jointlike];
			qMat = cat(3,qMat,PType);

			% update relevant q matrices
			for s=1:S
				PTypeS(:,:,s) = repmat(PType(:,s),[1 T 1]);
			end
			PTypel = reshape(PTypeS,[N*T*S 1]);

			% estimation
			estimateEMeqs

			criter = norm(PType(:)-oPType(:),Inf);
			disp(['Iteration is ',num2str(iteration)]);
			disp(['EM criterion is ',num2str(criter)]);
			iteration = iteration+1;
		end
	end
	save -v7.3 allBefDebuggingCCPs
end
est_FV_terms = false;
if est_FV_terms==true
load allBefDebuggingCCPs
%pooler = parpool(8);
Adj = zeros(N*T*S,J);
if Beta>0    
    disp('Starting in on the CCPs now');
    % ========================================================================
    % Construct CCP's from data and mlogit coefficients
    % ==========================================================================
	XS = [ones(N*T*S,2)];
	YS(~flagSl) = J;
    ZtempS = zeros(N*T*S,17,J,draws);
    shocker = zeros(N*T*S,J,draws);
    P = rand(N*T*S,J,draws);
    % Need to loop over Z matrix in flexible logit to cover each set of states
    tic
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov);
    end
    % Loop 1a: Pr({0,lprime} t+1 | {jprime,lprime} t and employed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives
        %parfor j=1:J
        for j=1:J
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockd1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
                if j==lp
                pim = 1-lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,1,1,experS(flagS)+jp,typeS(flagS),rho_hat_urate,1,shockd1,bLemp,bLunemp);
                size(pim);
                ZtempS(flagSl,:,j,:) = [pim zeros(sum(flagSl),1,draws) ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)+jp,typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ageS(flagS)+1,ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
                else
                pim = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,1,0,experS(flagS)+jp,typeS(flagS),rho_hat_urate,1,shockd1,bLemp,bLunemp);
                ZtempS(flagSl,:,j,:) = [pim zeros(sum(flagSl),1,draws) ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)+jp,typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ageS(flagS)+1,ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
                end
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ageS(flagS)+1,ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -Beta.*(pimatS(:,k)).*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
        else
            Adj(:,k) = -Beta.*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
        end
    end
    disp(['loop 1a finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');

    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov);
    end
    % Loop 1b: Pr({0,lprime} t+1 | {jprime,lprime} t and unemployed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives
        %parfor j=1:J 
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
                if j==lp
                pim = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,1,experS(flagS),typeS(flagS),rho_hat_urate,1,shockl1,bLemp,bLunemp);
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS),typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ageS(flagS)+1,zeros(sum(flagSl),1),ones(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
                else
                pim = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,0,experS(flagS),typeS(flagS),rho_hat_urate,1,shockl1,bLemp,bLunemp);
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS),typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ageS(flagS)+1,zeros(sum(flagSl),1),ones(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
                end
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ageS(flagS)+1,zeros(sum(flagSl),1),ones(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -Beta.*(1-pimatS(:,k)).*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
        end
    end
    disp(['loop 1b finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');

    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov2);
    end
    % Loop 2a: Pr({0,l} t+2 | {0,lprime} t+1, {jprime,lprime} t and employed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=nloc+lp (the t+1 decision)
        % if lp==l, then lambda
        % otherwise, lambda_u
        %parfor j=1:J
        for j=1:J
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
                entry1 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,1,experS(flagS)+jp,typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,0,experS(flagS)+jp,typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
				pim = entry1.*(homeLocWS(flagSl,:,:)==j)+entry2.*(homeLocWS(flagSl,:,:)~=j);				
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)+jp,typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -(Beta^2).*(pimatS(:,k)).*((1/draws)*sum(log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        else
            Adj(:,k) = -(Beta^2).*((1/draws)*sum(log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        end
    end
    disp(['loop 2a finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');

    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov2);
    end
    % Loop 2b: Pr({0,l} t+2 | {0,lprime} t+1, {jprime,lprime} t and unemployed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=nloc+lp (the t+1 decision)
        %parfor j=1:J
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
                entry1 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,1,experS(flagS),typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,0,experS(flagS),typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
				pim = entry1.*(homeLocWS(flagSl,:,:)==j)+entry2.*(homeLocWS(flagSl,:,:)~=j);
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS),typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -(Beta^2).*(1-pimatS(:,k)).*((1/draws)*sum(log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        end
    end
    disp(['loop 2b finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');
    
    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov);
    end
    for i=1:N*T*S %:numel(exper) % loop over individual-time observations
        if flagSl(i)==1
            LYSi = LYS(i);
            ti  = calyrS(i)-2003;
            l = LYSi-nloc*(LYSi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
			pimatS1w(i,1,:) = lambdervecRhoIntS(nloc,l,TmatS(i,:),TvecS(i),urate_lag55,0,1,experS(i),typeS(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
            for jj=1:nloc
                omegaSw(i,jj,:) = repmat(pimatS(i,jj),[1 1 draws])./pimatS1w(i,1,:);
            end
        end
    end
    % Loop 3a: Pr({jprime,l} t+1 | {0,l} t) -- weight by omega(:,k) [employment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives, fixing k=nloc+1 (the t decision)
        %parfor j=1:J 
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
        		entry1 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,1,experS(flagS),typeS(flagS),rho_hat_urate,1,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,0,experS(flagS),typeS(flagS),rho_hat_urate,1,shockl1,bLemp,bLunemp);
        		pim = entry1.*(homeLocWS(flagSl,:,:)==j)+entry2.*(homeLocWS(flagSl,:,:)~=j);
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS),typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+homeLocS(flagSl),j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+homeLocS(flagSl),j,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+homeLocS(flagSl),j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+homeLocS(flagSl),j,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = Beta.*((1/draws)*sum(omegaSw(:,k,:).*log(max(P.*LwS,[],2)),3)) + Adj(:,k);
        else
            Adj(:,k) = Beta.*((1/draws)*sum(log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        end
    end
    disp(['loop 3a finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');

    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov);
    end
    for i=1:N*T*S %:numel(exper) % loop over individual-time observations
        if flagSl(i)==1
            LYSi = LYS(i);
            ti  = calyrS(i)-2003;
            l = LYSi-nloc*(LYSi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
			pimatS1w(i,1,:) = lambdervecRhoIntS(nloc,l,TmatS(i,:),TvecS(i),urate_lag55,0,1,experS(i),typeS(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
            for jj=1:nloc
                omegaSw(i,jj,:) = repmat(pimatS(i,jj),[1 1 draws])./pimatS1w(i,1,:);
            end
        end
    end
    % Loop 3b: Pr({jprime,l} t+1 | {0,l} t) -- weight by 1-omega(:,k) [unemployment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives, fixing k=nloc+1 (the t decision)
        %parfor j=1:J 
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
        		entry1 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,1,experS(flagS),typeS(flagS),rho_hat_urate,1,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,0,experS(flagS),typeS(flagS),rho_hat_urate,1,shockl1,bLemp,bLunemp);
				pim = entry1.*(homeLocWS(flagSl,:,:)==j)+entry2.*(homeLocWS(flagSl,:,:)~=j);
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS),typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+homeLocS(flagSl),j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+homeLocS(flagSl),j,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+homeLocS(flagSl),j,ageS(flagS)+1,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+homeLocS(flagSl),j,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = Beta.*((1/draws)*sum((1-omegaSw(:,k,:)).*log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        end
    end
    disp(['loop 3b finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');

    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov2);
    end
    for i=1:N*T*S %:numel(exper) % loop over individual-time observations
        if flagSl(i)==1
            LYSi = LYS(i);
            ti  = calyrS(i)-2003;
            l = LYSi-nloc*(LYSi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
			pimatS1w(i,1,:) = lambdervecRhoIntS(nloc,l,TmatS(i,:),TvecS(i),urate_lag55,0,1,experS(i),typeS(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
            for jj=1:nloc
                omegaSw(i,jj,:) = repmat(pimatS(i,jj),[1 1 draws])./pimatS1w(i,1,:);
            end
        end
    end
    % Loop 4a: Pr({0,l} t+2 | {jprime,l} t+1, {0,l} t) [employment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
        %parfor j=1:J 
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockd1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
        		entry1 = 1-lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,1,1,experS(flagS)+jp,typeS(flagS),rho_hat_urate,2,shockd1,bLemp,bLunemp);
        		entry2 =   lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,1,0,experS(flagS)+jp,typeS(flagS),rho_hat_urate,2,shockd1,bLemp,bLunemp);
        		pim = entry1.*(homeLocWS(flagSl,:,:)==j)+entry2.*(homeLocWS(flagSl,:,:)~=j);
                ZtempS(flagSl,:,j,:) = [pim zeros(sum(flagSl),1,draws) ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)+jp,typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,distance,ageS(flagS)+2,ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,distance,ageS(flagS)+2,ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = (Beta^2).*((1/draws)*sum(pimatS1w(:,1,:).*omegaSw(:,k,:).*log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        else
            Adj(:,k) = (Beta^2).*((1/draws)*sum(log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        end
    end
    disp(['loop 4a finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');
    
    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov2);
    end
    for i=1:N*T*S %:numel(exper) % loop over individual-time observations
        if flagSl(i)==1
            LYSi = LYS(i);
            ti  = calyrS(i)-2003;
            l = LYSi-nloc*(LYSi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
			pimatS1w(i,1,:) = lambdervecRhoIntS(nloc,l,TmatS(i,:),TvecS(i),urate_lag55,0,1,experS(i),typeS(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
            for jj=1:nloc
                omegaSw(i,jj,:) = repmat(pimatS(i,jj),[1 1 draws])./pimatS1w(i,1,:);
            end
        end
    end
    % Loop 4b: Pr({0,l} t+2 | {jprime,l} t+1, {0,l} t) [unemployment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
        %parfor j=1:J 
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
        		entry1 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,1,experS(flagS),typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,0,experS(flagS),typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		pim = entry1.*(homeLocWS(flagSl,:,:)==j)+entry2.*(homeLocWS(flagSl,:,:)~=j);
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS),typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),ones(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+homeLocS(flagSl),j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),ones(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = (Beta^2).*((1/draws)*sum((1-pimatS1w(:,1,:)).*omegaSw(:,k,:).*log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        end
    end
    disp(['loop 4b finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');

    
    for dd=1:draws
        shocker(flagSl,:,dd) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov2);
    end
    for i=1:N*T*S %:numel(exper) % loop over individual-time observations
        if flagSl(i)==1
            LYSi = LYS(i);
            ti  = calyrS(i)-2003;
            l = LYSi-nloc*(LYSi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
			pimatS1w(i,1,:) = lambdervecRhoIntS(nloc,l,TmatS(i,:),TvecS(i),urate_lag55,0,1,experS(i),typeS(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
            for jj=1:nloc
                omegaSw(i,jj,:) = repmat(pimatS(i,jj),[1 1 draws])./pimatS1w(i,1,:);
            end
        end
    end
    % Loop 4c: Pr({0,l} t+2 | {0,l} t+1, {0,l} t) [NILF counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
        %parfor j=1:J 
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flagSl,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flagSl,(j-1)*2+2,:));
        		entry1 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,1,experS(flagS),typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),urate_lag55,0,0,experS(flagS),typeS(flagS),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		pim = entry1.*(homeLocWS(flagSl,:,:)==j)+entry2.*(homeLocWS(flagSl,:,:)~=j);
				ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),1,draws) pim ewagevecRhoIntS(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS),typeS(flagS),bwage,shockw1) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+homeLocS(flagSl),j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+homeLocS(flagSl),j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            else
                ZtempS(flagSl,:,j,:) = [zeros(sum(flagSl),3,draws) repmat(birthLocS(flagSl,j),[1 1 draws]) repmat(birthDivS(flagSl,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+homeLocS(flagSl),j,ageS(flagS)+2,typeS(flagS)),[1 1 draws]) repmat(movervecS(nloc,nloc+homeLocS(flagSl),j,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flagSl,:,dd) = pclogit(bHet,YS(flagSl),XS(flagSl,:),ZtempS(flagSl,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = (Beta^2).*((1/draws)*sum((1-omegaSw(:,k,:)).*log(max(P.*LhS,[],2)),3)) + Adj(:,k);
        end
    end
    disp(['loop 4c finished']);
    save(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');
    
    disp(['CCP loop took ',num2str(toc/3600),' hours']);
	[flagSl(1:10) LY(1:10) Adj(1:10,:)]
    % save(['FVFrict',money,num2str(nloc),'loc',sample,'.mat'],'Adj');
    summarize(Adj(flagSl,:));
    save -v7.3 allbefFinalFminunc
end
load(['FVFrictHet',money,num2str(nloc),'loc',sample,'.mat'],'Adj');
end
load allbefFinalFminunc
%==========================================================================
% Estimate structural model
%==========================================================================
%---------------------------------------------------------------------------
% prepare data
%---------------------------------------------------------------------------
% Xtilde = (1-Beta).*[ones(N*T,2)];
% Z = [prev_diff_emp(:) prev_diff_loc(:)];
%summarize(pimatS(flagSl,:));
%wageTest = nan(length(experS(:)),nloc);
%for j=1:nloc
%	wageTest(flagSl,j) = ewagevecRho(nloc,j,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,0,experS(flagS),bwage);
%end
%summarize(wageTest(flagSl,:));
tic;
Zloc = zeros(N*T*S,nloc,J);
Ztilde    = zeros(N*T*S,17,J);
Ztilde0eA = zeros(N*T*S,17,J);
Ztilde1eA = zeros(N*T*S,17,J);
Ztilde2eA = zeros(N*T*S,17,J);
Ztilde0uA = zeros(N*T*S,17,J);
Ztilde1uA = zeros(N*T*S,17,J);
Ztilde2uA = zeros(N*T*S,17,J);
Ztilde0B  = zeros(N*T*S,17,J);
Ztilde1eB = zeros(N*T*S,17,J);
Ztilde2eB = zeros(N*T*S,17,J);
Ztilde1uB = zeros(N*T*S,17,J);
Ztilde2uB = zeros(N*T*S,17,J);
Ztilde1nB = zeros(N*T*S,17,J);
Ztilde2nB = zeros(N*T*S,17,J);
shocker   = nan(N*T*S,2*nloc);
shocklee  = nan(N*T*S,1);
for dd=1:draws
    shocker(flagSl,:,1) = mvnrnd(zeros(sum(flagSl),2*nloc),ARcov);
    for i=1:N*T*S %:numel(exper) % loop over individual-time observations
        if flagSl(i)==1
            LYSi = LYS(i);
            ti  = calyrS(i)-2003;
            l = LYSi-nloc*(LYSi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
            for jj=1:nloc
                pimatS1(i,1) = squeeze(lambdervecRhoIntS(nloc,l,TmatS(i,:),TvecS(i),urate_lag55,0,1,experS(i),typeS(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp)); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
                omegaS(i,jj) = pimatS(i,jj)./pimatS1(i,1);
            end
            shocklee(i) = squeeze(shocker(i,(l-1)*2+1,:));
        end
    end
    for j=1:nloc
        lp = j;
        k=j+nloc;
        Ztilde0eA(flagSl,:,j) = repmat((  pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,j                    ,experS(flagS)) ewagevecRhoIntS(nloc,j                    ,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,0,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,j                          ) birthDivS(flagSl,j                          ) switchervecS(nloc,LYS(flagS)           ,j                    ,ageS(flagS)+0,typeS(flagS)) movervecS(nloc,LYS(flagS)           ,j                    ,distance,ageS(flagS)+0,              pvEmpS,            pvUnempS,typeS(flagS))];
        Ztilde0eA(flagSl,:,k) = repmat((  pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,k                    ,experS(flagS)) ewagevecRhoIntS(nloc,k                    ,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,0,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,k                          ) birthDivS(flagSl,k                          ) switchervecS(nloc,LYS(flagS)           ,k                    ,ageS(flagS)+0,typeS(flagS)) movervecS(nloc,LYS(flagS)           ,k                    ,distance,ageS(flagS)+0,              pvEmpS,            pvUnempS,typeS(flagS))];
        Ztilde0uA(flagSl,:,j) = repmat((1-pimatS(flagSl,j)),[1 17])                     .*[1*ones(sum(flagSl),1) benefits(nloc,k                    ,experS(flagS)) ewagevecRhoIntS(nloc,k                    ,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,0,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,j                          ) birthDivS(flagSl,j                          ) switchervecS(nloc,LYS(flagS)           ,j                    ,ageS(flagS)+0,typeS(flagS)) movervecS(nloc,LYS(flagS)           ,j                    ,distance,ageS(flagS)+0,              pvEmpS,            pvUnempS,typeS(flagS))];
        Ztilde0uA(flagSl,:,k) = repmat((1-pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,k                    ,experS(flagS)) ewagevecRhoIntS(nloc,k                    ,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,0,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,k                          ) birthDivS(flagSl,k                          ) switchervecS(nloc,LYS(flagS)           ,k                    ,ageS(flagS)+0,typeS(flagS)) movervecS(nloc,LYS(flagS)           ,k                    ,distance,ageS(flagS)+0,              pvEmpS,            pvUnempS,typeS(flagS))];
        Ztilde1eA(flagSl,:,j) = repmat((  pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,lp+nloc              ,experS(flagS)) ewagevecRhoIntS(nloc,lp              +nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)+1,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,lp+nloc                    ) birthDivS(flagSl,lp+nloc                    ) switchervecS(nloc,j                    ,lp+nloc              ,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,j                    ,lp+nloc              ,distance,ageS(flagS)+1, ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde1eA(flagSl,:,k) = repmat((  pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,lp+nloc              ,experS(flagS)) ewagevecRhoIntS(nloc,lp              +nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,lp+nloc                    ) birthDivS(flagSl,lp+nloc                    ) switchervecS(nloc,k                    ,lp+nloc              ,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,k                    ,lp+nloc              ,distance,ageS(flagS)+1, ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde1uA(flagSl,:,j) = repmat((1-pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,lp+nloc              ,experS(flagS)) ewagevecRhoIntS(nloc,lp              +nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,lp+nloc                    ) birthDivS(flagSl,lp+nloc                    ) switchervecS(nloc,j                    ,lp+nloc              ,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,j                    ,lp+nloc              ,distance,ageS(flagS)+1,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];
        Ztilde1uA(flagSl,:,k) = repmat((1-pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,lp+nloc              ,experS(flagS)) ewagevecRhoIntS(nloc,lp              +nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) birthLocS(flagSl,lp+nloc                    ) birthDivS(flagSl,lp+nloc                    ) switchervecS(nloc,k                    ,lp+nloc              ,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,k                    ,lp+nloc              ,distance,ageS(flagS)+1,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];
        Ztilde2eA(flagSl,:,j) = repmat((  pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)+1,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];
        Ztilde2eA(flagSl,:,k) = repmat((  pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];
        Ztilde2uA(flagSl,:,j) = repmat((1-pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];
        Ztilde2uA(flagSl,:,k) = repmat((1-pimatS(flagSl,j)),[1 17])                     .*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,lp+nloc              ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];

        Ztilde0B(flagSl,:,j)  =                                                           [0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,0,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,LYS(flagS)           ,k                    ,ageS(flagS)+0,typeS(flagS)) movervecS(nloc,LYS(flagS)           ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+0,              pvEmpS,            pvUnempS,typeS(flagS))];
        Ztilde0B(flagSl,:,k)  =                                                           [0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,0,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,LYS(flagS)           ,k                    ,ageS(flagS)+0,typeS(flagS)) movervecS(nloc,LYS(flagS)           ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+0,              pvEmpS,            pvUnempS,typeS(flagS))];
        Ztilde1eB(flagSl,:,j) = repmat((  pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)     ,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)     ,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,shocklee(flagSl)                    ) max(birthLocS(flagSl,:).*LwSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LwSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)     ,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)     ,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde1eB(flagSl,:,k) = repmat((  pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde1uB(flagSl,:,j) = repmat((1-pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[1*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)     ,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)     ,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde1uB(flagSl,:,k) = repmat((1-pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde1nB(flagSl,:,j) = repmat((1-                     omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde1nB(flagSl,:,k) = repmat((1-                     omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,1,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+1,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+1,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde2eB(flagSl,:,j) = repmat((  pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)+1,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)     ,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)     ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2, ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde2eB(flagSl,:,k) = repmat((  pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2, ones(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde2uB(flagSl,:,j) = repmat((1-pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)     ,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)     ,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];
        Ztilde2uB(flagSl,:,k) = repmat((1-pimatS1(flagSl,1)).*(omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1), ones(sum(flagSl),1),typeS(flagS))];
        Ztilde2nB(flagSl,:,j) = repmat((1-                     omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        Ztilde2nB(flagSl,:,k) = repmat((1-                     omegaS(flagSl,j)),[1 17]).*[0*ones(sum(flagSl),1) benefits(nloc,homeLocS(flagSl)+nloc,experS(flagS)) ewagevecRhoIntS(nloc,homeLocS(flagSl)+nloc,TmatS(flagSl,:),TvecS(flagSl),rho_hat_wage,2,experS(flagS)  ,typeS(flagS),bwage,squeeze(shocker(flagSl,(j-1)*2+1,:))) max(birthLocS(flagSl,:).*LhSS(flagSl,:),[],2) max(birthDivS(flagSl,:).*LhSS(flagSl,:),[],2) switchervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,ageS(flagS)+2,typeS(flagS)) movervecS(nloc,homeLocS(flagSl)+nloc,homeLocS(flagSl)+nloc,distance,ageS(flagS)+2,zeros(sum(flagSl),1),zeros(sum(flagSl),1),typeS(flagS))];
        for i=1:N*T*S
            if flagSl(i)==1
                LYSi = LYS(i);
                l = LYSi-nloc*(LYSi>nloc);
                if j~=l
                    Zloc(i,j,j) =  (1+Beta);
                    Zloc(i,l,j) = -(1+Beta);
                end
                if k~=nloc+l
                    Zloc(i,j,k) =  (1+Beta);
                    Zloc(i,l,k) = -(1+Beta);
                end
                if k==nloc+l
                    Zloc(i,:,[l l+nloc]) = zeros(1,nloc,2);
                end
            end
        end
    end
    ZtildeTemp = Ztilde0eA+Ztilde0uA-Ztilde0B + Beta*(Ztilde1eA+Ztilde1uA-Ztilde1eB-Ztilde1uB-Ztilde1nB) + Beta.^2*(Ztilde2eA+Ztilde2uA-Ztilde2eB-Ztilde2uB-Ztilde2nB);
    Ztilde = ZtildeTemp+Ztilde;
end
ZS = cat(2,Zloc,(1/draws)*Ztilde);
clear Ztilde*
YS = ChoiceS(:);
assert(size(Zloc,2)==nloc,'Zloc incorrectly constructed');
clear Zloc
disp(['Structural data loop took ',num2str(toc/60),' minutes']);
if nloc==3
    save ZtesterShorter ZS
    summarize(ZS(flagSl,:,1));
    summarize(ZS(flagSl,:,J));
    % return
end
save -v7.3 structZ baseAlt YS flagSl ZS PTypel Adj N T S J IDS flagS
%--------------------------------------------------------------------------
% estimate structural likelihood
%--------------------------------------------------------------------------
% restriction matrix
restrMat(1,:) = [baseAlt 0  0 0 0]; % amenities in base location
% restrMat = []; % amenities in base location

options = optimset('Disp','iter','LargeScale','on','maxiter',1e8,'maxfuneval',1e8,'TolX',1e-6,'Tolfun',1e-6,'DerivativeCheck','off','GradObj','on','FinDiffType','central');
o4Nu=optimset('Disp','Iter','LargeScale','off','MaxFunEvals',0,'MaxIter',0,'TolX',1e-6,'Tolfun',1e-6,'GradObj','off','DerivativeCheck','off','FinDiffType','central');
o4An=optimset('Disp','Iter','LargeScale','off','MaxFunEvals',0,'MaxIter',0,'TolX',1e-6,'Tolfun',1e-6,'GradObj','on' ,'DerivativeCheck','off','FinDiffType','central');

if Beta==0
    if exist(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'Beta0.mat'],'file')==2
        load(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'Beta0.mat'],'bstruc');
        startval = bstruc;
    else
        if nloc==3
            startval = .9*[0;.27;.44;-1;-.38;.08;-.8;-.06;.0006;.4;-.08;.036;-.22;.22];
        else
            startval = [rand(nloc,1);-1;2;1;-2;-.5;.0005;-4;-3;1.02;-.5;.0005];
        end
    end
else
    if exist(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'.mat'],'file')==2
        load(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'.mat'],'bstruc');
        startval = bstruc;
    end
end

derivative_checker = false;
if derivative_checker==true
	[bstruc0,lstruc,e,o,gNum]=fminunc('clogitAdj1',startval,o4Nu,restrMat,YS(flagSl),[],ZS(flagSl,:,:),baseAlt,[],Adj(flagSl,:),1);
	[bstruc0,lstruc,e,o,gAna]=fminunc('clogitAdj1',startval,o4An,restrMat,YS(flagSl),[],ZS(flagSl,:,:),baseAlt,[],Adj(flagSl,:),1);
	dlmwrite ('gradientChecker.csv',[gNum gAna gNum-gAna]);
	return
end
[bstruc,lstruc,~,~,~,hstruc] = fminunc('clogitAdj1',startval,options,restrMat,YS(flagSl),[],ZS(flagSl,:,:),baseAlt,PTypel(flagSl),Adj(flagSl,:),1);
[bstruc,invHstruc] = applyRestr(restrMat,bstruc,hstruc);
Pstruc = nan(N*T*S,J);
Pstruc(flagSl,:) = pclogitAdj1(bstruc,YS(flagSl),[],ZS(flagSl,:,:),baseAlt,Adj(flagSl,:),1);
if Beta==0
save(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'Beta0.mat'],'-v7.3','bstruc','Pstruc','YS','ZS','Adj','pimatS','flagSl');
dlmwrite(['estimatesFrict',money,'_',num2str(nloc),'loc',sample,'Beta0.csv'],[[bstruc sqrt(diag(invHstruc))];[-lstruc NaN];[sum(flagSl) NaN];[numel(unique(IDS(flagS))) NaN]]);
else
save(['strucFrictHet_',money,'_',num2str(nloc),'loc',sample,'.mat'],'-v7.3','bstruc','Pstruc','YS','ZS','Adj','pimatS','flagSl');
dlmwrite(['estimatesFrict',money,'_',num2str(nloc),'loc',sample,'.csv'],[[bstruc sqrt(diag(invHstruc))];[-lstruc NaN];[sum(flagSl) NaN];[numel(unique(IDS(flagS))) NaN]]);
end
[bstruc sqrt(diag(invHstruc))]
[[-lstruc NaN];[sum(flagSl) NaN];[numel(unique(IDS(flagS))) NaN]]

disp(['Estimation took ',num2str(toc/3600),' hours']);
diary off;
