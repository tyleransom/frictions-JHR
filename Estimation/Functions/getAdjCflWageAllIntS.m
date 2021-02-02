function [Adj] = getAdjCflWageAllIntS(b,draws,ARcov,ARcov2,yr,ager,LEi,LYi,experer,typer,birthLoc,birthDiv,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,shock_amt,distance,baseAlt,Beta)

corrw  = corrcov(ARcov(1:2:length(ARcov),1:2:length(ARcov)));

nloc = J/2;
shocker = nan(1,2*nloc,draws);
for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov);
end

Tmater = zeros(1,10);
Tmater(yr-2003)=1;
l = LYi;
pimat=nan(1,nloc);
pimat1w=nan(1,nloc,draws);
omegaw=nan(1,nloc,draws);
homeLocW= repmat(l,[1 1 draws]);

for k=1:nloc
    if LEi==0 && (k==LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = lambdervecRhoS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,1,experer,typer,rho_hat_urate  ,0,bLemp,bLunemp);
    elseif LEi==1 && (k==LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = 1-lambdervecRhoS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,1,experer,typer,rho_hat_urate   ,0,bLemp,bLunemp);
    elseif LEi==0 && (k~=LYi-nloc*(LYi>nloc));
       pimat(1,k)  = lambdervecRhoS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,0,experer,typer,rho_hat_urate,0,bLemp,bLunemp);
    elseif LEi==1 && (k~=LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = lambdervecRhoS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,0,experer,typer,rho_hat_urate,0,bLemp,bLunemp);
    end
end
pimat1w(1,1,:) = lambdervecRhoIntS(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp);
for jj=1:nloc 
    omegaw(1,jj,:) = repmat(pimat(1,jj),[1 1 draws])./pimat1w(1,1,:);
end

Adj = zeros(1,J);

% ========================================================================
% Construct CCP's from data and mlogit coefficients and Integrate!
% ==========================================================================
X = [ones(1,2)];
Ztemp = zeros(1,17,J,draws);
P = rand(1,J,draws);
% Need to loop over Z matrix in flexible logit to cover each set of states

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov);
end
% Loop 1a: Pr({0,lprime} t+1 | {jprime,lprime} t and employed in t)
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+1 alternatives
    for j=1:J
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockd1 = squeeze(shocker(1,(j-1)*2+2,:))';
            if j==lp
            pim = 1-lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,1,1,experer+jp,typer,rho_hat_urate,1,shockd1,bLemp,bLunemp);
            % size(pim)
            % size(ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,1,experer+jp,typer,bwage,shockw1))
            % size(repmat(movervecS(nloc,k,j,distance,ager+1,1,0),[1 1 draws]))
            Ztemp(1,:,j,:) = [pim zeros(1,1,draws) ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,1,experer+jp,typer,bwage,shockw1)  repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ager+1,1,0,typer),[1 1 draws])];
            else
            pim = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,1,0,experer+jp,typer,rho_hat_urate,1,shockd1,bLemp,bLunemp);
            Ztemp(1,:,j,:) = [pim zeros(1,1,draws) ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,1,experer+jp,typer,bwage,shockw1)  repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ager+1,1,0,typer),[1 1 draws])];
            end
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ager+1,1,0,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -Beta.*(pimat(:,k)).*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
    else
        Adj(:,k) = -Beta.*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
    end
end


for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov);
end
% Loop 1b: Pr({0,lprime} t+1 | {jprime,lprime} t and unemployed in t)
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+1 alternatives
    for j=1:J 
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockl1 = squeeze(shocker(1,(j-1)*2+2,:))';
            if j==lp
            pim = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,shockl1,bLemp,bLunemp);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,1,experer,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ager+1,0,1,typer),[1 1 draws])];
            else
            pim = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,typer,rho_hat_urate,1,shockl1,bLemp,bLunemp);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,1,experer,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ager+1,0,1,typer),[1 1 draws])];
            end
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,k,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,k,j,distance,ager+1,0,1,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -Beta.*(1-pimat(:,k)).*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
    end
end

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov2);
end
% Loop 2a: Pr({0,l} t+2 | {0,lprime} t+1, {jprime,lprime} t and employed in t)
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+2 alternatives, but now in the Z's, fix k=nloc+lp (the t+1 decision)
    for j=1:J 
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockl1 = squeeze(shocker(1,(j-1)*2+2,:))';
            entry1 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer+jp,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
    		entry2 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer+jp,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
			pim = entry1.*(homeLocW==j)+entry2.*(homeLocW~=j);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,2,experer+jp,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ager+2,0,0,typer),[1 1 draws])];
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ager+2,0,0,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -(Beta^2).*(pimat(:,k)).*((1/draws)*sum(log(P(:,l+nloc,:)),3)) + Adj(:,k);
    else
        Adj(:,k) = -(Beta^2).*((1/draws)*sum(log(P(:,l+nloc,:)),3)) + Adj(:,k);
    end
end

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov2);
end
% Loop 2b: Pr({0,l} t+2 | {0,lprime} t+1, {jprime,lprime} t and unemployed in t)
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+2 alternatives, but now in the Z's, fix k=nloc+lp (the t+1 decision)
    for j=1:J 
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockl1 = squeeze(shocker(1,(j-1)*2+2,:))';
            entry1 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
    		entry2 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
			pim = entry1.*(homeLocW==j)+entry2.*(homeLocW~=j);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,2,experer,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ager+2,0,0,typer),[1 1 draws])];
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+lp,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+lp,j,distance,ager+2,0,0,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -(Beta^2).*(1-pimat(:,k)).*((1/draws)*sum(log(P(:,l+nloc,:)),3)) + Adj(:,k);
    end
end

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov);
end
for jj=1:nloc
	pimat1w(1,1,:) = lambdervecRhoIntS(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
    omegaw(1,jj,:) = repmat(pimat(1,jj),[1 1 draws])./pimat1w(1,1,:);
end
% Loop 3a: Pr({jprime,l} t+1 | {0,l} t) -- weight by omega(:,k) [employment counterpath]
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+1 alternatives, fixing k=nloc+1 (the t decision)
    for j=1:J 
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockl1 = squeeze(shocker(1,(j-1)*2+2,:))';
        	entry1 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,shockl1,bLemp,bLunemp);
        	entry2 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,typer,rho_hat_urate,1,shockl1,bLemp,bLunemp);
        	pim = entry1.*(homeLocW==j)+entry2.*(homeLocW~=j);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,1,experer,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+l,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+l,j,distance,ager+1,0,0,typer),[1 1 draws])];
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+l,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+l,j,distance,ager+1,0,0,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = Beta.*((1/draws)*sum(omegaw(:,k,:).*log(P(:,l,:)),3)) + Adj(:,k);
    else
        Adj(:,k) = Beta.*((1/draws)*sum(log(P(:,l+nloc,:)),3)) + Adj(:,k);
    end
end

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov);
end
for jj=1:nloc
	pimat1w(1,1,:) = lambdervecRhoIntS(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
    omegaw(1,jj,:) = repmat(pimat(1,jj),[1 1 draws])./pimat1w(1,1,:);
end
% Loop 3b: Pr({jprime,l} t+1 | {0,l} t) -- weight by 1-omega(:,k) [unemployment counterpath]
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+1 alternatives, fixing k=nloc+1 (the t decision)
    for j=1:J 
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockl1 = squeeze(shocker(1,(j-1)*2+2,:))';
        	entry1 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,shockl1,bLemp,bLunemp);
        	entry2 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,typer,rho_hat_urate,1,shockl1,bLemp,bLunemp);
			pim = entry1.*(homeLocW==j)+entry2.*(homeLocW~=j);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,1,experer,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+l,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+l,j,distance,ager+1,0,0,typer),[1 1 draws])];
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+l,j,ager+1,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+l,j,distance,ager+1,0,0,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = Beta.*((1/draws)*sum((1-omegaw(:,k,:)).*log(P(:,l+nloc,:)),3)) + Adj(:,k);
    end
end

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov2);
end
for jj=1:nloc
	pimat1w(1,1,:) = lambdervecRhoIntS(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
    omegaw(1,jj,:) = repmat(pimat(1,jj),[1 1 draws])./pimat1w(1,1,:);
end
% Loop 4a: Pr({0,l} t+2 | {jprime,l} t+1, {0,l} t) [employment counterpath]
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
    for j=1:J
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockd1 = squeeze(shocker(1,(j-1)*2+2,:))';
        	entry1 = 1-lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,1,1,experer+jp,typer,rho_hat_urate,2,shockd1,bLemp,bLunemp);
        	entry2 =   lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,1,0,experer+jp,typer,rho_hat_urate,2,shockd1,bLemp,bLunemp);
        	pim = entry1.*(homeLocW==j)+entry2.*(homeLocW~=j);
            Ztemp(1,:,j,:) = [pim zeros(1,1,draws) ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,2,experer+jp,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+l,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+l,j,distance,ager+2,1,0,typer),[1 1 draws])];
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+l,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+l,j,distance,ager+2,1,0,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = (Beta^2).*((1/draws)*sum(pimat1w(:,1,:).*omegaw(:,k,:).*log(P(:,l+nloc,:)),3)) + Adj(:,k);
    else
        Adj(:,k) = (Beta^2).*((1/draws)*sum(log(P(:,l+nloc,:)),3)) + Adj(:,k);
    end
end

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov2);
end
for jj=1:nloc
	pimat1w(1,1,:) = lambdervecRhoIntS(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
    omegaw(1,jj,:) = repmat(pimat(1,jj),[1 1 draws])./pimat1w(1,1,:);
end
% Loop 4b: Pr({0,l} t+2 | {jprime,l} t+1, {0,l} t) [unemployment counterpath]
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
    for j=1:J 
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockd1 = squeeze(shocker(1,(j-1)*2+2,:))';
        	entry1 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
        	entry2 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
        	pim = entry1.*(homeLocW==j)+entry2.*(homeLocW~=j);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,2,experer,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+l,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+l,j,distance,ager+2,0,1,typer),[1 1 draws])];
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc*(k>nloc)+l,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc*(k>nloc)+l,j,distance,ager+2,0,1,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = (Beta^2).*((1/draws)*sum((1-pimat1w(:,1,:)).*omegaw(:,k,:).*log(P(:,l+nloc,:)),3)) + Adj(:,k);
    end
end

for dd=1:draws
    shocker(1,:,dd) = mvnrnd(zeros(1,2*nloc),ARcov2);
end
for jj=1:nloc
	pimat1w(1,1,:) = lambdervecRhoIntS(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
    omegaw(1,jj,:) = repmat(pimat(1,jj),[1 1 draws])./pimat1w(1,1,:);
end
% Loop 4c: Pr({0,l} t+2 | {0,l} t+1, {0,l} t) [NILF counterpath]
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
    for j=1:J 
        if j<=nloc
            shockw1 = squeeze(shocker(1,(j-1)*2+1,:))';
            shockd1 = squeeze(shocker(1,(j-1)*2+2,:))';
        	entry1 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
        	entry2 = lambdervecRhoIntS(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,typer,rho_hat_urate,2,shockl1,bLemp,bLunemp);
        	pim = entry1.*(homeLocW==j)+entry2.*(homeLocW~=j);
            Ztemp(1,:,j,:) = [zeros(1,1,draws) pim ewagevecRhoCflWageAllIntS(nloc,corrw(l,:),shock_amt,j,Tmater,yr-2003,rho_hat_wage,2,experer,typer,bwage,shockw1) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+l,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+l,j,distance,ager+2,0,0,typer),[1 1 draws])];
        else
            Ztemp(1,:,j,:) = [zeros(1,3,draws) repmat(birthLoc(1,j),[1 1 draws]) repmat(birthDiv(1,j),[1 1 draws]) repmat(switchervecS(nloc,nloc+l,j,ager+2,typer),[1 1 draws]) repmat(movervecS(nloc,nloc+l,j,distance,ager+2,0,0,typer),[1 1 draws])];
        end
    end
    % calculate CCPs
    for dd=1:draws
        P(1,:,dd) = pclogit1(b,l,X,Ztemp(:,:,:,dd),size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    end
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = (Beta^2).*((1/draws)*sum((1-omegaw(:,k,:)).*log(P(:,l+nloc,:)),3)) + Adj(:,k);
    end
end

end
