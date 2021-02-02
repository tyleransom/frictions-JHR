function [Adj] = estAdjFun(b,LY,baseAlt,N,T,draws,flag,flagw,empFT_lag,ARcov,ARcov2,J,nloc,Tmat,Tvec,urate_lag55,exper,rho_hat_urate,rho_hat_wage,bLemp,bLunemp,bwage,birthLoc,birthDiv,age,distance,Beta)

	Lh      = zeros(N*T,J);
	Lw      = zeros(N*T,J);
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
	homeLocW= repmat(homeLoc,[1 1 draws]);

	for k=1:nloc
		if flag && empFT_lag(:)==0 && (k==LY(:)-nloc*(LY(:)>nloc))
		pimat(1,k) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,1,exper,rho_hat_urate  ,0,bLemp,bLunemp);

        elseif flag && empFT_lag(:)==1 && (k==LY(:)-nloc*(LY(:)>nloc));
		pimat(1,k) = 1-lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,1,exper,rho_hat_urate   ,0,bLemp,bLunemp);

        elseif flag && empFT_lag(:)==0 && (k~=LY(:)-nloc*(LY(:)>nloc));
		pimat(1,k) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,0,exper,rho_hat_urate,0,bLemp,bLunemp);

        elseif flag && empFT_lag(:)==1 && (k~=LY(:)-nloc*(LY(:)>nloc));
		pimat(1,k) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,0,exper,rho_hat_urate,0,bLemp,bLunemp);
		
        elseif flag && (k==LY(:)-nloc*(LY(:)>nloc));
		pimat1(1,1) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,0,1,exper,rho_hat_urate,1,bLemp,bLunemp);
        end
	end
	
	X = [ones(sum(flag),2)];
	Ztemp = zeros(sum(flag),15,J,draws);
	P = rand(sum(flag),J,draws);

    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov);
    end
    % Loop 1a: Pr({0,lprime} t+1 | {jprime,lprime} t and employed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives
        for j=1:J
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockd1 = squeeze(shocker(flag,(j-1)*2+2,:));
                if j==lp
                pim = 1-lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,1,1,exper(flagw)+jp,rho_hat_urate,1,shockd1,bLemp,bLunemp);
                size(pim);
                Ztemp(flag,:,j,:) = [pim zeros(sum(flag),1,draws) ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)+jp,bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,k,j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,k,j,distance,age(flagw)+1,ones(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
                else
                pim = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,1,0,exper(flagw)+jp,rho_hat_urate,1,shockd1,bLemp,bLunemp);
                Ztemp(flag,:,j,:) = [pim zeros(sum(flag),1,draws) ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)+jp,bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,k,j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,k,j,distance,age(flagw)+1,ones(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
                end
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,k,j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,k,j,distance,age(flagw)+1,ones(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -Beta.*(pimat(:,k)).*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
        else
            Adj(:,k) = -Beta.*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
        end
    end

    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov);
    end
    % Loop 1b: Pr({0,lprime} t+1 | {jprime,lprime} t and unemployed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flag,(j-1)*2+2,:));
                if j==lp
                pim = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,1,exper(flagw),rho_hat_urate,1,shockl1,bLemp,bLunemp);
                Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw),bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,k,j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,k,j,distance,age(flagw)+1,zeros(sum(flag),1),ones(sum(flag),1)),[1 1 draws])];
                else
                pim = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,0,exper(flagw),rho_hat_urate,1,shockl1,bLemp,bLunemp);
                Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw),bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,k,j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,k,j,distance,age(flagw)+1,zeros(sum(flag),1),ones(sum(flag),1)),[1 1 draws])];
                end
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,k,j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,k,j,distance,age(flagw)+1,zeros(sum(flag),1),ones(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -Beta.*(1-pimat(:,k)).*((1/draws)*sum(log(P(:,nloc+lp,:)),3)) + Adj(:,k);
        end
    end

    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov2);
    end
    % Loop 2a: Pr({0,l} t+2 | {0,lprime} t+1, {jprime,lprime} t and employed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=nloc+lp (the t+1 decision)
        % if lp==l, then lambda
        % otherwise, lambda_u
        for j=1:J
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flag,(j-1)*2+2,:));
                entry1 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,1,exper(flagw)+jp,rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,0,exper(flagw)+jp,rho_hat_urate,2,shockl1,bLemp,bLunemp);
				pim = entry1.*(homeLocW(flag,:,:)==j)+entry2.*(homeLocW(flag,:,:)~=j);				
                Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)+jp,bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+lp,j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc+lp,j,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+lp,j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc+lp,j,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -(Beta^2).*(pimat(:,k)).*((1/draws)*sum(log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        else
            Adj(:,k) = -(Beta^2).*((1/draws)*sum(log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        end
    end

    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov2);
    end
    % Loop 2b: Pr({0,l} t+2 | {0,lprime} t+1, {jprime,lprime} t and unemployed in t)
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=nloc+lp (the t+1 decision)
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flag,(j-1)*2+2,:));
                entry1 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,1,exper(flagw),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,0,exper(flagw),rho_hat_urate,2,shockl1,bLemp,bLunemp);
				pim = entry1.*(homeLocW(flag,:,:)==j)+entry2.*(homeLocW(flag,:,:)~=j);
                Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw),bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+lp,j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc+lp,j,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+lp,j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc+lp,j,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = -(Beta^2).*(1-pimat(:,k)).*((1/draws)*sum(log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        end
    end
    
    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov);
    end
    for i=1:N*T %:numel(exper) % loop over individual-time observations
        if flag(i)==1
            LYi = LY(i);
            ti  = calyr(i)-2003;
            l = LYi-nloc*(LYi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
            for jj=1:nloc
                pimat1w(i,1,:) = lambdervecRhoInt(nloc,l,Tmat(i,:),Tvec(i),urate_lag55,0,1,exper(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
                omegaw(i,jj,:) = repmat(pimat(i,jj),[1 1 draws])./pimat1w(i,1,:);
            end
        end
    end
    % Loop 3a: Pr({jprime,l} t+1 | {0,l} t) -- weight by omega(:,k) [employment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives, fixing k=nloc+1 (the t decision)
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flag,(j-1)*2+2,:));
        		entry1 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,1,exper(flagw),rho_hat_urate,1,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,0,exper(flagw),rho_hat_urate,1,shockl1,bLemp,bLunemp);
        		pim = entry1.*(homeLocW(flag,:,:)==j)+entry2.*(homeLocW(flag,:,:)~=j);
                Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw),bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+homeLoc(flag),j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,nloc+homeLoc(flag),j,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+homeLoc(flag),j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,nloc+homeLoc(flag),j,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = Beta.*((1/draws)*sum(omegaw(:,k,:).*log(max(P.*Lw,[],2)),3)) + Adj(:,k);
        else
            Adj(:,k) = Beta.*((1/draws)*sum(log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        end
    end

    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov);
    end
    for i=1:N*T %:numel(exper) % loop over individual-time observations
        if flag(i)==1
            LYi = LY(i);
            ti  = calyr(i)-2003;
            l = LYi-nloc*(LYi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
            for jj=1:nloc
                pimat1w(i,1,:) = lambdervecRhoInt(nloc,l,Tmat(i,:),Tvec(i),urate_lag55,0,1,exper(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
                omegaw(i,jj,:) = repmat(pimat(i,jj),[1 1 draws])./pimat1w(i,1,:);
            end
        end
    end
    % Loop 3b: Pr({jprime,l} t+1 | {0,l} t) -- weight by 1-omega(:,k) [unemployment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+1 alternatives, fixing k=nloc+1 (the t decision)
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flag,(j-1)*2+2,:));
        		entry1 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,1,exper(flagw),rho_hat_urate,1,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,0,exper(flagw),rho_hat_urate,1,shockl1,bLemp,bLunemp);
				pim = entry1.*(homeLocW(flag,:,:)==j)+entry2.*(homeLocW(flag,:,:)~=j);
                Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw),bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+homeLoc(flag),j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,nloc+homeLoc(flag),j,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+homeLoc(flag),j,age(flagw)+1),[1 1 draws]) repmat(movervec(nloc,nloc+homeLoc(flag),j,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = Beta.*((1/draws)*sum((1-omegaw(:,k,:)).*log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        end
    end

    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov2);
    end
    for i=1:N*T %:numel(exper) % loop over individual-time observations
        if flag(i)==1
            LYi = LY(i);
            ti  = calyr(i)-2003;
            l = LYi-nloc*(LYi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
            for jj=1:nloc
                pimat1w(i,1,:) = lambdervecRhoInt(nloc,l,Tmat(i,:),Tvec(i),urate_lag55,0,1,exper(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
                omegaw(i,jj,:) = repmat(pimat(i,jj),[1 1 draws])./pimat1w(i,1,:);
            end
        end
    end
    % Loop 4a: Pr({0,l} t+2 | {jprime,l} t+1, {0,l} t) [employment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockd1 = squeeze(shocker(flag,(j-1)*2+2,:));
        		entry1 = 1-lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,1,1,exper(flagw)+jp,rho_hat_urate,2,shockd1,bLemp,bLunemp);
        		entry2 =   lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,1,0,exper(flagw)+jp,rho_hat_urate,2,shockd1,bLemp,bLunemp);
        		pim = entry1.*(homeLocW(flag,:,:)==j)+entry2.*(homeLocW(flag,:,:)~=j);
                Ztemp(flag,:,j,:) = [pim zeros(sum(flag),1,draws) ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)+jp,bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,distance,age(flagw)+2,ones(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,distance,age(flagw)+2,ones(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = (Beta^2).*((1/draws)*sum(pimat1w(:,1,:).*omegaw(:,k,:).*log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        else
            Adj(:,k) = (Beta^2).*((1/draws)*sum(log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        end
    end
    
    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov2);
    end
    for i=1:N*T %:numel(exper) % loop over individual-time observations
        if flag(i)==1
            LYi = LY(i);
            ti  = calyr(i)-2003;
            l = LYi-nloc*(LYi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
            for jj=1:nloc
                pimat1w(i,1,:) = lambdervecRhoInt(nloc,l,Tmat(i,:),Tvec(i),urate_lag55,0,1,exper(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
                omegaw(i,jj,:) = repmat(pimat(i,jj),[1 1 draws])./pimat1w(i,1,:);
            end
        end
    end
    % Loop 4b: Pr({0,l} t+2 | {jprime,l} t+1, {0,l} t) [unemployment counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flag,(j-1)*2+2,:));
        		entry1 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,1,exper(flagw),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,0,exper(flagw),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		pim = entry1.*(homeLocW(flag,:,:)==j)+entry2.*(homeLocW(flag,:,:)~=j);
                Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw),bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,distance,age(flagw)+2,zeros(sum(flag),1),ones(sum(flag),1)),[1 1 draws])];
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc*(k>nloc)+homeLoc(flag),j,distance,age(flagw)+2,zeros(sum(flag),1),ones(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = (Beta^2).*((1/draws)*sum((1-pimat1w(:,1,:)).*omegaw(:,k,:).*log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        end
    end

    
    for dd=1:draws
        shocker(flag,:,dd) = mvnrnd(zeros(sum(flag),2*nloc),ARcov2);
    end
    for i=1:N*T %:numel(exper) % loop over individual-time observations
        if flag(i)==1
            LYi = LY(i);
            ti  = calyr(i)-2003;
            l = LYi-nloc*(LYi>nloc);
            % Create pi_ikt and pi_it+1 for all locations
            for jj=1:nloc
                pimat1w(i,1,:) = lambdervecRhoInt(nloc,l,Tmat(i,:),Tvec(i),urate_lag55,0,1,exper(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
                omegaw(i,jj,:) = repmat(pimat(i,jj),[1 1 draws])./pimat1w(i,1,:);
            end
        end
    end
    % Loop 4c: Pr({0,l} t+2 | {0,l} t+1, {0,l} t) [NILF counterpath]
    % loop over all alternatives to store in FV matrix
    for k=1:J
        jp = k<=nloc;
        lp = k-nloc*(k>nloc);
        % loop over all t+2 alternatives, but now in the Z's, fix k=l or nloc+l (the t+1 decision)
        for j=1:J 
            if j<=nloc
                shockw1 = squeeze(shocker(flag,(j-1)*2+1,:));
                shockl1 = squeeze(shocker(flag,(j-1)*2+2,:));
        		entry1 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,1,exper(flagw),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		entry2 = lambdervecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),urate_lag55,0,0,exper(flagw),rho_hat_urate,2,shockl1,bLemp,bLunemp);
        		pim = entry1.*(homeLocW(flag,:,:)==j)+entry2.*(homeLocW(flag,:,:)~=j);
				Ztemp(flag,:,j,:) = [zeros(sum(flag),1,draws) pim ewagevecRhoInt(nloc,j,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw),bwage,shockw1) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+homeLoc(flag),j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc+homeLoc(flag),j,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            else
                Ztemp(flag,:,j,:) = [zeros(sum(flag),3,draws) repmat(birthLoc(flag,j),[1 1 draws]) repmat(birthDiv(flag,j),[1 1 draws]) repmat(switchervec(nloc,nloc+homeLoc(flag),j,age(flagw)+2),[1 1 draws]) repmat(movervec(nloc,nloc+homeLoc(flag),j,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1)),[1 1 draws])];
            end
        end
        % calculate CCPs
        for dd=1:draws
            P(flag,:,dd) = pclogit1(b,l,X(flag,:),Ztemp(flag,:,:,dd),baseAlt);
        end
        % store log CCPs in the Adjustment term for the period t alternative (k)
        if k<=nloc
            Adj(:,k) = (Beta^2).*((1/draws)*sum((1-omegaw(:,k,:)).*log(max(P.*Lh,[],2)),3)) + Adj(:,k);
        end
    end
    
end
