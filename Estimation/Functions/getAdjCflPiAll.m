function [Adj] = getAdjCflPiAll(b,yr,ager,LEi,LYi,experer,birthLoc,birthDiv,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,shock_amt,distance,baseAlt,Beta)

corru  = corrcov(ARcov(2:2:length(ARcov),2:2:length(ARcov)));

nloc = J/2;
Tmater = zeros(1,10);
Tmater(yr-2003)=1;
l = LYi;
pimat=nan(1,nloc);
omega=nan(1,nloc);

for k=1:nloc
    if LEi==0 && (k==LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = lambdervecRhoCflPiAll(nloc,k,Tmater,yr-2003,urate_lag55,LEi==0,1,experer,rho_hat_urate  ,0,bLemp,bLunemp,corru(l,k),shock_amt);
    elseif LEi==1 && (k==LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = 1-lambdervecRhoCflPiAll(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,1,experer,rho_hat_urate   ,0,bLemp,bLunemp,corru(l,k),shock_amt);
    elseif LEi==0 && (k~=LYi-nloc*(LYi>nloc));
       pimat(1,k)  = lambdervecRhoCflPiAll(nloc,k,Tmater,yr-2003,urate_lag55,LEi==0,0,experer,rho_hat_urate,0,bLemp,bLunemp,corru(l,k),shock_amt);
    elseif LEi==1 && (k~=LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = lambdervecRhoCflPiAll(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,0,experer,rho_hat_urate,0,bLemp,bLunemp,corru(l,k),shock_amt);
    end
    % no integration
    if (k==LYi-nloc*(LYi>nloc)); 
       pimat1(1,1) = lambdervecRhoCflPiAll(nloc,k,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,k),shock_amt);
    end
    omega(1,k) = pimat(1,k)./pimat1(1,1);
end

Adj = zeros(1,J);

X = [ones(1,2)];
Ztemp = zeros(1,13,J);
P = rand(1,J);
% Need to loop over Z matrix in flexible logit to cover each set of states
% Loop 1a: Pr({0,lprime} t+1 | {jprime,lprime} t and employed in t)
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+1 alternatives
    for j=1:J
        if j<=nloc
            if j==lp
            Ztemp(1,:,j) = [1-lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,1,1,experer+jp,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer+jp,bwage) (lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,1,1,experer+jp,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer+jp,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,k,j,ager+1) movervec(nloc,k,j,distance,ager+1)];
            else
            Ztemp(1,:,j) = [lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,1,0,experer+jp,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer+jp,bwage) lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,1,0,experer+jp,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer+jp,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,k,j,ager+1) movervec(nloc,k,j,distance,ager+1)];
            end
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,k,j,ager+1) movervec(nloc,k,j,distance,ager+1)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -Beta.*(pimat(:,k)).*log(P(:,nloc+lp)) + Adj(:,k);
    else
        Adj(:,k) = -Beta.*log(P(:,nloc+lp)) + Adj(:,k);
    end
end

% Loop 1b: Pr({0,lprime} t+1 | {jprime,lprime} t and unemployed in t)
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+1 alternatives
    for j=1:J 
        if j<=nloc
            if j==lp
            Ztemp(1,:,j) = [lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,k,j,ager+1) movervec(nloc,k,j,distance,ager+1)];
            else
            Ztemp(1,:,j) = [lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,k,j,ager+1) movervec(nloc,k,j,distance,ager+1)];
            end
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,k,j,ager+1) movervec(nloc,k,j,distance,ager+1)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -Beta.*(1-pimat(:,k)).*log(P(:,nloc+lp)) + Adj(:,k);
    end
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
            entry1 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer+jp,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            entry2 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer+jp,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            Ztemp(1,:,j) = [entry1.*(l==j)+entry2.*(l~=j) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer+jp,bwage) (entry1.*(l==j)+entry2.*(l~=j)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer+jp,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+lp,j,ager+2) movervec(nloc,nloc+lp,j,distance,ager+2)];
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+lp,j,ager+2) movervec(nloc,nloc+lp,j,distance,ager+2)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -(Beta^2).*(pimat(:,k)).*log(P(:,l)) + Adj(:,k);
    else
        Adj(:,k) = -(Beta^2).*log(P(:,l+nloc)) + Adj(:,k);
    end
end

% Loop 2b: Pr({0,l} t+2 | {0,lprime} t+1, {jprime,lprime} t and unemployed in t)
% loop over all alternatives to store in FV matrix
for k=1:J
    jp = k<=nloc;
    lp = k-nloc*(k>nloc);
    % loop over all t+2 alternatives, but now in the Z's, fix k=nloc+lp (the t+1 decision)
    for j=1:J 
        if j<=nloc
            entry1 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            entry2 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            Ztemp(1,:,j) = [entry1.*(l==j)+entry2.*(l~=j) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer,bwage) (entry1.*(l==j)+entry2.*(l~=j)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+lp,j,ager+2) movervec(nloc,nloc+lp,j,distance,ager+2)];
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+lp,j,ager+2) movervec(nloc,nloc+lp,j,distance,ager+2)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = -(Beta^2).*(1-pimat(:,k)).*log(P(:,l+nloc)) + Adj(:,k);
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
            entry1 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt);
            entry2 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt);
            Ztemp(1,:,j) = [entry1.*(l==j)+entry2.*(l~=j) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) (entry1.*(l==j)+entry2.*(l~=j)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+l,j,ager+1) movervec(nloc,nloc+l,j,distance,ager+1)];
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+l,j,ager+1) movervec(nloc,nloc+l,j,distance,ager+1)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = Beta.*(omega(:,k)).*log(P(:,l)) + Adj(:,k);
    else
        Adj(:,k) = Beta.*log(P(:,l+nloc)) + Adj(:,k);
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
            entry1 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt);
            entry2 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,rho_hat_urate,1,bLemp,bLunemp,corru(l,j),shock_amt);
            Ztemp(1,:,j) = [entry1.*(l==j)+entry2.*(l~=j) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) (entry1.*(l==j)+entry2.*(l~=j)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,1,experer,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+l,j,ager+1) movervec(nloc,nloc+l,j,distance,ager+1)];
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+l,j,ager+1) movervec(nloc,nloc+l,j,distance,ager+1)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = Beta.*(1-omega(:,k)).*log(P(:,l)) + Adj(:,k);
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
            entry1 = 1-lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,1,1,experer+jp,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            entry2 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,1,0,experer+jp,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            Ztemp(1,:,j) = [entry1.*(l==j)+entry2.*(l~=j) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer+jp,bwage) (entry1.*(l==j)+entry2.*(l~=j)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer+jp,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc*(k>nloc)+l,j,ager+2) movervec(nloc,nloc*(k>nloc)+l,j,distance,ager+2)];
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc*(k>nloc)+l,j,ager+2) movervec(nloc,nloc*(k>nloc)+l,j,distance,ager+2)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = (Beta^2).*(pimat1(:,1).*omega(:,k)).*log(P(:,l+nloc)) + Adj(:,k);
    else
        Adj(:,k) = (Beta^2).*log(P(:,l+nloc)) + Adj(:,k);
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
            entry1 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            entry2 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            Ztemp(1,:,j) = [entry1.*(l==j)+entry2.*(l~=j) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer,bwage) (entry1.*(l==j)+entry2.*(l~=j)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc*(k>nloc)+l,j,ager+2) movervec(nloc,nloc*(k>nloc)+l,j,distance,ager+2)];
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc*(k>nloc)+l,j,ager+2) movervec(nloc,nloc*(k>nloc)+l,j,distance,ager+2)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = (Beta^2).*((1-pimat1(:,1)).*omega(:,k)).*log(P(:,l+nloc)) + Adj(:,k);
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
            entry1 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            entry2 = lambdervecRhoCflPiAll(nloc,j,Tmater,yr-2003,urate_lag55,0,0,experer,rho_hat_urate,2,bLemp,bLunemp,corru(l,j),shock_amt);
            Ztemp(1,:,j) = [entry1.*(l==j)+entry2.*(l~=j) ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer,bwage) (entry1.*(l==j)+entry2.*(l~=j)).*ewagevecRho(nloc,j,Tmater,yr-2003,rho_hat_wage,2,experer,bwage) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+l,j,ager+2) movervec(nloc,nloc+l,j,distance,ager+2)];
        else
            Ztemp(1,:,j) = [zeros(1,3) birthLoc(1,j) birthDiv(1,j) switchervec(nloc,nloc+l,j,ager+2) movervec(nloc,nloc+l,j,distance,ager+2)];
        end
    end
    % calculate CCPs
    P = pclogit1(b,l,X,Ztemp,size(X,2),size(Ztemp,2),size(Ztemp,3),baseAlt);
    % store log CCPs in the Adjustment term for the period t alternative (k)
    if k<=nloc
        Adj(:,k) = (Beta^2).*(1-omega(:,k)).*log(P(:,l+nloc)) + Adj(:,k);
    end
end

end
