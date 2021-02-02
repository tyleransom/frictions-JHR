function [Z,pimat] = getZcflWageAll(yr,ager,LEi,LYi,experer,birthLoc,birthDiv,lambda_hat,delta_hat,lambda_u_hat,lambda_e_hat,J,rho_hat_lambda,rho_hat_lambda_u,rho_hat_lambda_e,rho_hat_delta,rho_hat_wage,bwage,shock_amt,distance,Beta)

corrw  = corrcov(ARcov(1:2:length(ARcov),1:2:length(ARcov)));

nloc = J/2;
Tmater = zeros(1,10);
Tmater(yr-2003)=1;
l = LYi;
pimat=nan(1,nloc);
omega=nan(1,nloc);
for k=1:nloc
    if (LEi==0 || LYi>nloc) && k==l
        pimat(1,k) = lambda_hat(k,yr-2003);
    end
    if LEi==1 && k==l
        pimat(1,k) = 1-delta_hat(k,yr-2003);
    end
    if (LEi==0 || LYi>nloc) && k~=l
        pimat(1,k) = lambda_u_hat(k,yr-2003);
    end
    if LEi==1 && k~=l
        pimat(1,k) = lambda_e_hat(k,yr-2003);
    end
end
% no integration
pimat1(1,1) = piPred(lambda_hat(l,yr-2003),rho_hat_lambda,l); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
% integration
% pimat1(i,1) = piInt(lambda_hat(l,ti),rho_hat_lambda,sig_hat_lambda,l,10); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
for k=1:nloc
    omega(1,k) = pimat(1,k)./pimat1(1,1);
end

Zloc = zeros(1,nloc,J);
Ztilde0eA = zeros(1,13,J);
Ztilde1eA = zeros(1,13,J);
Ztilde2eA = zeros(1,13,J);
Ztilde0uA = zeros(1,13,J);
Ztilde1uA = zeros(1,13,J);
Ztilde2uA = zeros(1,13,J);
Ztilde0B  = zeros(1,13,J);
Ztilde1eB = zeros(1,13,J);
Ztilde2eB = zeros(1,13,J);
Ztilde1uB = zeros(1,13,J);
Ztilde2uB = zeros(1,13,J);
Ztilde1nB = zeros(1,13,J);
Ztilde2nB = zeros(1,13,J);
for j=1:nloc
    lp = j;
    k=j+nloc;
    Ztilde0eA(1,:,j) = repmat((  pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,j      ,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,j      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage) birthLoc(1,j      ) birthDiv(1,j      ) switchervec(nloc,LYi    ,j      ,ager) movervec(nloc,LYi    ,j      ,distance,ager)];
    Ztilde0eA(1,:,k) = repmat((  pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage) birthLoc(1,k      ) birthDiv(1,k      ) switchervec(nloc,LYi    ,k      ,ager) movervec(nloc,LYi    ,k      ,distance,ager)];
    Ztilde0uA(1,:,j) = repmat((1-pimat(1,j)),[1 13])               .*[1*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage) birthLoc(1,j      ) birthDiv(1,j      ) switchervec(nloc,LYi    ,j      ,ager) movervec(nloc,LYi    ,j      ,distance,ager)];
    Ztilde0uA(1,:,k) = repmat((1-pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage) birthLoc(1,k      ) birthDiv(1,k      ) switchervec(nloc,LYi    ,k      ,ager) movervec(nloc,LYi    ,k      ,distance,ager)];
    Ztilde1eA(1,:,j) = repmat((  pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer+1,bwage) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,j      ,lp+nloc,ager) movervec(nloc,j      ,lp+nloc,distance,ager)];
    Ztilde1eA(1,:,k) = repmat((  pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,k      ,lp+nloc,ager) movervec(nloc,k      ,lp+nloc,distance,ager)];
    Ztilde1uA(1,:,j) = repmat((1-pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,j      ,lp+nloc,ager) movervec(nloc,j      ,lp+nloc,distance,ager)];
    Ztilde1uA(1,:,k) = repmat((1-pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,k      ,lp+nloc,ager) movervec(nloc,k      ,lp+nloc,distance,ager)];
    Ztilde2eA(1,:,j) = repmat((  pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer+1,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,lp+nloc,l +nloc,ager) movervec(nloc,lp+nloc,l +nloc,distance,ager)];
    Ztilde2eA(1,:,k) = repmat((  pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,lp+nloc,l +nloc,ager) movervec(nloc,lp+nloc,l +nloc,distance,ager)];
    Ztilde2uA(1,:,j) = repmat((1-pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,lp+nloc,l +nloc,ager) movervec(nloc,lp+nloc,l +nloc,distance,ager)];
    Ztilde2uA(1,:,k) = repmat((1-pimat(1,j)),[1 13])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,lp+nloc,l +nloc,ager) movervec(nloc,lp+nloc,l +nloc,distance,ager)];

    Ztilde0B(1,:,j)  =                                               [0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,LYi    ,k      ,ager) movervec(nloc,LYi    ,l +nloc,distance,ager)];
    Ztilde0B(1,:,k)  =                                               [0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,LYi    ,k      ,ager) movervec(nloc,LYi    ,l +nloc,distance,ager)];
    Ztilde1eB(1,:,j) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l      ,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l      ,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l      ,ager) movervec(nloc,l +nloc,l      ,distance,ager)];
    Ztilde1eB(1,:,k) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    Ztilde1uB(1,:,j) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 13]).*[1*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l      ,ager) movervec(nloc,l +nloc,l      ,distance,ager)];
    Ztilde1uB(1,:,k) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    Ztilde1nB(1,:,j) = repmat((1-               omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    Ztilde1nB(1,:,k) = repmat((1-               omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    Ztilde2eB(1,:,j) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer+1,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l      ,l +nloc,ager) movervec(nloc,l      ,l +nloc,distance,ager)];
    Ztilde2eB(1,:,k) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    Ztilde2uB(1,:,j) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l      ,l +nloc,ager) movervec(nloc,l      ,l +nloc,distance,ager)];
    Ztilde2uB(1,:,k) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    Ztilde2nB(1,:,j) = repmat((1-               omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    Ztilde2nB(1,:,k) = repmat((1-               omega(1,j)),[1 13]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageAll(nloc,corrw(l,:),shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l +nloc,ager) movervec(nloc,l +nloc,l +nloc,distance,ager)];
    if j~=l
        Zloc(1,j,j) =  (1+Beta);
        Zloc(1,l,j) = -(1+Beta);
    end
    if k~=nloc+l
        Zloc(1,j,k) =  (1+Beta);
        Zloc(1,l,k) = -(1+Beta);
    end
    if k==nloc+l
        Zloc(1,:,[l l+nloc]) = zeros(1,nloc,2);
    end
end
Ztilde = Ztilde0eA+Ztilde0uA-Ztilde0B + Beta*(Ztilde1eA+Ztilde1uA-Ztilde1eB-Ztilde1uB-Ztilde1nB) + Beta.^2*(Ztilde2eA+Ztilde2uA-Ztilde2eB-Ztilde2uB-Ztilde2nB);
Z = cat(2,Zloc,Ztilde);

end
