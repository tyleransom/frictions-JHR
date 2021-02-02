function [Z,pimat] = getZcflWageInt(draws,ARcov,ARcov2,yr,ager,LEi,LYi,experer,birthLoc,birthDiv,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,shock_amt,distance,Beta)



nloc = J/2;
Tmater = zeros(1,10);
Tmater(yr-2003)=1;
l = LYi;
pvEmp    = LYi<nloc & LEi==1;
pvUnemp  = LYi<nloc & LEi==0;
pimat=nan(1,nloc);
pimat1=nan(1,nloc);
omega=nan(1,nloc);

for k=1:nloc
    if LEi==0 && (k==LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = lambdervecRho(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,1,experer,rho_hat_urate  ,0,bLemp,bLunemp);
    elseif LEi==1 && (k==LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = 1-lambdervecRho(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,1,experer,rho_hat_urate   ,0,bLemp,bLunemp);
    elseif LEi==0 && (k~=LYi-nloc*(LYi>nloc));
       pimat(1,k)  = lambdervecRho(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,0,experer,rho_hat_urate,0,bLemp,bLunemp);
    elseif LEi==1 && (k~=LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = lambdervecRho(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,0,experer,rho_hat_urate,0,bLemp,bLunemp);
    end
end

Zloc = zeros(1,nloc,J);
Ztilde    = zeros(1,15,J);
Ztilde0eA = zeros(1,15,J);
Ztilde1eA = zeros(1,15,J);
Ztilde2eA = zeros(1,15,J);
Ztilde0uA = zeros(1,15,J);
Ztilde1uA = zeros(1,15,J);
Ztilde2uA = zeros(1,15,J);
Ztilde0B  = zeros(1,15,J);
Ztilde1eB = zeros(1,15,J);
Ztilde2eB = zeros(1,15,J);
Ztilde1uB = zeros(1,15,J);
Ztilde2uB = zeros(1,15,J);
Ztilde1nB = zeros(1,15,J);
Ztilde2nB = zeros(1,15,J);
shocker   = nan(1,2*nloc);
for dd=1:draws
    shocker(1,:,1) = mvnrnd(zeros(1,2*nloc),ARcov);
	for jj=1:nloc
        pimat1(1,1) = squeeze(lambdervecRhoInt(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp));
        omega(1,jj) = pimat(1,jj)./pimat1(1,1);
	end
	for j=1:nloc
        lp = j;
        k=j+nloc;
        Ztilde0eA(1,:,j) = repmat((  pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,j      ,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,j      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,j      ) birthDiv(1,j      ) switchervec(nloc,LYi    ,j      ,ager+0) movervec(nloc,LYi    ,j      ,distance,ager+0,pvEmp,pvUnemp)];
        Ztilde0eA(1,:,k) = repmat((  pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,k      ) birthDiv(1,k      ) switchervec(nloc,LYi    ,k      ,ager+0) movervec(nloc,LYi    ,k      ,distance,ager+0,pvEmp,pvUnemp)];
        Ztilde0uA(1,:,j) = repmat((1-pimat(1,j)),[1 15])               .*[1*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,j      ) birthDiv(1,j      ) switchervec(nloc,LYi    ,j      ,ager+0) movervec(nloc,LYi    ,j      ,distance,ager+0,pvEmp,pvUnemp)];
        Ztilde0uA(1,:,k) = repmat((1-pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,k      ) birthDiv(1,k      ) switchervec(nloc,LYi    ,k      ,ager+0) movervec(nloc,LYi    ,k      ,distance,ager+0,pvEmp,pvUnemp)];
        Ztilde1eA(1,:,j) = repmat((  pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer+1,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,j      ,lp+nloc,ager+1) movervec(nloc,j      ,lp+nloc,distance,ager+1,1,0)];
        Ztilde1eA(1,:,k) = repmat((  pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,k      ,lp+nloc,ager+1) movervec(nloc,k      ,lp+nloc,distance,ager+1,1,0)];
        Ztilde1uA(1,:,j) = repmat((1-pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,j      ,lp+nloc,ager+1) movervec(nloc,j      ,lp+nloc,distance,ager+1,0,1)];
        Ztilde1uA(1,:,k) = repmat((1-pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervec(nloc,k      ,lp+nloc,ager+1) movervec(nloc,k      ,lp+nloc,distance,ager+1,0,1)];
        Ztilde2eA(1,:,j) = repmat((  pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer+1,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,lp+nloc,l +nloc,ager+2) movervec(nloc,lp+nloc,l +nloc,distance,ager+2,0,0)];
        Ztilde2eA(1,:,k) = repmat((  pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,lp+nloc,l +nloc,ager+2) movervec(nloc,lp+nloc,l +nloc,distance,ager+2,0,0)];
        Ztilde2uA(1,:,j) = repmat((1-pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,lp+nloc,l +nloc,ager+2) movervec(nloc,lp+nloc,l +nloc,distance,ager+2,0,0)];
        Ztilde2uA(1,:,k) = repmat((1-pimat(1,j)),[1 15])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,lp+nloc,l +nloc,ager+2) movervec(nloc,lp+nloc,l +nloc,distance,ager+2,0,0)];

        Ztilde0B(1,:,j)  =                                               [0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,LYi    ,k      ,ager+0) movervec(nloc,LYi    ,l +nloc,distance,ager+0,pvEmp,pvUnemp)];
        Ztilde0B(1,:,k)  =                                               [0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,0,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,LYi    ,k      ,ager+0) movervec(nloc,LYi    ,l +nloc,distance,ager+0,pvEmp,pvUnemp)];
        Ztilde1eB(1,:,j) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l      ,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l      ,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l      ) birthDiv(1,l      ) switchervec(nloc,l +nloc,l      ,ager+1) movervec(nloc,l +nloc,l      ,distance,ager+1,0,0)];
        Ztilde1eB(1,:,k) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+1) movervec(nloc,l +nloc,l +nloc,distance,ager+1,0,0)];
        Ztilde1uB(1,:,j) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 15]).*[1*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l      ,ager+1) movervec(nloc,l +nloc,l      ,distance,ager+1,0,0)];
        Ztilde1uB(1,:,k) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+1) movervec(nloc,l +nloc,l +nloc,distance,ager+1,0,0)];
        Ztilde1nB(1,:,j) = repmat((1-               omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+1) movervec(nloc,l +nloc,l +nloc,distance,ager+1,0,0)];
        Ztilde1nB(1,:,k) = repmat((1-               omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+1) movervec(nloc,l +nloc,l +nloc,distance,ager+1,0,0)];
        Ztilde2eB(1,:,j) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer+1,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l      ,l +nloc,ager+2) movervec(nloc,l      ,l +nloc,distance,ager+2,1,0)];
        Ztilde2eB(1,:,k) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+2) movervec(nloc,l +nloc,l +nloc,distance,ager+2,1,0)];
        Ztilde2uB(1,:,j) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l      ,l +nloc,ager+2) movervec(nloc,l      ,l +nloc,distance,ager+2,0,1)];
        Ztilde2uB(1,:,k) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+2) movervec(nloc,l +nloc,l +nloc,distance,ager+2,0,1)];
        Ztilde2nB(1,:,j) = repmat((1-               omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+2) movervec(nloc,l +nloc,l +nloc,distance,ager+2,0,0)];
        Ztilde2nB(1,:,k) = repmat((1-               omega(1,j)),[1 15]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoCflWageInt(nloc,LYi,shock_amt,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervec(nloc,l +nloc,l +nloc,ager+2) movervec(nloc,l +nloc,l +nloc,distance,ager+2,0,0)];
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
    ZtildeTemp = Ztilde0eA+Ztilde0uA-Ztilde0B + Beta*(Ztilde1eA+Ztilde1uA-Ztilde1eB-Ztilde1uB-Ztilde1nB) + Beta.^2*(Ztilde2eA+Ztilde2uA-Ztilde2eB-Ztilde2uB-Ztilde2nB);
    Ztilde = ZtildeTemp+Ztilde;
end
Z = cat(2,Zloc,(1/draws)*Ztilde);

end
