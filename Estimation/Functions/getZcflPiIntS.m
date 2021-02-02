function [Z,pimat] = getZcflPiIntS(draws,ARcov,ARcov2,yr,ager,LEi,LYi,experer,typer,birthLoc,birthDiv,urate_lag55,J,rho_hat_urate,bLemp,bLunemp,rho_hat_wage,bwage,shock_amt,distance,Beta)



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
       pimat(1,k)  = lambdervecRhoCflPiS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,1,experer,typer,rho_hat_urate  ,0,bLemp,bLunemp,l,shock_amt);
    elseif LEi==1 && (k==LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = 1-lambdervecRhoCflPiS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,1,experer,typer,rho_hat_urate   ,0,bLemp,bLunemp,l,shock_amt);
    elseif LEi==0 && (k~=LYi-nloc*(LYi>nloc));
       pimat(1,k)  = lambdervecRhoCflPiS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,0,experer,typer,rho_hat_urate,0,bLemp,bLunemp,l,shock_amt);
    elseif LEi==1 && (k~=LYi-nloc*(LYi>nloc)); 
       pimat(1,k)  = lambdervecRhoCflPiS(nloc,k,Tmater,yr-2003,urate_lag55,LEi==1,0,experer,typer,rho_hat_urate,0,bLemp,bLunemp,l,shock_amt);
    end
end

Zloc = zeros(1,nloc,J);
Ztilde    = zeros(1,17,J);
Ztilde0eA = zeros(1,17,J);
Ztilde1eA = zeros(1,17,J);
Ztilde2eA = zeros(1,17,J);
Ztilde0uA = zeros(1,17,J);
Ztilde1uA = zeros(1,17,J);
Ztilde2uA = zeros(1,17,J);
Ztilde0B  = zeros(1,17,J);
Ztilde1eB = zeros(1,17,J);
Ztilde2eB = zeros(1,17,J);
Ztilde1uB = zeros(1,17,J);
Ztilde2uB = zeros(1,17,J);
Ztilde1nB = zeros(1,17,J);
Ztilde2nB = zeros(1,17,J);
shocker   = nan(1,2*nloc);
for dd=1:draws
    shocker(1,:,1) = mvnrnd(zeros(1,2*nloc),ARcov);
	for jj=1:nloc
        pimat1(1,1) = squeeze(lambdervecRhoCflPiIntS(nloc,l,Tmater,yr-2003,urate_lag55,0,1,experer,typer,rho_hat_urate,1,squeeze(shocker(1,(l-1)*2+2,:))',bLemp,bLunemp,l,shock_amt));
        omega(1,jj) = pimat(1,jj)./pimat1(1,1);
	end
	for j=1:nloc
        lp = j;
        k=j+nloc;
        Ztilde0eA(1,:,j) = repmat((  pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,j      ,experer) ewagevecRhoIntS(nloc,j      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,j      ) birthDiv(1,j      ) switchervecS(nloc,LYi    ,j      ,ager+0,typer) movervecS(nloc,LYi    ,j      ,distance,ager+0,pvEmp,pvUnemp,typer)];
        Ztilde0eA(1,:,k) = repmat((  pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoIntS(nloc,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,k      ) birthDiv(1,k      ) switchervecS(nloc,LYi    ,k      ,ager+0,typer) movervecS(nloc,LYi    ,k      ,distance,ager+0,pvEmp,pvUnemp,typer)];
        Ztilde0uA(1,:,j) = repmat((1-pimat(1,j)),[1 17])               .*[1*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoIntS(nloc,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,j      ) birthDiv(1,j      ) switchervecS(nloc,LYi    ,j      ,ager+0,typer) movervecS(nloc,LYi    ,j      ,distance,ager+0,pvEmp,pvUnemp,typer)];
        Ztilde0uA(1,:,k) = repmat((1-pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,k      ,experer) ewagevecRhoIntS(nloc,k      ,Tmater,yr-2003,rho_hat_wage,0,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,k      ) birthDiv(1,k      ) switchervecS(nloc,LYi    ,k      ,ager+0,typer) movervecS(nloc,LYi    ,k      ,distance,ager+0,pvEmp,pvUnemp,typer)];
        Ztilde1eA(1,:,j) = repmat((  pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoIntS(nloc,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer+1,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervecS(nloc,j      ,lp+nloc,ager+1,typer) movervecS(nloc,j      ,lp+nloc,distance,ager+1,1,0,typer)];
        Ztilde1eA(1,:,k) = repmat((  pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoIntS(nloc,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervecS(nloc,k      ,lp+nloc,ager+1,typer) movervecS(nloc,k      ,lp+nloc,distance,ager+1,1,0,typer)];
        Ztilde1uA(1,:,j) = repmat((1-pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoIntS(nloc,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervecS(nloc,j      ,lp+nloc,ager+1,typer) movervecS(nloc,j      ,lp+nloc,distance,ager+1,0,1,typer)];
        Ztilde1uA(1,:,k) = repmat((1-pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,lp+nloc,experer) ewagevecRhoIntS(nloc,lp+nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,lp+nloc) birthDiv(1,lp+nloc) switchervecS(nloc,k      ,lp+nloc,ager+1,typer) movervecS(nloc,k      ,lp+nloc,distance,ager+1,0,1,typer)];
        Ztilde2eA(1,:,j) = repmat((  pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer+1,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,lp+nloc,l +nloc,ager+2,typer) movervecS(nloc,lp+nloc,l +nloc,distance,ager+2,0,0,typer)];
        Ztilde2eA(1,:,k) = repmat((  pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,lp+nloc,l +nloc,ager+2,typer) movervecS(nloc,lp+nloc,l +nloc,distance,ager+2,0,0,typer)];
        Ztilde2uA(1,:,j) = repmat((1-pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,lp+nloc,l +nloc,ager+2,typer) movervecS(nloc,lp+nloc,l +nloc,distance,ager+2,0,0,typer)];
        Ztilde2uA(1,:,k) = repmat((1-pimat(1,j)),[1 17])               .*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(j-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,lp+nloc,l +nloc,ager+2,typer) movervecS(nloc,lp+nloc,l +nloc,distance,ager+2,0,0,typer)];

        Ztilde0B(1,:,j)  =                                               [0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,0,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,LYi    ,k      ,ager+0,typer) movervecS(nloc,LYi    ,l +nloc,distance,ager+0,pvEmp,pvUnemp,typer)];
        Ztilde0B(1,:,k)  =                                               [0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,0,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,LYi    ,k      ,ager+0,typer) movervecS(nloc,LYi    ,l +nloc,distance,ager+0,pvEmp,pvUnemp,typer)];
        Ztilde1eB(1,:,j) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l      ,experer) ewagevecRhoIntS(nloc,l      ,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l      ) birthDiv(1,l      ) switchervecS(nloc,l +nloc,l      ,ager+1,typer) movervecS(nloc,l +nloc,l      ,distance,ager+1,0,0,typer)];
        Ztilde1eB(1,:,k) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+1,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+1,0,0,typer)];
        Ztilde1uB(1,:,j) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 17]).*[1*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l      ,ager+1,typer) movervecS(nloc,l +nloc,l      ,distance,ager+1,0,0,typer)];
        Ztilde1uB(1,:,k) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+1,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+1,0,0,typer)];
        Ztilde1nB(1,:,j) = repmat((1-               omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+1,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+1,0,0,typer)];
        Ztilde1nB(1,:,k) = repmat((1-               omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,1,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+1,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+1,0,0,typer)];
        Ztilde2eB(1,:,j) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer+1,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l      ,l +nloc,ager+2,typer) movervecS(nloc,l      ,l +nloc,distance,ager+2,1,0,typer)];
        Ztilde2eB(1,:,k) = repmat((  pimat1(1,1)).*(omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+2,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+2,1,0,typer)];
        Ztilde2uB(1,:,j) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l      ,l +nloc,ager+2,typer) movervecS(nloc,l      ,l +nloc,distance,ager+2,0,1,typer)];
        Ztilde2uB(1,:,k) = repmat((1-pimat1(1,1)).*(omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+2,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+2,0,1,typer)];
        Ztilde2nB(1,:,j) = repmat((1-               omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+2,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+2,0,0,typer)];
        Ztilde2nB(1,:,k) = repmat((1-               omega(1,j)),[1 17]).*[0*ones(1,1) benefits(nloc,l +nloc,experer) ewagevecRhoIntS(nloc,l +nloc,Tmater,yr-2003,rho_hat_wage,2,experer  ,typer,bwage,squeeze(shocker(1,(l-1)*2+1,:))) birthLoc(1,l +nloc) birthDiv(1,l +nloc) switchervecS(nloc,l +nloc,l +nloc,ager+2,typer) movervecS(nloc,l +nloc,l +nloc,distance,ager+2,0,0,typer)];
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
