function [Z,pimat] = constructZ(N,T,J,nloc,draws,flag,flagw,ARcov,ARcov2,calyr,Tmat,Tvec,urate_lag55,rho_hat_urate,rho_hat_wage,bLemp,bLunemp,bwage,exper,age,distance,birthLoc,birthDiv,LY,homeLoc,pvEmp,pvUnemp,empFT_lag,Beta)
	%seed = 1234;
	%rng(seed,'twister');

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

	for k=1:nloc
		if flag & empFT_lag(:)==0 & (k==LY(:)-nloc*(LY(:)>nloc))
		pimat(1,k) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,1,exper,rho_hat_urate  ,0,bLemp,bLunemp);

        elseif flag & empFT_lag(:)==1 & (k==LY(:)-nloc*(LY(:)>nloc));
		pimat(1,k) = 1-lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,1,exper,rho_hat_urate   ,0,bLemp,bLunemp);

        elseif flag & empFT_lag(:)==0 & (k~=LY(:)-nloc*(LY(:)>nloc));
		pimat(1,k) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,0,exper,rho_hat_urate,0,bLemp,bLunemp);

        elseif flag & empFT_lag(:)==1 & (k~=LY(:)-nloc*(LY(:)>nloc));
		pimat(1,k) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,empFT_lag==1,0,exper,rho_hat_urate,0,bLemp,bLunemp);
		
        elseif flag & (k==LY(:)-nloc*(LY(:)>nloc));
		pimat1(1,1) =   lambdervecRho(nloc,k,Tmat,Tvec,urate_lag55,0,1,exper,rho_hat_urate,1,bLemp,bLunemp);
        end
	end


	Zloc = zeros(N*T,nloc,J);
	Ztilde    = zeros(N*T,15,J);
	Ztilde0eA = zeros(N*T,15,J);
	Ztilde1eA = zeros(N*T,15,J);
	Ztilde2eA = zeros(N*T,15,J);
	Ztilde0uA = zeros(N*T,15,J);
	Ztilde1uA = zeros(N*T,15,J);
	Ztilde2uA = zeros(N*T,15,J);
	Ztilde0B  = zeros(N*T,15,J);
	Ztilde1eB = zeros(N*T,15,J);
	Ztilde2eB = zeros(N*T,15,J);
	Ztilde1uB = zeros(N*T,15,J);
	Ztilde2uB = zeros(N*T,15,J);
	Ztilde1nB = zeros(N*T,15,J);
	Ztilde2nB = zeros(N*T,15,J);
	shocker   = nan(N*T,2*nloc);
	shocklee  = nan(N*T,1);
	for dd=1:draws
		shocker(flag,:,1) = mvnrnd(zeros(sum(flag),2*nloc),ARcov);
		for i=1:N*T %:numel(exper) % loop over individual-time observations
		    if flag(i)==1
		        LYi = LY(i);
		        ti  = calyr(i)-2003;
		        l = LYi-nloc*(LYi>nloc);
		        % Create pi_ikt and pi_it+1 for all locations
		        for jj=1:nloc
		            pimat1(i,1) = squeeze(lambdervecRhoInt(nloc,l,Tmat(i,:),Tvec(i),urate_lag55,0,1,exper(i),rho_hat_urate,1,squeeze(shocker(i,(l-1)*2+2,:))',bLemp,bLunemp)); % f(rho,lambda_hat(l)); %not sure how log-odds AR(1) interfaces here
		            omega(i,jj) = pimat(i,jj)./pimat1(i,1);
		        end
		        shocklee(i) = squeeze(shocker(i,(l-1)*2+1,:));
		    end
		end
		for j=1:nloc
		    lp = j;
		    k=j+nloc;
		    Ztilde0eA(flag,:,j) = repmat((  pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,j                 ,exper(flagw)) ewagevecRhoInt(nloc,j                 ,Tmat(flag,:),Tvec(flag),rho_hat_wage,0,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,j                       ) birthDiv(flag,j                       ) switchervec(nloc,LY(flagw)         ,j                 ,age(flagw)+0) movervec(nloc,LY(flagw)         ,j                 ,distance,age(flagw)+0,             pvEmp,           pvUnemp)];
		    Ztilde0eA(flag,:,k) = repmat((  pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,k                 ,exper(flagw)) ewagevecRhoInt(nloc,k                 ,Tmat(flag,:),Tvec(flag),rho_hat_wage,0,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,k                       ) birthDiv(flag,k                       ) switchervec(nloc,LY(flagw)         ,k                 ,age(flagw)+0) movervec(nloc,LY(flagw)         ,k                 ,distance,age(flagw)+0,             pvEmp,           pvUnemp)];
		    Ztilde0uA(flag,:,j) = repmat((1-pimat(flag,j)),[1 15])                  .*[1*ones(sum(flag),1) benefits(nloc,k                 ,exper(flagw)) ewagevecRhoInt(nloc,k                 ,Tmat(flag,:),Tvec(flag),rho_hat_wage,0,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,j                       ) birthDiv(flag,j                       ) switchervec(nloc,LY(flagw)         ,j                 ,age(flagw)+0) movervec(nloc,LY(flagw)         ,j                 ,distance,age(flagw)+0,             pvEmp,           pvUnemp)];
		    Ztilde0uA(flag,:,k) = repmat((1-pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,k                 ,exper(flagw)) ewagevecRhoInt(nloc,k                 ,Tmat(flag,:),Tvec(flag),rho_hat_wage,0,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,k                       ) birthDiv(flag,k                       ) switchervec(nloc,LY(flagw)         ,k                 ,age(flagw)+0) movervec(nloc,LY(flagw)         ,k                 ,distance,age(flagw)+0,             pvEmp,           pvUnemp)];
		    Ztilde1eA(flag,:,j) = repmat((  pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,lp+nloc           ,exper(flagw)) ewagevecRhoInt(nloc,lp           +nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)+1,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,lp+nloc                 ) birthDiv(flag,lp+nloc                 ) switchervec(nloc,j                 ,lp+nloc           ,age(flagw)+1) movervec(nloc,j                 ,lp+nloc           ,distance,age(flagw)+1, ones(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde1eA(flag,:,k) = repmat((  pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,lp+nloc           ,exper(flagw)) ewagevecRhoInt(nloc,lp           +nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,lp+nloc                 ) birthDiv(flag,lp+nloc                 ) switchervec(nloc,k                 ,lp+nloc           ,age(flagw)+1) movervec(nloc,k                 ,lp+nloc           ,distance,age(flagw)+1, ones(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde1uA(flag,:,j) = repmat((1-pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,lp+nloc           ,exper(flagw)) ewagevecRhoInt(nloc,lp           +nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,lp+nloc                 ) birthDiv(flag,lp+nloc                 ) switchervec(nloc,j                 ,lp+nloc           ,age(flagw)+1) movervec(nloc,j                 ,lp+nloc           ,distance,age(flagw)+1,zeros(sum(flag),1), ones(sum(flag),1))];
		    Ztilde1uA(flag,:,k) = repmat((1-pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,lp+nloc           ,exper(flagw)) ewagevecRhoInt(nloc,lp           +nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) birthLoc(flag,lp+nloc                 ) birthDiv(flag,lp+nloc                 ) switchervec(nloc,k                 ,lp+nloc           ,age(flagw)+1) movervec(nloc,k                 ,lp+nloc           ,distance,age(flagw)+1,zeros(sum(flag),1), ones(sum(flag),1))];
		    Ztilde2eA(flag,:,j) = repmat((  pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)+1,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1), ones(sum(flag),1))];
		    Ztilde2eA(flag,:,k) = repmat((  pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1), ones(sum(flag),1))];
		    Ztilde2uA(flag,:,j) = repmat((1-pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1), ones(sum(flag),1))];
		    Ztilde2uA(flag,:,k) = repmat((1-pimat(flag,j)),[1 15])                  .*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,lp+nloc           ,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1), ones(sum(flag),1))];

		    Ztilde0B(flag,:,j)  =                                                     [0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,0,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,LY(flagw)         ,k                 ,age(flagw)+0) movervec(nloc,LY(flagw)         ,homeLoc(flag)+nloc,distance,age(flagw)+0,             pvEmp,           pvUnemp)];
		    Ztilde0B(flag,:,k)  =                                                     [0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,0,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,LY(flagw)         ,k                 ,age(flagw)+0) movervec(nloc,LY(flagw)         ,homeLoc(flag)+nloc,distance,age(flagw)+0,             pvEmp,           pvUnemp)];
		    Ztilde1eB(flag,:,j) = repmat((  pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)     ,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)     ,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,shocklee(flag)                    ) max(birthLoc(flag,:).*LwS(flag,:),[],2) max(birthDiv(flag,:).*LwS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)     ,age(flagw)+1) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)     ,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde1eB(flag,:,k) = repmat((  pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+1) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde1uB(flag,:,j) = repmat((1-pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[1*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)     ,age(flagw)+1) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)     ,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde1uB(flag,:,k) = repmat((1-pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+1) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde1nB(flag,:,j) = repmat((1-                  omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+1) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde1nB(flag,:,k) = repmat((1-                  omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,1,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+1) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+1,zeros(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde2eB(flag,:,j) = repmat((  pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)+1,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)     ,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,homeLoc(flag)     ,homeLoc(flag)+nloc,distance,age(flagw)+2, ones(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde2eB(flag,:,k) = repmat((  pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+2, ones(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde2uB(flag,:,j) = repmat((1-pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)     ,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,homeLoc(flag)     ,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1), ones(sum(flag),1))];
		    Ztilde2uB(flag,:,k) = repmat((1-pimat1(flag,1)).*(omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1), ones(sum(flag),1))];
		    Ztilde2nB(flag,:,j) = repmat((1-                  omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1))];
		    Ztilde2nB(flag,:,k) = repmat((1-                  omega(flag,j)),[1 15]).*[0*ones(sum(flag),1) benefits(nloc,homeLoc(flag)+nloc,exper(flagw)) ewagevecRhoInt(nloc,homeLoc(flag)+nloc,Tmat(flag,:),Tvec(flag),rho_hat_wage,2,exper(flagw)  ,bwage,squeeze(shocker(flag,(j-1)*2+1,:))) max(birthLoc(flag,:).*LhS(flag,:),[],2) max(birthDiv(flag,:).*LhS(flag,:),[],2) switchervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,age(flagw)+2) movervec(nloc,homeLoc(flag)+nloc,homeLoc(flag)+nloc,distance,age(flagw)+2,zeros(sum(flag),1),zeros(sum(flag),1))];
		    for i=1:N*T
		        if flag(i)==1
		            LYi = LY(i);
		            l = LYi-nloc*(LYi>nloc);
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
	Z = cat(2,Zloc,(1/draws)*Ztilde);
	assert(size(Zloc,2)==nloc,'Zloc incorrectly constructed');
end
