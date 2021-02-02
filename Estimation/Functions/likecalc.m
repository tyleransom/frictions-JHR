function [full] = likecalc(dat,parms)

N      = dat.N;
T      = dat.T;
S      = dat.S;
J      = dat.J;

ID = repmat([1:N]',[1 T S]);
time = repmat(1:T,[N 1 S]);
typer = zeros(size(dat.typeS));
for s=1:S
	typer(:,:,s) = s;
end

% Initialize matrices
Pe     = ones(N,T,S);
Pu     = ones(N,T,S);
wRes   = zeros(N,T,S);

% Fill in each likelihood component:
Pe = dat.Pemp;
Pu = dat.Punemp;
wResTemp = dat.wage(:) - dat.Xw*parms.wageBeta;
wRes = reshape(wResTemp,[N T S]);
wRes(~dat.wageflagStemp) = -999;

%for j=1:J
%	tempj  = dat.Pc(:,j);
%	PcTemp = ones(N,T,S);
%	PcTemp(dat.flagS==1) = tempj;
%	if j==1
%		Pc = PcTemp;
%	else
%		Pc = cat(4,Pc,PcTemp);
%	end
%end

Pc = zeros(N,T,S,J);
Choicer = zeros(N,T,S,J);
Flagger = zeros(N,T,S,J);
for j=1:J
	Pc(:,:,:,j) = reshape(dat.Pc(:,j),[N T S 1]);
    Choicer(:,:,:,j) = dat.ChoiceS==j;
    Flagger(:,:,:,j) = dat.flagS;
end
Pc(isnan(Pc))=0;

emp_like = ones(N,S);
for s=1:S
	emp_like(:,s)   = prod((squeeze(Pe(:,:,s)).^(dat.empFTS(:,:,s)==1).*((1-squeeze(Pe(:,:,s))).^(1-(dat.empFTS(:,:,s)==1)))).^(dat.inlfS(:,:,s)==1 & dat.empFT_lagS(:,:,s)==1 & dat.flagS(:,:,s)==1),2);
end

unemp_like = ones(N,S);
for s=1:S
	unemp_like(:,s) = prod((squeeze(Pu(:,:,s)).^(dat.empFTS(:,:,s)==1).*((1-squeeze(Pu(:,:,s))).^(1-(dat.empFTS(:,:,s)==1)))).^(dat.inlfS(:,:,s)==1 & dat.empFT_lagS(:,:,s)==0 & dat.flagS(:,:,s)==1),2);
end

wage_like = ones(N,S);
for s=1:S
	wage_like(:,s) = prod(normpdf(wRes(:,:,s),0,parms.wageSig).^(dat.wageflagStemp(:,:,s)==1),2);
end

choice_like = ones(N,S);
for s=1:S
	choice_like(:,s) = prod(prod(Pc(:,:,s,:).^(Choicer(:,:,s,:)==1 & Flagger(:,:,s,:)==1),4),2);
end

full = emp_like.*unemp_like.*wage_like.*choice_like;

end
