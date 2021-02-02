function lambda = lambdervecRhoCflPi(nloc,j,Tmat,Tvec,u,empStat,here,e,r,t,be,bu,shock_loc,shock_amt)
% OLD INPUTS: function lambda = lambdervecRho(nloc,j,y,e,l,r,t)
% This function calculates employment probabilities
if isscalar(j)
   j=kron(ones(size(e)),j); %Nx1 vector containing j in each element; otherwise it's Nx1 vector of different j's depending on the person
end
e(isnan(e))=0;

% input: Nx1 vector of row indices, Nx1 vector of column indices
% output: Nx1 vector of (row,column) elements
Umat  = zeros(size(Tvec));
locer = zeros(size(e,1),nloc-1);
r1     = r(1,j)';
r2     = r(2,j)';
for jj=1:nloc
	if jj>1
		locer(j==jj,jj-1) = 1;
	end
	for tt=1:size(Tmat,2)
		flag = j==jj & Tmat(:,tt)==1;
		Umat(flag) = (r1(flag).*(t==1)+r1(flag).*(1+r2(flag)).*(t==2))+(r2(flag).^t).*(u(jj,tt)+(jj==shock_loc)*shock_amt);
	end
end
if all(j>nloc)
	lambda = zeros(size(e));
else
	% NEED TO FIGURE OUT:
	% how to vectorize the u(j,t) reference [over the j's]
	% how to vectorize the r(x,j) reference [over the j's]
	
	if empStat==1
		if here==1
			X = [ones(size(e)) locer 100*Umat zeros(size(e)) e e.^2./100];
			p = 1./(1+exp(X*be));
		else
			X = [ones(size(e)) locer 100*Umat ones(size(e)) e e.^2./100];
			p = exp(X*be)./(1+exp(X*be));
		end
	else
		if here==1
			X = [ones(size(e)) locer 100*Umat zeros(size(e)) e e.^2./100];
			p = exp(X*bu)./(1+exp(X*bu));
		else
			X = [ones(size(e)) locer 100*Umat ones(size(e)) e e.^2./100];
			p = exp(X*bu)./(1+exp(X*bu));
		end
	end
    lambda = p;
end

end