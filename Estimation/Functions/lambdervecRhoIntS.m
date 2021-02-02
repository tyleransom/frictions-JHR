function lambda = lambdervecRhoIntS(nloc,j,Tmat,Tvec,u,empStat,here,e,type,r,t,err,be,bu)
% OLD INPUTS: function lambda = lambdervecRho(nloc,j,y,e,l,r,t)
% This function calculates employment probabilities
if isscalar(j)
   j=kron(ones(size(e)),j); %Nx1 vector containing j in each element; otherwise it's Nx1 vector of different j's depending on the person
end
e(isnan(e))=0;
N = length(e);
D = size(err,2);
erep = repmat(e,[1 1 D]);
typerep = repmat(type,[1 1 D]);

% input: Nx1 vector of row indices, Nx1 vector of column indices
% output: Nx1 vector of (row,column) elements
Umat  = zeros(size(Tvec,1),size(Tvec,2),D);
locer = zeros(size(e,1),nloc-1,D);
p     = zeros(size(e,1),1,D);
r1    = r(1,j)';
r2    = r(2,j)';
if N>1
    for jj=1:nloc
        if jj>1
            locer(j==jj,jj-1,:) = 1;
        end
        for tt=1:size(Tmat,2)
            flag = j==jj & Tmat(:,tt)==1;
            Umat(flag,1,:) = repmat((r1(flag).*(t==1)+r1(flag).*(1+r2(flag)).*(t==2))+(r2(flag).^t).*u(jj,tt),[1 1 D])+reshape(err(flag,:),[sum(flag) 1 D]);
        end
    end
elseif N==1
    for jj=1:nloc
        if jj>1
            locer(j==jj,jj-1,:) = 1;
        end
        for tt=1:size(Tmat,2)
            flag = j==jj & Tmat(:,tt)==1;
            if flag
                Umat(flag,1,:) = repmat((r1(flag).*(t==1)+r1(flag).*(1+r2(flag)).*(t==2))+(r2(flag).^t).*u(jj,tt),[1 1 D])+reshape(err(flag,:),[sum(flag) 1 D]);
            end
        end
    end
end
if all(j>nloc)
	lambda = zeros(N,1,D);
else
	% NEED TO FIGURE OUT:
	% how to vectorize the u(j,t) reference [over the j's]
	% how to vectorize the r(x,j) reference [over the j's]
	
	if empStat==1
        K = length(be);
		if here==1
			X = cat(2,ones(N,1,D),locer,100*Umat,zeros(N,1,D),erep,erep.^2./100,typerep);
			linpred = reshape(reshape(permute(X,[2 1 3]),K,[]).'*be,N,1,[]);
            p = 1./(1+exp(linpred));
		else
			X = cat(2,ones(N,1,D),locer,100*Umat,ones(N,1,D),erep,erep.^2./100,typerep);
			linpred = reshape(reshape(permute(X,[2 1 3]),K,[]).'*be,N,1,[]);
            p = exp(linpred)./(1+exp(linpred));
		end
    else
        K = length(bu);
		if here==1
			X = cat(2,ones(N,1,D),locer,100*Umat,zeros(N,1,D),erep,erep.^2./100,typerep);
			linpred = reshape(reshape(permute(X,[2 1 3]),K,[]).'*bu,N,1,[]);
            p = exp(linpred)./(1+exp(linpred));
		else
			X = cat(2,ones(N,1,D),locer,100*Umat,ones(N,1,D),erep,erep.^2./100,typerep);
			linpred = reshape(reshape(permute(X,[2 1 3]),K,[]).'*bu,N,1,[]);
            p = exp(linpred)./(1+exp(linpred));
		end
	end
    lambda = p;
end

end
