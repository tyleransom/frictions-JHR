function lambda = lambdervecRhoCflPiInt(nloc,j,Tmat,Tvec,u,empStat,here,e,r,t,err,be,bu,shock_loc,shock_amt)
% OLD INPUTS: function lambda = lambdervecRho(nloc,j,y,e,l,r,t)
% This function calculates employment probabilities
if isscalar(j)
   j=kron(ones(size(e)),j); %Nx1 vector containing j in each element; otherwise it's Nx1 vector of different j's depending on the person
end
e(isnan(e))=0;
N = length(e);
D = size(err,2);
erep = repmat(e,[1 1 D]);

% input: Nx1 vector of row indices, Nx1 vector of column indices
% output: Nx1 vector of (row,column) elements
Umat  = zeros(size(Tvec,1),size(Tvec,2),D);
locer = zeros(size(e,1),nloc-1,D);
p     = zeros(size(e,1),1,D);
r1    = r(1,j)';
r2    = r(2,j)';
for jj=1:nloc
	if jj>1
		locer(j==jj,jj-1,:) = 1;
	end
	for tt=1:size(Tmat,2)
		flag = j==jj & Tmat(:,tt)==1;
        if isscalar(flag)
            if flag
                Umat(1,1,:) = repmat((r1(flag).*(t==1)+r1(flag).*(1+r2(flag)).*(t==2))+(r2(flag).^t).*(u(jj,tt)+(jj==shock_loc)*shock_amt),[1 1 D])+reshape(err(flag,:),[sum(flag) 1 D]);
            end
        else
            Umat(flag,1,:) = repmat((r1(flag).*(t==1)+r1(flag).*(1+r2(flag)).*(t==2))+(r2(flag).^t).*(u(jj,tt)+(jj==shock_loc)*shock_amt),[1 1 D])+reshape(err(flag,:),[sum(flag) 1 D]);
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
		if here==1
			X = cat(2,ones(N,1,D),locer,100*Umat,zeros(N,1,D),erep,erep.^2./100);
			for d=1:D
				p(:,:,d) = 1./(1+exp(squeeze(X(:,:,d))*be));
			end
		else
			X = cat(2,ones(N,1,D),locer,100*Umat,ones(N,1,D),erep,erep.^2./100);
			for d=1:D
				p(:,:,d) = exp(squeeze(X(:,:,d))*be)./(1+exp(squeeze(X(:,:,d))*be));
			end
		end
	else
		if here==1
			X = cat(2,ones(N,1,D),locer,100*Umat,zeros(N,1,D),erep,erep.^2./100);
			for d=1:D
				p(:,:,d) = exp(squeeze(X(:,:,d))*bu)./(1+exp(squeeze(X(:,:,d))*bu));
			end
		else
			X = cat(2,ones(N,1,D),locer,100*Umat,ones(N,1,D),erep,erep.^2./100);
			for d=1:D
				p(:,:,d) = exp(squeeze(X(:,:,d))*bu)./(1+exp(squeeze(X(:,:,d))*bu));
			end
		end
	end
    lambda = p;
end

end
