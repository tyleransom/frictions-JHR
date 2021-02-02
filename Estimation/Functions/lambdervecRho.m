function lambda = lambdervecRho(nloc,j,Tmat,Tvec,u,empStat,here,e,r,t,be,bu)
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
		Umat(flag) = (r1(flag).*(t==1)+r1(flag).*(1+r2(flag)).*(t==2))+(r2(flag).^t).*u(jj,tt);
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


% % Wage function
% if isscalar(j)
%    j=kron(ones(size(e)),j); %Nx1 vector containing j in each element; otherwise it's Nx1 vector of different j's depending on the person
% end
% e(isnan(e))=0;
% 
% if all(j>nloc)
%     w=zeros(size(e));
% else
%     r1 = r(1,j)';
%     r2 = r(2,j)';
%     if ~isscalar(r2);
%         r2temp = r2*ones(1,9); 
%     else
%         r2temp = r2;
%     end
%     locer = zeros(length(e),nloc-1);
%     timer = (r2temp.^t).*Tmat(:,2:end);
%     loctimer = zeros(length(e),(nloc-1)*9);
%     if numel(unique(j))==1
%         jj = mean(j);
%         if j>1
%             locer(:,jj-1) = ones(size(e));
%             for k=2:size(Tmat,2)
%                 loctimer(:,(jj-2)*(size(Tmat,2)-1)+k-1) = (r2.^t).*(ones(size(e)) & Tvec==k);
%             end
%         end
%     else
%         for jj=setdiff(1:nloc,1);
%             locer(:,jj-1) = logical(ones(size(e)).*(j==jj));
%         end
%         for jj=setdiff(1:nloc,1); %baseLoc = 1
%             for k=2:size(Tmat,2)
%                 loctimer(:,(jj-2)*(size(Tmat,2)-1)+k-1) = (r2.^t).*(logical(ones(size(e)).*(j==jj)) & Tvec==k);
%             end
%         end
%     end
% 
%     if nloc==55 || nloc==50 || nloc==29
%         loctimer(:,end-17:end)=[];
%     end
%     if nloc==54 || nloc==49 || nloc==28
%         loctimer(:,end-8:end)=[];
%     end
%     w = [ones(size(e)) locer timer loctimer e e.^2./100]*b;
%     w = w+(t==1).*r1 + (t==2).*(1+r2).*r1;
% end
% end
