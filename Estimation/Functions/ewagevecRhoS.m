function w = ewagevecRhoS(nloc,j,Tmat,Tvec,r,t,e,type,b)
% This function computes expected wages given state variables
if isscalar(j)
   j=kron(ones(size(e)),j); 
end
e(isnan(e))=0;

if all(j>nloc)
    w=zeros(size(e));
else
    r1 = r(1,j)';
    r2 = r(2,j)';
    if ~isscalar(r2);
        r2temp = r2*ones(1,9); 
    else
        r2temp = r2;
    end
    locer = zeros(length(e),nloc-1);
    timer = (r2temp.^t).*Tmat(:,2:end);
    loctimer = zeros(length(e),(nloc-1)*9);
    if numel(unique(j))==1
        jj = mean(j);
        if j>1
            locer(:,jj-1) = ones(size(e));
            for k=2:size(Tmat,2)
                loctimer(:,(jj-2)*(size(Tmat,2)-1)+k-1) = (r2.^t).*(ones(size(e)) & Tvec==k);
            end
        end
    else
        for jj=setdiff(1:nloc,1);
            locer(:,jj-1) = logical(ones(size(e)).*(j==jj));
        end
        for jj=setdiff(1:nloc,1); %baseLoc = 1
            for k=2:size(Tmat,2)
                loctimer(:,(jj-2)*(size(Tmat,2)-1)+k-1) = (r2.^t).*(logical(ones(size(e)).*(j==jj)) & Tvec==k);
            end
        end
    end

    if nloc==55 || nloc==50 || nloc==29
        loctimer(:,end-17:end)=[];
    end
    if nloc==54 || nloc==49 || nloc==28
        loctimer(:,end-8:end)=[];
    end
    w = [ones(size(e)) locer timer loctimer e e.^2./100 type==1]*b;
    w = w+(t==1).*r1 + (t==2).*(1+r2).*r1;
end
end
