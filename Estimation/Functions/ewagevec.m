function w = ewagevec(nloc,j,Tmat,Tvec,e,b)
% This function computes expected wages given state variables
if isscalar(j)
   j=kron(ones(size(e)),j); 
end
e(isnan(e))=0;
if all(j>nloc)
    w=zeros(size(e));
else
    locer = zeros(length(e),nloc-1);
    timer = Tmat(:,2:end);
    loctimer = zeros(length(e),(nloc-1)*9);
    if numel(unique(j))==1
        jj = mean(j);
        if j>1
            locer(:,jj-1) = ones(size(e));
            for k=2:size(Tmat,2)
                loctimer(:,(jj-2)*(size(Tmat,2)-1)+k-1) = ones(size(e)) & Tvec==k;
            end
        end
    else
        for jj=setdiff(1:nloc,1);
            locer(:,jj-1) = logical(ones(size(e)).*(j==jj));
        end
        for jj=setdiff(1:nloc,1); %baseLoc = 1
            for k=2:size(Tmat,2)
                loctimer(:,(jj-2)*(size(Tmat,2)-1)+k-1) = logical(ones(size(e)).*(j==jj)) & Tvec==k;
            end
        end
    end

    if nloc==55 || nloc==50 || nloc==29
        loctimer(:,end-17:end)=[];
    end
    if nloc==54 || nloc==49 || nloc==28
        loctimer(:,end-8:end)=[];
    end
    w = [ones(size(e)) locer timer loctimer e e.^2./100]*b;
end
end