function [like,grad] = mlogit_restrict_base_alt_exclZ(b,restrMat,Y,X,Z,nb,nb2,nc,baseAlt,exZcol,exZalt)
b = applyRestr(restrMat,b);
if nargin==8
    baseAlt = nc;
    exZcol = [];
    exZalt = [];
end
if ~isempty(Z);
    %setting up the parameters that are the same across choices
    N    = size(X,1);
    b2   = b(nb*(nc-1)+1:nb*(nc-1)+nb2);
    num  = zeros(N,1);
    dem  = zeros(N,1);
    k    = 1;
    
    %sets the base alternative to be the one that is normalized to zero
    for j=setdiff(1:nc,baseAlt)
        if ismember(j,exZalt)
            flag = setdiff(1:nb2,exZcol);
            num=(Y==j).*(X*b((k-1)*nb+1:k*nb)+Z(:,flag,j)*b2(flag))+num;
            %First dimension of Z is individuals, second dim is characteristics, third is choice alternatives
            dem=exp(X*b((k-1)*nb+1:k*nb)+Z(:,flag,j)*b2(flag))+dem;
            k = k+1;
        else
            num=(Y==j).*(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)+num;
            %First dimension of Z is individuals, second dim is characteristics, third is choice alternatives
            dem=exp(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)+dem;
            k = k+1;
        end
    end
    dem=dem+1;
    
    like=sum(log(dem)-num);
    
    numg = zeros(nb2,1);
    demg = zeros(nb2,1);
    
    grad = zeros(size(b));
    k    = 1;
    for j=setdiff(1:nc,baseAlt)
        if ismember(j,exZalt)
            flag  = setdiff(1:nb2,exZcol);
            grad((k-1)*nb+1:k*nb)=-X'*((Y==j)-exp(X*b((k-1)*nb+1:k*nb)+Z(:,flag,j)*b2(flag))./dem);
            k = k+1;
        else
            grad((k-1)*nb+1:k*nb)=-X'*((Y==j)-exp(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)./dem);
            k = k+1;
        end
    end
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        if ismember(j,exZalt)
            numg1 = zeros(nb2,1);
            demg1 = zeros(nb2,1);
            flag  = setdiff(1:nb2,exZcol);
            numg1(flag) = -Z(:,flag,j)'*(Y==j);
            demg1(flag) = -Z(:,flag,j)'*(exp(X*b((k-1)*nb+1:k*nb)+Z(:,flag,j)*b2(flag))./dem);
            numg  = numg1+numg;
            demg  = demg1+demg;
            k = k+1;
        else
            numg = -Z(:,:,j)'*(Y==j)+numg;
            demg = -Z(:,:,j)'*(exp(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)./dem)+demg;
            k = k+1;
        end
    end
    grad(nb*(nc-1)+1:nb*(nc-1)+nb2)=numg-demg;
    grad = applyRestrGrad(restrMat,grad);
else
    if nb2~=0
        error('Z is empty but nb2~=0!!!')
    end
    %setting up the parameters that are the same across choices
    N    = size(X,1);
    num  = zeros(N,1);
    dem  = zeros(N,1);
    k    = 1;
    
    %sets the last alternative to be the one that is normalized to zero
    for j=setdiff(1:nc,baseAlt)
        num=(Y==j).*(X*b((k-1)*nb+1:k*nb))+num;
        %First dimension of Z is individuals, second dim is characteristics, third is choice alternatives
        dem=exp(X*b((k-1)*nb+1:k*nb))+dem;
        k  = k+1;
    end
    dem=dem+1;
    
    like=sum(log(dem)-num);
    
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        P(:,j) = exp(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)./dem;
        k = k+1;
    end
    
    grad = zeros(size(b));
    k    = 1;
    for j=setdiff(1:nc,baseAlt)
        grad((k-1)*nb+1:k*nb)=-X'*((Y==j)-exp(X*b((k-1)*nb+1:k*nb))./dem);
        k = k+1;
    end
    grad = applyRestrGrad(restrMat,grad);
end