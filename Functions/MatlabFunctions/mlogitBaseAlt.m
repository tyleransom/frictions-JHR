function [like,grad] = mlogitBaseAlt(b,Y,X,Z,nb,nb2,nc,baseAlt)
if nargin==7
   baseAlt = nc;
end
if ~isempty(Z);
    %setting up the parameters that are the same across choices
    N    = size(X,1);
    b2   = b(nb*(nc-1)+1:nb*(nc-1)+nb2);
    num  = zeros(N,1);
    dem  = zeros(N,1);
    k    = 1;

    %sets the last alternative to be the one that is normalized to zero
    for j=setdiff(1:nc,baseAlt)
	num=(Y==j).*(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)+num;
	%First dimension of Z is individuals, second dim is characteristics, third is choice alternatives
	dem=exp(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)+dem;
    k  = k+1;
    end
    dem=dem+1;

    like=sum(log(dem)-num);


    numg = zeros(nb2,1);
    demg = zeros(nb2,1);

    grad = zeros(size(b));
    k    = 1;
    for j=setdiff(1:nc,baseAlt)
	grad((k-1)*nb+1:k*nb)=-X'*((Y==j)-exp(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)./dem);
    k = k+1;
    end
    k = 1;
    for j=setdiff(1:nc,baseAlt)
	numg = -Z(:,:,j)'*(Y==j)+numg;
	demg = -Z(:,:,j)'*(exp(X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2)./dem)+demg;
    k    = k+1;
    end
    grad(nb*(nc-1)+1:nb*(nc-1)+nb2)=numg-demg;
else
    if nb2~=0
	error('Z is empty but nb2~=0!!!')
    end
    %setting up the parameters that are the same across choices
    N    = size(X,1);
    num  = zeros(N,1);
    dem  = zeros(N,1);
    k    = 1;

    %sets BASEALT to be normalized to zero
    for j=setdiff(1:nc,baseAlt)
	num=(Y==j).*(X*b((k-1)*nb+1:k*nb))+num;
	%First dimension of Z is individuals, second dim is characteristics, third is choice alternatives
	dem=exp(X*b((k-1)*nb+1:k*nb))+dem;
    k = k+1;
    end
    dem=dem+1;

    like=sum(log(dem)-num);

    grad = zeros(size(b));
    k    = 1;
    for j=setdiff(1:nc,baseAlt)
	grad((k-1)*nb+1:k*nb)=-X'*((Y==j)-exp(X*b((k-1)*nb+1:k*nb))./dem);
    k = k+1;
    end
end