function [like,grad]=mlogit(b,Y,X,Z,nb,nb2,nc)
%setting up the parameters that are the same across choices
N    = size(X,1);
b2   = b(nb*(nc-1)+1:nb*(nc-1)+nb2);
num  = zeros(N,1);
dem  = zeros(N,1);

%sets the last alternative to be the one that is normalized to zero
for j=1:nc-1
    num=(Y==j).*(X*b((j-1)*nb+1:j*nb)+Z(:,:,j)*b2)+num;
	%First dimension of Z is individuals, second dim is characteristics, third is choices
    dem=exp(X*b((j-1)*nb+1:j*nb)+Z(:,:,j)*b2)+dem;
end
dem=dem+1;

like=sum(log(dem)-num);


numg = zeros(nb2,1);
demg = zeros(nb2,1);

grad = zeros(size(b));
for j=1:nc-1
    grad((j-1)*nb+1:j*nb)=-X'*((Y==j)-exp(X*b((j-1)*nb+1:j*nb)+Z(:,:,j)*b2)./dem);
end
for j=1:nc-1
	numg = -Z(:,:,j)'*(Y==j)+numg;
    demg = -Z(:,:,j)'*(exp(X*b((j-1)*nb+1:j*nb)+Z(:,:,j)*b2)./dem)+demg;
end
grad(nb*(nc-1)+1:nb*(nc-1)+nb2)=numg-demg;