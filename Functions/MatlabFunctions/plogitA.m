function P=plogit(b,X,Z,nb,nb2,nc)
%setting up the parameters that are the same across choices
b2=b(nb*(nc-1)+1:nb*(nc-1)+nb2);

num=zeros(size(X,1),1);
dem=zeros(size(X,1),1);

%sets the last alternative to be the one that is normalized to zero
for j=1:nc-1
    temp=X*b((j-1)*nb+1:j*nb)+Z(:,:,j)*b2;
    num(:,j)=exp(temp);
    dem=exp(temp)+dem;
end
num=[num ones(size(X,1),1)];
dem=dem+1;

P=num./(dem*ones(1,nc));