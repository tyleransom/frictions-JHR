function P=plogitBaseAlt(b,X,Z,nb,nb2,nc,baseAlt)
if nargin==6
   baseAlt = nc; 
end
%setting up the parameters that are the same across choices
b2=b(nb*(nc-1)+1:nb*(nc-1)+nb2);

num=zeros(size(X,1),nc);
dem=zeros(size(X,1),1);

if numel(b)~=nb*(nc-1)+nb2
    error('Number of Covariates and Choice Alternatives Doesn''t Match Up!');
end

if isempty(Z)
    %normalizes BASEALT to zero
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        temp=X*b((k-1)*nb+1:k*nb);
        num(:,j)=exp(temp);
        dem=exp(temp)+dem;
        k=k+1;
    end
else
    %normalizes BASEALT to zero
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        temp=X*b((k-1)*nb+1:k*nb)+Z(:,:,j)*b2;
        num(:,j)=exp(temp);
        dem=exp(temp)+dem;
        k=k+1;
    end
end
num(:,baseAlt)=ones(size(X,1),1);
dem=dem+1;

P=num./(dem*ones(1,nc));