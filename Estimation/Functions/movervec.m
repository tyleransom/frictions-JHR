function MC = movervec(nloc,k,j,dist,a,prevEmp,prevUnemp)
% This function creates moving costs interacted with a quadratic in both age and distance
if isscalar(j)
   j=kron(ones(size(a)),j); 
end
if isscalar(k)
   k=kron(ones(size(a)),k); 
end
if ~isequal(size(k),size(a))
   error('LY and age need to have same dimension');
end
if ~isequal(size(k),size(j))
   error('LY and current choice need to have same dimension');
end
a(isnan(a))=0;
dist = dist./1e3;
disttricky = zeros(110110,1);
for ii=1:2*nloc
    for jjj=1:2*nloc
        disttricky(1000*ii+jjj)=dist(ii,jjj);
    end;
end
k(k==0)=1;
j(j==0)=1;
% switch happens only if lp~=l
kp = k-nloc*(k>nloc);
jp = j-nloc*(j>nloc);
MC = [(kp~=jp) disttricky(1000*kp+jp,1) disttricky(1000*kp+jp,1).^2 a.*(kp~=jp) (a.^2).*(kp~=jp) prevEmp.*(kp~=jp) prevUnemp.*(kp~=jp)];

end
