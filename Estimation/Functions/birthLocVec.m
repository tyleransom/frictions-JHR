function out = birthLocVec(nloc,data,j)
% This function creates moving costs interacted with a quadratic in both age and distance
if isscalar(j)
   j=kron(ones(size(data)),j); 
end
data(isnan(data))=0;
disttricky = zeros(600000,1);
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
MC = [(kp~=jp) disttricky(1000*kp+jp,1) disttricky(1000*kp+jp,1).^2 a.*(kp~=jp) (a.^2).*(kp~=jp)];

end
