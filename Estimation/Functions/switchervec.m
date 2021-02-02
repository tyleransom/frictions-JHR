function SC = switchervec(nloc,k,j,a)
% This function creates switching costs interacted with a quadratic in age
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
% switch happens only if kp~=jp (LY(:)> nloc)
kp = k<=nloc;
jp = j<=nloc;
SC = [(kp~=jp) a.*(kp~=jp) (a.^2).*(kp~=jp)];

end