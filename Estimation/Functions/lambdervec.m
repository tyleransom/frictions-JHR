function lambda = lambdervec(nloc,j,t,e,l)
% This function plugs in predicted frictions into Z matrix
% assert(isequal(sum(Tmat,2),ones(size(Tmat,1),1)),'Error in Tmat.')
if j>nloc
	lambda = zeros(size(e));
else
	lambda = l(j,t)';
end

end