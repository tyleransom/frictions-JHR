function y = wmean(x,wgt)
%WMEAN   Weighted average or mean value.
%   S = WMEAN(X,W) is the weighted mean value of the elements in X if X is a vector. 

y = sum(wgt.*x)./sum(wgt);

end
