function [like,grad]=normalMLE(b,restrMat,Y,X,d)

% error checking
assert(size(X,1)==size(Y,1),'X and Y must be the same length');
if nargin==4
	d = ones(size(Y));
elseif nargin==5 && isempty(d)
	d = ones(size(Y));
end
J = numel(unique(d));
assert(min(d)==1 && max(d)==J,'d should contain integers numbered consecutively from 1 through J');

% apply restrictions as defined in restrMat
if ~isempty(restrMat)
	b = applyRestr(restrMat,b);
end

% slice parameter vector
beta      = b(1:end-J);
wagesigma = b(end-(J-1):end);
n         = length(Y);

% log likelihood
likemat = zeros(n,J);
dmat    = zeros(n,J);
for j=1:J
	dmat(:,j)    = d==j;
	likemat(:,j) = -.5*(log(2*pi)+log(wagesigma(j)^2)+((Y-X*beta)./wagesigma(j)).^2);
end
like = -sum(sum(dmat.*likemat));

% analytical gradient
grad = zeros(size(b));
for j=1:J
	grad(1:end-J) = -X'*((d==j).*(Y-X*beta)./(wagesigma(j).^2)) + grad(1:end-J);
end
for j=1:J
	k=length(b)-(J-1)+j-1;
	temp = 1./wagesigma(j)-((Y-X*beta).^2)./(wagesigma(j).^3);
	grad(k) = sum((d==j).*temp);
end

end
