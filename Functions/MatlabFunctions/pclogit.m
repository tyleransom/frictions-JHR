function P=pclogit(b,Y,X,Z,baseAlt)
%PCLOGIT1 returns fitted probabilities from a conditional logit model
%   P = PCLOGIT1(B,Y,X,Z,NB,NB2,NC,BASEALT) 
%   returns fitted probabilites from a parameter vector estimated by the 
%   conditional logistic regression model. There are J choice
%   alternatives.
%   
%   P is an N x J matrix of fitted probabilities
%   B is the parameter vector, with (J-1)*K1 + K2 elements
%   Y is an N x 1 vector of integers 1 through J indicating which
%   alternative was chosen. 
%   X is an N x K1 matrix of individual-specific covariates.
%   Z is an N x K2 x J array of covariates that are alternative-specific.
%   BASEALT is the integer number of the category in Y that will be used
%   as the reference alternative. Alternative J is the default.
%   
%   This function does *not* automatically include a column of ones in X.
%   It also does *not* automatically drop NaNs
%   
%   Parameters should be ordered as follows: {X parameters for alternative 1,
%   ..., X parameters for alternative J, Z parameters}

% Copyright 2014 Jared Ashworth and Tyler Ransom, Duke University
% Special thanks to Vladi Slanchev and StataCorp's asclogit command
% Revision History: 
%   July 29, 2014
%     Created
%==========================================================================

% error checking
assert((~isempty(X) || ~isempty(Z)) && ~isempty(Y) && ~isempty(b),'You must supply parameters and data to the model');

N  = size(Y,1);
nb = size(X,2);
nb2= size(Z,2);
nc = numel(unique(Y));

assert(ndims(Y)==2 && size(Y,2)==1 ,'Y must be a vector');
assert(  min(Y)==1 && max(Y)==nc   ,'Y should contain integers numbered consecutively from 1 through J');
assert(ndims(b)==2 && size(b,2)==1 ,'b must be a vector');
if ~isempty(X)
	assert(ndims(X)==2  ,'X must be a 2-dimensional matrix');
	assert(size(X,1)==N ,'The 1st dimension of X should equal the number of observations in Y');
end
if ~isempty(Z)
	assert(ndims(Z)==3  ,'Z must be a 3-dimensional array');
	assert(size(Z,1)==N ,'The 1st dimension of Z should equal the number of observations in Y');
	assert(size(Z,3)==nc,'The 3rd dimension of Z should equal the number of alternatives in Y');
end
assert(size(b,1)==(nc-1)*nb+nb2,'b should be of dimension ((J-1)*K1 + K2) x 1');

if nargin==4 || isempty(baseAlt)
	baseAlt = nc;
end
b2   = b(nb*(nc-1)+1:nb*(nc-1)+nb2);
num  = zeros(N,nc);
dem  = zeros(N,1);

if nb2>0 && nb>0
	% sets BASEALT to be the alternative that is normalized to zero
	k = 1;
	for j=setdiff(1:nc,baseAlt)
		temp=X*b((k-1)*nb+1:k*nb)+(Z(:,:,j)-Z(:,:,baseAlt))*b2;
		num(:,j)=exp(temp);
		dem=exp(temp)+dem;
		k = k+1;
	end
elseif nb>0 && nb2==0
	% sets BASEALT to be the alternative that is normalized to zero
	k = 1;
	for j=setdiff(1:nc,baseAlt)
		temp=X*b((k-1)*nb+1:k*nb);
		num(:,j)=exp(temp);
		dem=exp(temp)+dem;
		k = k+1;
	end
elseif nb==0 && nb2>0
	% sets BASEALT to be the alternative that is normalized to zero
	k = 1;
	for j=setdiff(1:nc,baseAlt)
		temp=(Z(:,:,j)-Z(:,:,baseAlt))*b2;
		num(:,j)=exp(temp);
		dem=exp(temp)+dem;
		k = k+1;
	end
end
num(:,baseAlt)=ones(N,1);
dem=dem+1;

P=num./(dem*ones(1,nc));
