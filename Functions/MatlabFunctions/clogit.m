function [like,grad] = clogit(b,restrMat,Y,X,Z,baseAlt)
%CLOGIT estimates a conditional logit model
%   LIKE = CLOGIT(B,RESTRMAT,Y,X,Z,BASEALT)
%   estimates McFadden's choice modle, which is a special case of the
%   conditional logistic regression model. There are J choice
%   alternatives. Parameter restrictions are constructed in RESTRMAT.
%   
%   For estimation without restrictions, set RESTRMAT to be an empty matrix.
%   
%   Y is an N x 1 vector of integers 1 through J indicating which
%   alternative was chosen.
%   X is an N x K1 matrix of individual-specific covariates
%   Z is an N x K2 x J array of covariates that are alternative-specific.
%   BASEALT is the integer number of the category in Y that will be used as
%   the reference alternative. Alternative J is the default.
%   B is the parameter vector, with (J-1)*K1 + K2 elements
%   
%   CLOGIT can estimate one of three possible models:
%   1. Multinomial logit model: Z is empty
%   2. Conditional logit model: X is empty
%   3. Alternative-specific conditional logit model: X and Z both non-empty
%   
%   This function does *not* automatically include a column of ones in X
%   It also does *not* automatically drop NaNs
%   
%   Paramters are ordered as follows: {X parameters for alternative 1,
%   ..., X parameters for alternative J, Z parameters}
%   
%   Reference: McFadden, D. L. 1974. "Conditional Logit Analysis of
%   Qualitative Choice Behavior." in Frontiers in Econometrics, ed.
%   P. Zarembka. 105-142. New York: Academic Press.

% Copyright 2014 Jared Ashworth and Tyler Ransom, Duke University
% Special thanks to Vladi Slanchev and StataCorp's asclogit command
% Revision History:
%   July 28, 2014
%     Created
%==========================================================================

% error checking
assert((~isempty(X) || ~isempty(Z)) && ~isempty(Y),'You must supply data to the model');

N  = size(Y,1);
nb = size(X,2);
nb2= size(Z,2);
nc = numel(unique(Y));

assert(ndims(Y)==2 && size(Y,2)==1,'Y must be a vector');
assert(  min(Y)==1 && max(Y)==nc  ,'Y should contain integers numbered consecutively from 1 through J');
if ~isempty(X)
	assert(ndims(X)==2  ,'X must be a 2-dimensional matrix');
    assert(size(X,1)==N ,'The 1st dimension of X should equal the number of observations in Y');
end
if ~isempty(Z)
    assert(ndims(Z)==3  ,'Z must be a 3-dimensional array');
    assert(size(Z,1)==N ,'The 1st dimension of Z should equal the number of observations in Y');
    assert(size(Z,3)==nc,'The 3rd dimension of Z should equal the number of alternatives in Y');
end

% apply restrictions as defined in restrMat
if ~isempty(restrMat)
    b = applyRestr(restrMat,b);
end

if nargin==5 || isempty(baseAlt)
    baseAlt = nc;
end
b2   = b(nb*(nc-1)+1:nb*(nc-1)+nb2);
num  = zeros(N,1);
dem  = zeros(N,1);

if nb2>0 && nb>0
    %sets BASEALT to be the one that is normalized to zero
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        temp=X*b((k-1)*nb+1:k*nb)+(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        num=(Y==j).*temp+num;
        dem=exp(temp)+dem;
        k = k+1;
    end
    dem=dem+1;
    
    like=sum(log(dem)-num);
    
    % analytical gradient
    numg = zeros(nb2,1);
    demg = zeros(nb2,1);
    
    k = 1;
    grad = zeros(size(b));
    for j=setdiff(1:nc,baseAlt)
        temp=X*b((k-1)*nb+1:k*nb)+(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        grad((k-1)*nb+1:k*nb)=-X'*((Y==j)-exp(temp)./dem);
        k = k+1;
    end
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        temp = X*b((k-1)*nb+1:k*nb)+(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        numg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(Y==j)+numg;
        demg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(exp(temp)./dem)+demg;
        k = k+1;
    end
    grad(nb*(nc-1)+1:nb*(nc-1)+nb2)=numg-demg;
    if ~isempty(restrMat)
        grad = applyRestrGrad(restrMat,grad);
    end
elseif nb>0 && nb2==0
    %sets BASEALT to be the one that is normalized to zero
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        temp=X*b((k-1)*nb+1:k*nb);
        num=(Y==j).*temp+num;
        dem=exp(temp)+dem;
        k = k+1;
    end
    dem=dem+1;
    
    like=sum(log(dem)-num);
    
    % analytical gradient
    k = 1;
    grad = zeros(size(b));
    for j=setdiff(1:nc,baseAlt)
        temp=X*b((k-1)*nb+1:k*nb);
        grad((k-1)*nb+1:k*nb)=-X'*((Y==j)-exp(temp)./dem);
        k = k+1;
    end
    if ~isempty(restrMat)
        grad = applyRestrGrad(restrMat,grad);
    end
elseif nb==0 && nb2>0
    %sets BASEALT to be the one that is normalized to zero
    k = 1;
    for j=setdiff(1:nc,baseAlt)
        temp=(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        num=(Y==j).*temp+num;
        dem=exp(temp)+dem;
        k = k+1;
    end
    dem=dem+1;
    
    like=sum(log(dem)-num);
    
    % analytical gradient
    numg = zeros(nb2,1);
    demg = zeros(nb2,1);
    
    k = 1;
    grad = zeros(size(b));
    for j=setdiff(1:nc,baseAlt)
        temp =  (Z(:,:,j)-Z(:,:,baseAlt))*b2;
        numg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(Y==j)+numg;
        demg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(exp(temp)./dem)+demg;
        k = k+1;
    end
    grad(nb*(nc-1)+1:nb*(nc-1)+nb2)=numg-demg;
    if ~isempty(restrMat)
        grad = applyRestrGrad(restrMat,grad);
    end
end