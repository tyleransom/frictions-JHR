function [like,grad,hess] = clogitHess(b,restrMat,Y,X,Z,baseAlt)
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
K1 = size(X,2);
K2= size(Z,2);
J = numel(unique(Y));

assert(isvector(Y),'Y must be a vector');
assert(  min(Y)==1 && max(Y)==J  ,'Y should contain integers numbered consecutively from 1 through J');
if ~isempty(X)
	assert(ismatrix(X)  ,'X must be a 2-dimensional matrix');
    assert(size(X,1)==N ,'The 1st dimension of X should equal the number of observations in Y');
end
if ~isempty(Z)
    assert(ndims(Z)==3  ,'Z must be a 3-dimensional array');
    assert(size(Z,1)==N ,'The 1st dimension of Z should equal the number of observations in Y');
    assert(size(Z,3)==J,'The 3rd dimension of Z should equal the number of alternatives in Y');
end

% apply restrictions as defined in restrMat
if ~isempty(restrMat)
    b = applyRestr(restrMat,b);
end

if nargin==5 || isempty(baseAlt)
    baseAlt = J;
end
b2   = b(K1*(J-1)+1:K1*(J-1)+K2);
num  = zeros(N,1);
dem  = zeros(N,1);
numer= ones(N,J);

if K2>0 && K1>0
    %sets BASEALT to be the one that is normalized to zero
    k = 1;
    for j=setdiff(1:J,baseAlt)
        temp=X*b((k-1)*K1+1:k*K1)+(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        num=(Y==j).*temp+num;
        dem=exp(temp)+dem;
        numer(:,j)=exp(temp);
        k = k+1;
    end
    dem=dem+1;
    P=numer./(dem*ones(1,J));
    
    like=sum(log(dem)-num);
    
    % analytical gradient
    numg = zeros(K2,1);
    demg = zeros(K2,1);
    
    k = 1;
    grad = zeros(size(b));
    for j=setdiff(1:J,baseAlt)
        temp=X*b((k-1)*K1+1:k*K1)+(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        grad((k-1)*K1+1:k*K1)=-X'*((Y==j)-exp(temp)./dem);
        k = k+1;
    end
    k = 1;
    for j=setdiff(1:J,baseAlt)
        temp = X*b((k-1)*K1+1:k*K1)+(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        numg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(Y==j)+numg;
        demg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(exp(temp)./dem)+demg;
        k = k+1;
    end
    grad(K1*(J-1)+1:K1*(J-1)+K2)=numg-demg;
    
    % analytical hessian
    k = 1;
    hess = zeros(length(b));
    for j=setdiff(1:J,baseAlt)
        kk = 1;
        for jj=setdiff(1:J,baseAlt)
            hess((k-1)*K1+1:k*K1,(kk-1)*K1+1:kk*K1)             = (X.*(P(:,j)*ones(1,K1)))'*(X.*((j==jj)-(P(:,jj)*ones(1,K1))));
            hess((k-1)*K1+1:k*K1,K1*(J-1)+1:K1*(J-1)+K2)        = (X.*(P(:,j)*ones(1,K1)))'*((Z(:,:,jj)-Z(:,:,baseAlt)).*((j==jj)-(P(:,jj)*ones(1,K2)))) + hess((k-1)*K1+1:k*K1,K1*(J-1)+1:K1*(J-1)+K2);
            hess(K1*(J-1)+1:K1*(J-1)+K2,(kk-1)*K1+1:kk*K1)      = ((Z(:,:,j)-Z(:,:,baseAlt)).*(P(:,j)*ones(1,K2)))'*(X.*((j==jj)-(P(:,jj)*ones(1,K1 )))) + hess(K1*(J-1)+1:K1*(J-1)+K2,(kk-1)*K1+1:kk*K1);
            hess(K1*(J-1)+1:K1*(J-1)+K2,K1*(J-1)+1:K1*(J-1)+K2) = ((Z(:,:,j)-Z(:,:,baseAlt)).*(P(:,j)*ones(1,K2)))'*((Z(:,:,jj)-Z(:,:,baseAlt)).*((j==jj)-(P(:,jj)*ones(1,K2)))) + hess(K1*(J-1)+1:K1*(J-1)+K2,K1*(J-1)+1:K1*(J-1)+K2);
            kk = kk+1;
        end
        k = k+1;
    end
    
    if ~isempty(restrMat)
        [grad,hess] = applyRestrGradHess(restrMat,grad,hess);
    end
elseif K1>0 && K2==0
    %sets BASEALT to be the one that is normalized to zero
    k = 1;
    for j=setdiff(1:J,baseAlt)
        temp=X*b((k-1)*K1+1:k*K1);
        num=(Y==j).*temp+num;
        dem=exp(temp)+dem;
        numer(:,j)=exp(temp);
        k = k+1;
    end
    dem=dem+1;
    P=numer./(dem*ones(1,J));
    
    like=sum(log(dem)-num);
    
    % analytical gradient
    k = 1;
    grad = zeros(size(b));
    for j=setdiff(1:J,baseAlt)
        temp=X*b((k-1)*K1+1:k*K1);
        grad((k-1)*K1+1:k*K1)=-X'*((Y==j)-exp(temp)./dem);
        k = k+1;
    end
    
    % analytical hessian
    k = 1;
    hess = zeros(length(b));
    for j=setdiff(1:J,baseAlt)
        kk = 1;
        for jj=setdiff(1:J,baseAlt)
            hess((k-1)*K1+1:k*K1,(kk-1)*K1+1:kk*K1) = (X.*(P(:,j)*ones(1,K1)))'*(X.*((j==jj)-(P(:,jj)*ones(1,K1))));
            kk = kk+1;
        end
        k = k+1;
    end
    
    if ~isempty(restrMat)
        [grad,hess] = applyRestrGradHess(restrMat,grad,hess);
    end
elseif K1==0 && K2>0
    %sets BASEALT to be the one that is normalized to zero
    k = 1;
    for j=setdiff(1:J,baseAlt)
        temp=(Z(:,:,j)-Z(:,:,baseAlt))*b2;
        num=(Y==j).*temp+num;
        dem=exp(temp)+dem;
        numer(:,j)=exp(temp);
        k = k+1;
    end
    dem=dem+1;
    P=numer./(dem*ones(1,J));
    
    like=sum(log(dem)-num);
    
    % analytical gradient
    numg = zeros(K2,1);
    demg = zeros(K2,1);
    
    k = 1;
    grad = zeros(size(b));
    for j=setdiff(1:J,baseAlt)
        temp =  (Z(:,:,j)-Z(:,:,baseAlt))*b2;
        numg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(Y==j)+numg;
        demg = -(Z(:,:,j)-Z(:,:,baseAlt))'*(exp(temp)./dem)+demg;
        k = k+1;
    end
    grad(K1*(J-1)+1:K1*(J-1)+K2)=numg-demg;
    
    % analytical hessian
    k = 1;
    hess = zeros(length(b));
    for j=setdiff(1:J,baseAlt)
        kk = 1;
        for jj=setdiff(1:J,baseAlt)
            hess = ((Z(:,:,j)-Z(:,:,baseAlt)).*(P(:,j)*ones(1,K2)))'*((Z(:,:,jj)-Z(:,:,baseAlt)).*((j==jj)-(P(:,jj)*ones(1,K2)))) + hess;
            kk = kk+1;
        end
        k = k+1;
    end
    
    if ~isempty(restrMat)
        [grad,hess] = applyRestrGradHess(restrMat,grad,hess);
    end
end
