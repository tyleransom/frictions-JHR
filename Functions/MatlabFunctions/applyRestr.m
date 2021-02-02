function [bRestr,invHrestr] = applyRestr(restrMat,b,H)
%APPLYRESTR applies restrictions to a model
%   B = APPLYRESTR(RESTRMAT,B) implements restrictions on the vector B 
%   of parameters according to the specifications found in RESTRMAT. 
%   Two types of restrictions are supported:
%
%     Type 1  Restricting one parameter ("parmA") to equal a fixed value
%     Type 2  Restricting one parameter, parmA, to equal another ("parmB"),
%             potentially multiplied by some real number q and addd to
%             some constant m, e.g. parmA = m + q*parmB.
%   
%   RESTRMAT follows a very specific format. It is an R-by-5 matrix, 
%   where R is the number of restrictions. The role of each of the four 
%   columns is as follows
%
% 	  Column 1  The index of parmA
% 	  Column 2  The index of parmB (zero if type 1 restriction)
% 	  Column 3  Binary vector where 0 indciates a type 1 restriction (parmA
%               set equal to fixed value) and 1 indicates a type 2 
%               restriction (parmA set equal to parmB)
% 	  Column 4  If a type 1 restriction, 0. If a type 2 restriction, any 
%               real number q such that parmA = q*parmB.
% 	  Column 5  If a type 1 restriction, the fixed value. If a type 2
%               restriction, any real number m such that parmA = m+q*parmB.
%
%   APPLYRESTR does not allow for any combination of restrictions. If 
%   two parameters are to be restricted to the same fixed value, they
%   should both be type 1 restrictions rather than a type 1 restriction
%   and a type 2 restriction
%
%   [B,INVH] = APPLYRESTR(RESTRMAT,B,H) takes as input a hessian matrix H
%   (typically from an optimization routine) and returns an inverted
%   hessian INVH where restrictions have been applied. A type 1 restriction
%   results in a row and a column of zeroes in the hessian at the index for
%   that paramter. A type 2 restriction duplicates
%   the rows and columns for parmA and parmB. This implies that the 
%   covariance between parmA and parmB is set equal to their variance.
%   Moreover, the covariance of parmA and parmB with other parameters is
%   restricted to be equal.
%
%   hessian with H
% The size of your parameter vector never need change. However, inverting
%  a hessian with restrictions will create undesirable outcomes. Thus,
%  a 'correct' inverse hessian matrix is created.

% This function allows easy application of restrictions to fminunc and
%  other optimaztion routines.


% Copyright 2013 Jared Ashworth and Tyler Ransom, Duke University
%  .... oh, and Vladi, well, he was pretty important too. :)
% Last Revision: July 2, 2013
%==========================================================================
bRestr=b;
sortrows(restrMat,1);
R = size(restrMat,1);
if R>0
	for r=1:R
		if restrMat(r,3)==0
			bRestr(restrMat(r,1))=restrMat(r,5);
		elseif restrMat(r,3)==1
			bRestr(restrMat(r,1))=restrMat(r,5)+restrMat(r,4)*bRestr(restrMat(r,2));
		end
	end
end

if nargin==3
	if R>0
		for r=R:-1:1
			H(restrMat(r,1),:)=[];
			H(:,restrMat(r,1))=[];
		end
	end
	invH = full(H)\eye(size(H));
	if R>0
        for r=1:R
			invH= [ invH(1:restrMat(r,1)-1,:); zeros(1,size(invH,1)); invH(restrMat(r,1):end,:)];
			invH= [ invH(:,1:restrMat(r,1)-1)  zeros(size(invH,1),1)  invH(:,restrMat(r,1):end)];
        end
        for r=1:R
            if restrMat(r,3)==1
				invH(restrMat(r,1),:)=restrMat(r,4)*invH(restrMat(r,2),:);
				invH(:,restrMat(r,1))=restrMat(r,4)*invH(:,restrMat(r,2));
            end
        end
	end
    invHrestr = invH;
end

end