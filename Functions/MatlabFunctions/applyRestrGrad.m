function [gRestr] = applyRestrGrad(restrMat,grad)
%APPLYRESTR applies restrictions to a model
%   G = APPLYRESTR(RESTRMAT,G) implements restrictions on the gradient
%   vector G of parameters according to the specifications found in
%   RESTRMAT. 
%   Two types of restrictions are supported:
%
%     Type 1  Restricting one parameter ("parmA") to equal a fixed value
%     Type 2  Restricting one parameter, parmA, to equal another ("parmB"),
%             potentially multiplied by some real number q.
%   
%   RESTRMAT follows a very specific format. It is an R-by-4 matrix, 
%   where R is the number of restrictions. The role of each of the four 
%   columns is as follows
%
% 	  Column 1  The index of parmA
% 	  Column 2  Either the fixed value for parmA the index of parmB
% 	  Column 3  Binary vector where 0 indciates a type 1 restriction (parmA
%               set equal to fixed value) and 1 indicates a type 2 
%               restriction (parmA set equal to parmB)
% 	  Column 4  Any real number q such that parmA = q*parmB. This should
%               be zero for a type 1 restriction.
%
%   APPLYRESTR does not allow for any combination of restrictions. If 
%   two parameters are to be restricted to the same fixed value, they
%   should both be type 1 restrictions rather than a type 1 restriction
%   and a type 2 restriction
%

% Copyright 2013 Jared Ashworth and Tyler Ransom, Duke University
%  .... oh, and Vladi, well, he was pretty important too. :)
% Last Revision: July 3, 2013
%==========================================================================
gRestr=grad;
sortrows(restrMat,1);
R = size(restrMat,1);
if R>0
	for r=1:R
        i = restrMat(r,1);
        h = restrMat(r,2);
        gRestr(i)=0;
		if restrMat(r,3)==1
			gRestr(h)=gRestr(h)+restrMat(r,4)*grad(i);
		end
	end
end


end