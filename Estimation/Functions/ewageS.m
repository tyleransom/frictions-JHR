function w = ewageS(nloc,j,t,e,type,b)
% This function computes expected wages given state variables
locer = zeros(1,nloc-1);
timer = zeros(1,9);
if nloc==55
    loctimer = zeros(1,(nloc-3)*9);
else
    loctimer = zeros(1,(nloc-1)*9);
end
if j>1
	locer(j-1) = 1;
end
if t>1
    timer(t-1) = 1;
end
if nloc==55
    if j>1 && t>1 && j<nloc-1
        loctimer((j-2)*size(timer,2)+t-1) = 1;
    end
else
    if j>1 && t>1
        loctimer((j-2)*size(timer,2)+t-1) = 1;
    end
end
if j>nloc
	w = 0;
else
	w = [1 locer timer loctimer e e.^2./100 type==1]*b;
end

end
