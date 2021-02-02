function [prior,PType,jointlike] = typeprob(prior,base)

N = size(base,1);
S = size(base,2);

for s=1:S
	PType(:,s) = prior(s)*base(:,s)./(base*prior');
end

prior = mean(PType(base(:,1)~=1,:));
jointlike = sum(log(base*prior'));

end
