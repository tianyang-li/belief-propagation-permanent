function [s] = logsumexp(a)
% a - vector of real values
% s = log(sum(exp(a(i))))

max_a = max(a);

s = max_a;

s = s + log(sum(exp(a - max_a)));

end

