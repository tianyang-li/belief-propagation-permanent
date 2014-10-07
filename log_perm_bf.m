function [log_Z] = log_perm_bf(A, B)
%
% A - matrix m by n + 1 (log weights)
% B - matrix n by m + 1 (log weights)
%
% A has m elements that may or may not be assigned to one in B
% B has n elelemts that may or may not be assigned to one in A
%

[m n] = size(A);
n = n - 1;

if ~isequal([n, m + 1], size(B))
	throw(MException('LogPermBF:unmatchedDim', 'dimensionnot matched'));
end

log_Z = sum(A(:, n + 1)) + sum(B(:, m + 1));

maxMatch = min(m, n);

for i = 1:maxMatch
	invPerm = zeros(1, i);
	aC = nchoosek(1:m, i);
	bC = nchoosek(1:n, i);
	for a = 1:size(aC, 1)
		for b = 1:size(bC, 1)
			p = perms(1:i);
			tmp1 = sum(A(setdiff(1:m, aC(a, :)), n + 1)) ...
				+ sum(B(setdiff(1:n, bC(b, :)), m + 1));
			for j = 1:size(p, 1);
				invPerm(p(j, :)) = 1:i;
				tmp2 = sum(A(sub2ind(size(A), aC(a, :), bC(b, p(j, :))))) ...
					+ sum(B(sub2ind(size(B), bC(b, :), aC(a, invPerm)))) ...
					+ tmp1;
				tmp3 = logsumexp([log_Z tmp2]);
				log_Z = tmp3;
			end
		end
	end
end

end



