function [Z] = perm_bf(A, B)
%
% A - matrix m by n + 1
% B - matrix n by m + 1
%
% A has m elements that may or may not be assigned to one in B
% B has n elelemts that may or may not be assigned to one in A
%

[m n] = size(A);
n = n - 1;

if ~isequal([n, m + 1], size(B))
	throw(MException('PermBF:unmatchedDim', 'dimensionnot matched'));
end

Z = prod(A(:, n + 1)) * prod(B(:, m + 1));

maxMatch = min(m, n);

for i = 1:maxMatch
	invPerm = zeros(1, i);
	aC = nchoosek(1:m, i);
	bC = nchoosek(1:n, i);
	for a = 1:size(aC, 1)
		for b = 1:size(bC, 1)
			p = perms(1:i);
			tmp1 = prod(A(setdiff(1:m, aC(a, :)), n + 1)) ...
				* prod(B(setdiff(1:n, bC(b, :)), m + 1));
			for j = 1:size(p, 1);
				invPerm(p(j, :)) = 1:i;
				Z = Z + prod(A(sub2ind(size(A), aC(a, :), bC(b, p(j, :))))) ...
					* prod(B(sub2ind(size(B), bC(b, :), aC(a, invPerm)))) ...
					* tmp1;
			end
		end
	end
end

end



