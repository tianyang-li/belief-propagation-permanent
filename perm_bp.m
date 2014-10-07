function [Z] = perm_bp(A, B)
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
	throw(MException('PermBP:unmatchedDim', 'dimensionnot matched'));
end

% TODO: tune parameters?
MAX_ITER = 1000;
ERR = 1e-5;

mAB = ones(m, n); % A -> B
tmpA = ones(m, n);

mBA = ones(n, m); % B -> A
tmpB = ones(n, m);

iter = 0;

A_sp = zeros(m, 1); % temp variable for sum product for A
A_pf = zeros(m, n); % temp variable for product of messages with factor for A

B_sp = zeros(n, 1); % temp variable for sum product for B
B_pf = zeros(n, m); % temp variable for product of messages with factor for B

m_ones = ones(1, m);
n_ones = ones(1, n);

while 1
	tmpA = mAB;
	tmpB = mBA;
	
	A_pf = mBA.' .* A(:, 1:n);
	A_sp = A(:, n + 1) + sum(A_pf, 2);
	mAB = A(:, 1:n) ./ (A_sp(:, n_ones) - A_pf);
	
	B_pf = tmpA.' .* B(:, 1:m);
	B_sp = B(:, m + 1) + sum(B_pf, 2);
	mBA = B(:, 1:m) ./ (B_sp(:, m_ones) - B_pf);
	
	tmpA = mAB - tmpA;
	tmpB = mBA - tmpB;
	err = max([max(abs(tmpA(:))), max(abs(tmpB(:)))]);
	
	iter = iter + 1;	
	if iter > MAX_ITER || err < ERR
		break;
	end
end

Z = 0;

% factor * message
A_fm = [mBA.' .* A(:, 1:n), A(:, n + 1)];
A_fm_sum = sum(A_fm, 2);

B_fm = [mAB.' .* B(:, 1:m), B(:, m + 1)];
B_fm_sum = sum(B_fm, 2);

b_A = A_fm ./ A_fm_sum(:, [n_ones 1]);
b_B = B_fm ./ B_fm_sum(:, [m_ones 1]);

Z = Z + (n - 1) * sum(sum(x_log_x(b_A), 2));
Z = Z + (m - 1) * sum(sum(x_log_x(b_B), 2));

Z = Z + sum(sum(b_A .* log(A), 2));
Z = Z + sum(sum(b_B .* log(B), 2));

for a = 1:m
	for b = 1:n
		b_ab = B_fm(b, :).' * A_fm(a, :);
		b_ab(a, :) = zeros(1, n + 1);
		b_ab(:, b) = zeros(m + 1, 1);
		b_ab(a, b) = A(a, b) * B(b, a);
		Z_ab = sum(b_ab(:));
		b_ab = b_ab / Z_ab;
		b_ab = x_log_x(b_ab);
		Z = Z - sum(b_ab(:));
	end
end

Z = exp(Z);

end

function [y] = x_log_x(x)

y = zeros(size(x));
y(x > 0) = x(x > 0) .* log(x(x > 0));

end


