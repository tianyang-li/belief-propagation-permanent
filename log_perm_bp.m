function [log_Z] = log_perm_bp(A, B)
%
% A - matrix m by n + 1 of log weights
% B - matrix n by m + 1 of log weight
%
% A has m elements that may or may not be assigned to one in B
% B has n elelemts that may or may not be assigned to one in A
%

[m n] = size(A);
n = n - 1;

if ~isequal([n, m + 1], size(B))
	throw(MException('LogPermBP:unmatchedDim', 'dimension not matched'));
end

% TODO: tune parameters?
MAX_ITER = 1000;
ERR = 1e-5;

mAB = zeros(m, n); % A -> B
tmpA = zeros(m, n);

mBA = zeros(n, m); % B -> A
tmpB = zeros(n, m);

iter = 0;

% temp variable for sum product for A
A_sp = zeros(m, n); 
% temp variable for product of messages with factor for A
A_pf = zeros(m, n); 

% temp variable for sum product for B
B_sp = zeros(n, m); 
% temp variable for product of messages with factor for B
B_pf = zeros(n, m); 

m_ones = ones(1, m);
n_ones = ones(1, n);

while 1
	tmpA = mAB;
	tmpB = mBA;
	
	A_pf = mBA.' + A(:, 1:n);
	for i = 1:m
		for j = 1:n
			A_sp(i, j) = logsumexp([A(i, n + 1) ...
				A_pf(i, [1:(j - 1) (j + 1):n])]);
		end
	end
	mAB = A(:, 1:n) - A_sp;
	
	B_pf = tmpA.' + B(:, 1:m);
	for i = 1:n
		for j = 1:m
			B_sp(i, j) = logsumexp([B(i, m + 1) ...
				B_pf(i, [1:(j - 1) (j + 1):m])]);
		end
	end
	mBA = B(:, 1:m) - B_sp;
	
	tmpA = mAB - tmpA;
	tmpB = mBA - tmpB;
	err = max([max(abs(tmpA(:))), max(abs(tmpB(:)))]);
	
	iter = iter + 1;	
	if iter > MAX_ITER || err < ERR
		break;
	end
end

log_b_A = zeros(m, n + 1);
log_b_A(:, 1:n) = A(:, 1:n) + mBA.';
log_b_A(:, n + 1) = A(:, n + 1);

tmp_A_1 = log_b_A;

for i = 1:m
	log_sum = logsumexp(log_b_A(i, :));
	log_b_A(i, :) = log_b_A(i, :) - log_sum;
end

log_b_B = zeros(n, m + 1);
log_b_B(:, 1:m) = B(:, 1:m) + mAB.';
log_b_B(:, m + 1) = B(:, m + 1);

tmp_B_1 = log_b_B;

for i = 1:n
	log_sum = logsumexp(log_b_B(i, :));
	log_b_B(i, :) = log_b_B(i, :) - log_sum;
end

log_Z = (n - 1) * sum(sum(log_b_A .* exp(log_b_A), 2)) ...
	+ (m - 1) * sum(sum(log_b_B .* exp(log_b_B), 2)) ...
	+ sum(sum(exp(log_b_A) .* A, 2)) ...
	+ sum(sum(exp(log_b_B) .* B, 2));

log_b_AB = zeros(m + 1, n + 1);

for i = 1:m
	for j = 1:n
		log_b_AB = tmp_B_1(ones(n + 1, 1) * j, :).' ...
			+ tmp_A_1(ones(m + 1, 1) * i, :);

		tmp_AB_1 = log_b_AB(1:(i - 1), 1:(j - 1));
		tmp_AB_2 = log_b_AB(1:(i - 1), (j + 1):(n + 1));
		tmp_AB_3 = log_b_AB((i + 1):(m + 1), 1:(j - 1));
		tmp_AB_4 = log_b_AB((i + 1):(m + 1), (j + 1):(n + 1));
		tmp_AB_ji = A(i, j) + B(j, i);
		
		tmp_AB = [tmp_AB_1(:); tmp_AB_2(:); ...
			tmp_AB_3(:); tmp_AB_4(:); tmp_AB_ji];
		
		log_sum = logsumexp(tmp_AB);
		
		tmp_AB = tmp_AB - log_sum;
		
		log_Z = log_Z - sum(tmp_AB .* exp(tmp_AB));
	end
end

end



