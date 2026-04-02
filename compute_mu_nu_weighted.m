function [mu, nu, Pbar] = compute_mu_nu_weighted(A, lambda)
% Computes:
%   p_i        = ||a_i||^2 / ||A||_F^2
%   P_{lambda,i}    = (a_i a_i^T) / (||a_i||^2 + lambda)
%   Pbar       = sum_i p_i P_{lambda,i}
%   mu          = lambda_min^+(Pbar)
%   ν          = lambda_max( sum_i p_i (Pbar^{-1/2} P_{lambda,i} Pbar^{-1/2})^2 )

[m,n] = size(A);

% --- Sampling probabilities p_i = ||a_i||^2 / ||A||_F^2 ---
row_norm_sq = vecnorm(A,2,2).^2;      % m×1 vector of row norms squared
Z = sum(row_norm_sq);                 % ||A||_F^2
p = row_norm_sq / Z;                  % probabilities

% --- Build all P_{lambda,i} ---
P = cell(m,1);
for i = 1:m
    ai = A(i,:)';                     % n×1 column
    denom = norm(ai)^2 + lambda;
    P{i} = (ai * ai.') / denom;       % rank-1 PSD
end

% --- Compute Pbar = sum_i p_i P_i ---
Pbar = zeros(n);
for i = 1:m
    Pbar = Pbar + p(i) * P{i};
end

% --- Eigen-decomposition of Pbar ---
[V,D] = eig(Pbar);
d = diag(D);

% --- Pbar^{-1/2} ---
d_inv_sqrt = zeros(size(d));
positive = (d > 1e-14);
d_inv_sqrt(positive) = 1 ./ sqrt(d(positive));

Pbar_inv_sqrt = V * diag(d_inv_sqrt) * V';

% --- Compute second moment matrix M = sum_i p_i (T_i^2) ---
M = zeros(n);
for i = 1:m
    T = Pbar_inv_sqrt * P{i} * Pbar_inv_sqrt;
    M = M + p(i) * (T * T);
end

% --- mu = lambda_min^+(Pbar) ---
positive_eigs = d(positive);
if isempty(positive_eigs)
    mu = 0;
else
    mu = min(positive_eigs);
end

% --- nu = lambda_max(M) ---
nu = max(eig(M));

end
