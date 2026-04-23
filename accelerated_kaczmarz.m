function fe = accelerated_kaczmarz(test_suite, xs, lambda, iter, refine)
% Implements:
%   x_k   = alpha*v_k + (1-alpha)*y_k
%   w_k   = regularized projection of x_k onto {x: Ax = b}
%   y_{k+1} = x_k - w_k
%   v_{k+1} = beta*v_k + (1-beta)*x_k - gamma*w_k
%   Optimal value of mu and nu are used
%   Probabilities are regularized
%   During iterative refinement, y and v are reset to 0
%
% Inputs:
%     test_suite     - matrices in a cell - can directly pass in test_suite
%                      generated from create_test_suite
%     xs             - vectors x in a cell - length matches len(test_suite)
%     lambda         - regularization parameter. If lambda < 0, then 
%                      lambda is set to norm(A, "fro")^2 / n
%     iter           - number of total iterations
%     refine         - refinement is performed whenever mod(cur_iter,
%                      refine) = 0. Can set refine > iter to skip
%                      refinement.
%   Outputs:
%     fe             - forward errors (scaled by |x|)

rng(124192032)
fe = zeros(numel(test_suite), 1000); % forward errors
%be = zeros(numel(test_suite), 1000); % backward errors
step = round(iter/1000);

tic
for index = 1:numel(test_suite)
    A = test_suite{index};
    x = xs{index};
    if lambda < 0
        lambda = norm(A, "fro")^2 / size(A, 2);
    end
    [mu, nu] = compute_mu_nu_weighted(A, lambda);
    beta = single(1 - sqrt(mu / nu));
    gamma = single(1 / sqrt(mu*nu));
    alpha = single(1 / (1 + gamma*nu));
    
    [m, n] = size(A);
    b = A*x;
    A_s = single(A);
    b_s = single(b);
    b2 = b_s;
    x_s = single(x);
    
    % Initialize
    xx = zeros(n, 1, 'single');
    v = zeros(n,1, 'single');
    y = zeros(n,1, 'single');
    
    squared_row_norms = sum(A_s.^2, 2) + lambda; % squared row norms, regularized
    prob = squared_row_norms / sum(squared_row_norms);  
    cdf = cumsum(prob);              % cumulative probability
    
    count = 1;
    for i = 1:iter
        r = single(rand);
        k = find(cdf >= r, 1);  
        if isempty(k)
            k = m;
        end

        xx2 = alpha * v + (1 - alpha) * y;
        ak = A_s(k,:)';                    % n×1
        w = (ak*(ak' * xx2 - b2(k))) / (squared_row_norms(k));
        y_new = xx2 - w;
        v_new = beta * v + (1 - beta) * xx2 - gamma * w;
        
        if mod(i,step) == 0
            fe(index, count) = norm(double(xx + xx2) - x) / norm(x);
            %be(index, count) = norm(A * double(xx + xx2) - b) / (norm(A) * norm(double(xx + xx2)));
            count = count + 1;
            if mod(i, refine) == 0
                xx = xx + xx2;    
                y_new = zeros(n, 1, 'single');
                v_new = zeros(n, 1, 'single');
                xx2 = zeros(n, 1, 'single');
                b2 = b_s - A_s*xx;
                disp(norm(double(xx + xx2) - x) / norm(x))
            end
        end
        y = y_new;
        v = v_new;
    end

    toc
    save('data/prelim_accel.mat', 'fe', 'index')
end
end
