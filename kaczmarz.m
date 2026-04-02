function [fe, ferr] = kaczmarz(test_suite, xs, iter, refine)
%   Inputs:
%     test_suite     - matrices in a cell - can directly pass in test_suite
%                      generated from create_test_suite
%     xs             - vectors x in a cell - length matches len(test_suite)
%     iter           - number of total iterations
%     refine         - refinement is performed whenever mod(cur_iter,
%                      refine) = 0. Can set refine > iter to skip
%                      refinement.
%   Outputs:
%     fe             - forward errors (scaled by |x|)
%     ferr           - baseline forward error from A\b (scaled by |x|)
rng(124192032)
num = length(test_suite);
precision = 'single';
step = iter/1000;
fe = zeros(num, iter/step);
be = zeros(num, iter/step);
ferr = zeros(num, 1);
%berr = zeros(num, 1);


for index = 1:numel(test_suite)
    disp("working on matrix " + index)
    A = test_suite{index};
    n = size(A, 1);
    x = xs{index};
    
    b = A*x;

    A_s = single(A); x_s = single(x);
    b_s = single(b);
    ferr(index) = norm(A_s\b_s - x_s) / norm(x_s);
    %berr(index) = norm(A_s * (A_s\b_s) - b_s) / (norm(A_s) * norm(x_s));
    % Optional normalize rows of A
    % rownorms = sqrt(sum(A_s.^2, 2));
    % A_s = bsxfun(@rdivide, A_s, rownorms);
    % b_s = bsxfun(@rdivide, b_s, rownorms);

    count = 1;
    
    xx = zeros(n,1, precision);   % current iterate
    xx2 = zeros(n,1, precision);   % current iterate
    b2 = b_s;
    squared_row_norms = sum(A_s.^2, 2);       % squared row norms
    prob = squared_row_norms / sum(squared_row_norms);  
    cdf = cumsum(prob);              % cumulative probability

    for i = 1:iter
        r = rand;
        k = find(cdf >= r, 1);  
        if isempty(k)
            k = n;
        end
        % --- Kaczmarz step ---
        update = (((b2(k) - A_s(k,:)*xx2)) / squared_row_norms(k)) * A_s(k,:)';
        xx2 = xx2 + update;
     
        % --- Update history ---
        if mod(i,step) == 0
            fe(index, count) = norm(double(xx + xx2) - x) / norm(x);
            be(index, count) = norm(A * double(xx + xx2) - b) / (norm(A) * norm(double(xx + xx2)));
            if mod(i, refine) == 0
                xx = xx + xx2;  
                xx2 = zeros(n, 1, precision);   
                b2 = b_s - A_s*xx;
            end
        count = count + 1;
        end
    end
end
end
%%
function A = rand_with_evals(eigvals)
    n = length(eigvals);
    [Q, ~] = qr(randn(n));  
    Lambda = diag(eigvals);
    A = Q * Lambda / Q;
end
