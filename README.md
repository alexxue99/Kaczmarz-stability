# Kaczmarz-stability
Tests stability of (accelerated) Kaczmarz, with and without iterative refinement, on a test suite.

## Example code flow
```matlab
run('create_test_suite.m') % modify the condition number and size of matrix inside the file itself

for i = 1:5
x = randn(n, 1); x = x / norm(x);
xs{i} = x;
end

% Optional: if you want to focus on just one matrix in the test suite
% test_suite = {A1};
% xs = {xs{1}}; 

iter = 1e9;
[fe, ferr] = kaczmarz(test_suite, xs, iter, iter + 1); % no refinement
feR = kaczmarz(test_suite, xs, iter, iter/2); % refine at halfway mark
feA = accelerated_kaczmarz(test_suite, xs, 1e-6, iter, iter + 1); % no refinement
feAR = accelerated_kaczmarz(test_suite, xs, 1e-6, iter, iter/10); % iterate 10 times

y = 1e6:1e6:1e9;

graph_test_suite(n, cond_num, ferr, y, false, fe, feR, feA, feAR); % replace false with true if you want to save figure

% Optional: save into .mat file
% save('data\' + n + '_' + cond_num_str + '.mat', 'fe', 'feR', 'feA', 'feAR', 'ferr', 'xs')