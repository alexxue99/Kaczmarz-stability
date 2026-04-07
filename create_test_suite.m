rng(123);
cond_num = 1e5;
n = 500;
A1 = make_test_matrix_cond(n, n, 'poly', [], cond_num);
A2 = make_test_matrix_cond(n, n, 'exp', [], cond_num);
%A3 = make_test_matrix_cond(n, n, 'cluster', .1, cond_num);
%A4 = make_test_matrix_cond(n, n, 'cluster', .2, cond_num);
%A3 = make_test_matrix_cond(n, n, 'lowrank', .1, cond_num);
%A6 = make_test_matrix_cond(n, n, 'lowrank', .2, cond_num);
A3 = make_test_matrix_cond(n, n, 'lowrank', .9, cond_num); %highrank
%A5 = make_test_matrix_cond(n, n, 'linear', [], cond_num);
A4 = make_test_matrix_cond(n, n, 'harmonic', [], cond_num);

test_suite = {A1, A2, A3, A4};

function A = make_test_matrix_cond(m, n, type, param, cond_num)
%MAKE_TEST_MATRIX_COND  Generate m-by-n test matrix whose decay law yields cond_num
%
%   A = MAKE_TEST_MATRIX_COND(m, n, type, param, cond_num)
%
%   Inputs:
%     m, n     - matrix sizes
%     type     - 'poly' | 'exp' | 'cluster' | 'lowrank' | 'linear' |
%     'harmonic'
%     param    - shape parameter OR []/NaN to auto-compute param from cond_num:
%                - poly:  alpha (if given). If param==[] then alpha = log(cond_num)/log(r).
%                - exp:   beta  (if given). If param==[] then beta  = log(cond_num)/(r-1).
%                - cluster: frac_high (fraction of large singular values). If param==[] uses 0.1.
%                - lowrank: rank_frac (fraction of non-negligible singulars). If param==[] uses 0.1.
%                - linear: linear decay from 1 to 1/cond_num. Param unused.
%     cond_num - desired condition number (>=1)
%
%
%   Output:
%     A - m-by-n matrix with singular values in descending order and condition number cond_num
%
%   Example:
%     A = make_test_matrix_cond(300,200,'poly',[],1e6); % auto alpha to get cond=1e6

    if nargin < 5
        error('Usage: make_test_matrix_cond(m,n,type,param,cond_num)');
    end

    if cond_num < 1
        error('cond_num must be >= 1');
    end

    r = min(m, n);
    k = (1:r)';

    switch lower(type)
        case 'poly'
            % s_k = k^{-alpha}
            if isempty(param) || isnan(param)
                % solve alpha from cond_num = 1 / r^{-alpha} = r^{alpha}
                % => alpha = log(cond_num)/log(r)
                alpha = log(cond_num) / log(r);
            else
                alpha = param;
            end
            s = k .^ (-alpha);

            left = 0;
            right = alpha;

            while true
                m = (left + right) / 2;
                s = k.^(-m);
                c = norm(s) / min(s);

                if abs(c - cond_num) < 1e-3
                    break
                elseif c > cond_num
                    right = m;
                else
                    left = m;
                end
            end

        case 'exp'
            % s_k = exp(-beta*(k-1)); raw cond = exp(beta*(r-1)) => beta = log(cond)/ (r-1)
            if isempty(param) || isnan(param)
                if r == 1
                    beta = 0;
                else
                    beta = log(cond_num) / (r - 1);
                end
            else
                beta = param;
            end
            s = exp(-beta * (k - 1));
            left = 0;
            right = beta;

            while true
                m = (left + right) / 2;
                s = exp(-m*(k-1));
                c = norm(s) / min(s);

                if abs(c - cond_num) < 1e-3
                    break
                elseif c > cond_num
                    right = m;
                else
                    left = m;
                end
            end

        case 'cluster'
            % param = fraction of "large" singular values (e.g. 0.1)
            frac = param;
            n_high = max(1, round(frac * r));
            % choose three clusters: high, mid, low. Set low plateau so that cond holds.
            % Let high = 1, low = 1/cond_num. mid = geometric mean between.
            n_mid = max(0, round(0.25 * (r - n_high)));
            n_low = r - n_high - n_mid;
            high_vals = cond_num * ones(n_high, 1);
            if n_mid > 0
                mid_val = sqrt(cond_num); % geometric mean
                mid_vals = mid_val * ones(n_mid, 1);
            else
                mid_vals = [];
            end
            low_vals = ones(n_low, 1);
            s = [high_vals; mid_vals; low_vals];
            left = 1;
            right = cond_num;

            while true
                m = (left + right) / 2;
                high_vals = m * ones(n_high, 1);
                mid_vals = sqrt(m) * ones(n_mid, 1);
                s = [high_vals; mid_vals; low_vals];
                c = norm(s) / min(s);

                if abs(c - cond_num) < 1e-3
                    break
                elseif c > cond_num
                    right = m;
                else
                    left = m;
                end
            end

        case 'lowrank'
            r_eff = max(1, round(param * r));
            % interpolate first r_eff values from 1 down to 1/cond_num (so cond holds)
            if r_eff == 1
                leading = 1;
            else
                leading = linspace(1, 1/cond_num, r_eff).';
            end
            % remaining singular values are set to the tiny tail = 1/cond_num
            tail_count = r - r_eff;
            tail_vals = (1/cond_num) * ones(tail_count, 1);
            s = [leading; tail_vals];

            left = 1;
            right = cond_num;

            while true
                m = (left + right) / 2;
                leading = linspace(1, 1/m, r_eff).';
                tail_vals = 1/m * ones(tail_count, 1);
                s = [leading; tail_vals];
                c = norm(s) / min(s);

                if abs(c - cond_num) < 1e-3
                    break
                elseif c > cond_num
                    right = m;
                else
                    left = m;
                end
            end
        
        case 'linear'
            s = linspace(1, 1/cond_num, r);
            left = 1;
            right = cond_num;

            while true
                m = (left + right) / 2;
                s = linspace(1, 1/m, r);
                c = norm(s) / min(s);

                if abs(c - cond_num) < 1e-3
                    break
                elseif c > cond_num
                    right = m;
                else
                    left = m;
                end
            end

        case 'harmonic'
            s = linspace(1, cond_num, r);
            for i = 1:numel(s)
                s(i) = 1 / s(i);
            end
            left = 1;
            right = cond_num;

            while true
                m = (left + right) / 2;
                s = linspace(1, m, r);
                for i = 1:numel(s)
                    s(i) = 1 / s(i);
                end
                c = norm(s) / min(s);

                if abs(c - cond_num) < 1e-3
                    break
                elseif c > cond_num
                    right = m;
                else
                    left = m;
                end
            end

        otherwise
            error('Unknown type ''%s''', type);
    end

    % construct random orthonormal U,V (economy QR)
    [U, ~] = qr(randn(n, r), 0);
    [V, ~] = qr(randn(n, r), 0);

    A = U * diag(s) * V';
end