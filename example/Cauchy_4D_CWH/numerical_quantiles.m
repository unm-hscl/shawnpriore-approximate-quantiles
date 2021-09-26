%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantile approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up
by = 5e-6; % step size
approx_to = .9995;
approx_from = .9;
tolerance = 1e-1;
linearize_from = 0.92;

lb_approx = 1 - approx_to;
ub_approx = 1 - linearize_from;

dist_scale = 1e-8;

% piecewise linear numerical approximation of cauchy quantile in [approx_from, approx_to]
c_x_0 = 0;  % known x_p

c_f = @(x) 1 / (pi * (1 + x^2)) ;
c_d_f = @(x) -2*x / (pi * (1 + x^2)^2) ;
c_dd_f = @(x) -(2 * (-3*x^2 +1)) / (pi * (1 + x^2)^3);
c_ddd_f = @(x) 24*x * (-x^2 + 1) / (pi * (1 + x^2)^4);

[q_c_m, q_c_c] = quantile_approx(1, by, approx_to, 0.5, c_x_0, ...
    c_f, c_d_f, c_dd_f, c_ddd_f, tolerance, linearize_from);
q_c_m = q_c_m .* dist_scale;
q_c_c = q_c_c .* dist_scale;


% piecewise linear numerical approximation of squared norm of cauchy quantile in [alpha, norm_c_upper_bound]
norm_c_x_0 = 160.447638798;  % known x_p

norm_c_f = @(x)         2/(pi*sqrt(1+x)*(2+x));
norm_c_d_f = @(x)      -2/(pi*sqrt(1+x)*(2+x)^2) - 1/(pi*sqrt(1+x)^3*(2+x));
norm_c_dd_f = @(x)      4/(pi*sqrt(1+x)*(2+x)^3) + 2/(pi*sqrt(1+x)^3*(2+x)^2) + 3/(2*pi*sqrt(1+x)^5*(2+x));
norm_c_ddd_f = @(x)   -12/(pi*sqrt(1+x)*(2+x)^4) - 6/(pi*sqrt(1+x)^3*(2+x)^3) - 9/(2*pi*sqrt(1+x)^5*(2+x)^2) - 15/(4*pi*sqrt(1+x)^7*(2+x));

[q_norm_c_m, q_norm_c_c] = quantile_approx(1, by, approx_to, approx_from, norm_c_x_0, ...
    norm_c_f, norm_c_d_f, norm_c_dd_f, norm_c_ddd_f, tolerance, linearize_from);
q_norm_c_m = q_norm_c_m .* dist_scale;
q_norm_c_c = q_norm_c_c .* dist_scale;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;

% input initilizations
U_a_p = U_default;
U_b_p = U_default;
U_c_p = U_default;

solve

%%%%%%%%%%%%%%%%%%%%%%%%%%
% print some useful information
%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s \n', cvx_status);
fprintf('Computation time (sec): %f \n', total_time);
fprintf('Total Cost: %f \n', total_cost(k+1));
fprintf('Slack Cost: %f \n', lambda_sum(k+1));
fprintf('Input Cost: %f \n', input_cost(k+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify we meet the constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
    verify_nq(1).initial_condition = x_0_a;
    verify_nq(2).initial_condition = x_0_b;
    verify_nq(3).initial_condition = x_0_c;

    verify_nq(1).input = U_a;
    verify_nq(2).input = U_b;
    verify_nq(3).input = U_c;

    verify_nq(1).target_set = target_set_a;
    verify_nq(2).target_set = target_set_b;
    verify_nq(3).target_set = target_set_c;

    verification_nq = verify(10e5, sys, time_horizon, "L2", r, verify_nq, 4);
    verification_nq
else
    return
end