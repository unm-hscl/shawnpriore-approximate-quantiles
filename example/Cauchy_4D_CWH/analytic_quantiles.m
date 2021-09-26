%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantile approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up
by = 5e-6; % step size
approx_to = .9995;
approx_from = .92;
tolerance = 1e-1;

lb_approx = 1 - approx_to;
ub_approx = 1 - approx_from;

dist_scale = 1e-8;

% piecewise linear analytical approximation of cauchy quantile in [approx_from, approx_to]
q_c = @(x) tan((x-0.5)*pi);

[q_c_m, q_c_c] = analytical_approx(1, by, approx_to, approx_from, q_c, tolerance);
q_c_m = q_c_m .* dist_scale;
q_c_c = q_c_c .* dist_scale;


% piecewise linear analytical approximation of squared norm of cauchy quantile in [alpha, norm_c_upper_bound]
q_norm_c = @(x) (tan((1+x)*pi/4)).^2 -1;

[q_norm_c_m, q_norm_c_c] = analytical_approx(1, by, approx_to, approx_from, q_norm_c, tolerance);
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
    verify_aq(1).initial_condition = x_0_a;
    verify_aq(2).initial_condition = x_0_b;
    verify_aq(3).initial_condition = x_0_c;

    verify_aq(1).input = U_a;
    verify_aq(2).input = U_b;
    verify_aq(3).input = U_c;

    verify_aq(1).target_set = target_set_a;
    verify_aq(2).target_set = target_set_b;
    verify_aq(3).target_set = target_set_c;

    verification_aq = verify(10e5, sys, time_horizon, "L2", r, verify_aq, 4);
    verification_aq
else
    return
end