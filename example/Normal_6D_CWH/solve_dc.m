%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantile approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up
h = 5e-6; % step size
approx_ub = .9995;
lb_approx = 1 - approx_ub;
ub_approx = 1 - alpha;

% piecewise linear approximation of normal quantile in [0.5, norm_upper_bound]
norm_x_0 = 0;  % known x_p

norm_f = @(x) exp(-x^2/2) / sqrt(2*pi);
norm_d_f = @(x) -x * exp(-x^2/2) / sqrt(2*pi);
norm_dd_f = @(x) (x^2 -1) * exp(-x^2/2) / sqrt(2*pi);
norm_ddd_f = @(x) (3*x - x^3) * exp(-x^2/2) / sqrt(2*pi);

[qnorm_m, qnorm_c] = quantile_approx(1, h, approx_ub, 0.5, norm_x_0, ...
    norm_f, norm_d_f, norm_dd_f, norm_ddd_f, 1e-3, alpha);

% piecewise linear approximation of chi(3) quantile in [alpha, chi_upper_bound]
chi_x_0 = 2.1479;  % known x_p

chi_f = @(x) sqrt(2/pi) * x^2 * exp(-x^2/2);
chi_d_f = @(x) sqrt(2/pi) * (2 * x * exp(-x^2/2) - x^3 * exp(-x^2/2));
chi_dd_f = @(x) sqrt(2/pi) * (2 * exp(-x^2/2) - 5 * x^2 * exp(-x^2/2) + x^4 * exp(-x^2/2));
chi_ddd_f = @(x) sqrt(2/pi) * (-12 * x * exp(-x^2/2) + 9 * x^3 * exp(-x^2/2) - x^5 * exp(-x^2/2));

[qchi_m, qchi_c] = quantile_approx(1, h, approx_ub, alpha, chi_x_0, ...
    chi_f, chi_d_f, chi_dd_f, chi_ddd_f, 1e-3, alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

% input initilizations
U_a_p = zeros(input_dim * time_horizon, 1);
U_b_p = zeros(input_dim * time_horizon, 1);
U_c_p = zeros(input_dim * time_horizon, 1);

% max iterations
kmax = 100;
k = 1;

% convergence perameters
epsilon_dc = 1e-8; % convergence in cost
epsilon_lambda = 1e-8; % convergence of sum of slack variables to zero

% cost of slack variable
tau_max = 100;
gamma = 1.5;
tau = min(tau_max * ones(kmax,1), gamma.^(0:(kmax-1))');

% storage initial cost for convergence check
input_cost = [1e10; zeros(kmax,1)];
lambda_sum = [1e10; zeros(kmax,1)];
total_cost = [1e20; zeros(kmax,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = tic;
while k <= kmax 

    % update collision avoid probabilities and gradient
    if k ~= 1
        [g_ab, del_g_ab] = update_g(mean_X_a, mean_X_b, Cu, time_horizon, 0, 6);
        [g_ac, del_g_ac] = update_g(mean_X_a, mean_X_c, Cu, time_horizon, 0, 6);
        [g_bc, del_g_bc] = update_g(mean_X_b, mean_X_c, Cu, time_horizon, 0, 6);
    else
        [g_ab, del_g_ab] = update_g(mean_X_a_no_input + Cu * U_a_p, mean_X_b_no_input + Cu * U_b_p, Cu, time_horizon, 0, 6);
        [g_ac, del_g_ac] = update_g(mean_X_a_no_input + Cu * U_a_p, mean_X_c_no_input + Cu * U_c_p, Cu, time_horizon, 0, 6);
        [g_bc, del_g_bc] = update_g(mean_X_b_no_input + Cu * U_b_p, mean_X_c_no_input + Cu * U_c_p, Cu, time_horizon, 0, 6);
    end

    cvx_begin quiet
        variable U_a(input_dim * time_horizon,1);
        variable U_b(input_dim * time_horizon,1);
        variable U_c(input_dim * time_horizon,1);

        variable mean_X_a(state_dim * time_horizon, 1);
        variable mean_X_b(state_dim * time_horizon, 1);
        variable mean_X_c(state_dim * time_horizon, 1);

        variable lambda_ab(time_horizon, 1);
        variable lambda_ac(time_horizon, 1);
        variable lambda_bc(time_horizon, 1);
        
        variable gamma_ab(time_horizon, 1);
        variable gamma_ac(time_horizon, 1);
        variable gamma_bc(time_horizon, 1);
        
        variable chiapprox_ab(time_horizon, 1);
        variable chiapprox_ac(time_horizon, 1);
        variable chiapprox_bc(time_horizon, 1);

        variable delta_a(n_lin_state_a, 1);
        variable delta_b(n_lin_state_b, 1);
        variable delta_c(n_lin_state_c, 1);

        variable normapprox_a(n_lin_state_a, 1);
        variable normapprox_b(n_lin_state_b, 1);
        variable normapprox_c(n_lin_state_c, 1);

        minimize (tau(k)*(sum(lambda_ab) + sum(lambda_ac) + sum(lambda_bc)) + U_a'*U_a + U_b'*U_b + U_c'*U_c)
        subject to
            %----------------------------
            % linear equations defining the state
            %----------------------------
            mean_X_a == mean_X_a_no_input + Cu * U_a;
            mean_X_b == mean_X_b_no_input + Cu * U_b;
            mean_X_c == mean_X_c_no_input + Cu * U_c; 

            %----------------------------
            % u \in \mathcal{U} 
            %----------------------------
            input_space_A * U_a <= input_space_b;
            input_space_A * U_b <= input_space_b; 
            input_space_A * U_c <= input_space_b;

            %----------------------------
            % colission avoidance constraint
            %----------------------------

            % difference of convex function representation of 
            % ||x_a - x_b||^2 >= (r + Raylinv(\alpha)*min_eig(Sigma_w(k))^1/2)^2 - slack
            % slack variables added for feasibility.
            lambda_ab >= 0;
            lambda_ac >= 0;
            lambda_bc >= 0;
            
            for gamma_indx = 1:time_horizon
                chiapprox_ab(gamma_indx) >= qchi_m.* gamma_ab(gamma_indx) + qchi_c;
            end
            for gamma_indx = 1:time_horizon
                chiapprox_ac(gamma_indx) >= qchi_m.* gamma_ac(gamma_indx) + qchi_c;
            end
            for gamma_indx = 1:time_horizon
                chiapprox_bc(gamma_indx) >= qchi_m.* gamma_bc(gamma_indx) + qchi_c;
            end

            gamma_ab >= lb_approx;
            gamma_ac >= lb_approx;
            gamma_bc >= lb_approx;

            gamma_ab <= ub_approx;
            gamma_ac <= ub_approx;
            gamma_bc <= ub_approx;
            
            g_ab + del_g_ab * [U_a - U_a_p;U_b - U_b_p] >= ...
                (r + chiapprox_ab.*sigma_norm_lb).^2 - lambda_ab;
            g_ac + del_g_ac * [U_a - U_a_p;U_c - U_c_p] >= ...
                (r + chiapprox_ac.*sigma_norm_lb).^2 - lambda_ac;
            g_bc + del_g_bc * [U_b - U_b_p;U_c - U_c_p] >= ...
                (r + chiapprox_bc.*sigma_norm_lb).^2 - lambda_bc;

            %----------------------------
            % terminal state constraint
            %----------------------------

            % approximation of inverse normal in convex region
            for delta_indx = 1:n_lin_state_a
                normapprox_a(delta_indx) >= qnorm_m.* delta_a(delta_indx) + qnorm_c;
            end
            for delta_indx = 1:n_lin_state_b
                normapprox_b(delta_indx) >= qnorm_m.* delta_b(delta_indx) + qnorm_c;
            end
            for delta_indx = 1:n_lin_state_c
                normapprox_c(delta_indx) >= qnorm_m.* delta_c(delta_indx) + qnorm_c;
            end

            % \mu_v in target shrunk by \beta
            target_set_a.A * mean_X_a(end-5:end) + scaled_sigma_a_vec.* normapprox_a <= target_set_a.b;
            target_set_b.A * mean_X_b(end-5:end) + scaled_sigma_b_vec.* normapprox_b <= target_set_b.b;
            target_set_c.A * mean_X_c(end-5:end) + scaled_sigma_c_vec.* normapprox_c <= target_set_c.b;

            % \delta_i,v not infinity
            delta_a >= lb_approx;
            delta_b >= lb_approx;
            delta_c >= lb_approx;

            % \delta_i,v in convex region
            delta_a <= ub_approx;
            delta_b <= ub_approx;
            delta_c <= ub_approx;

            %----------------------------
            % overall safety
            %----------------------------
            sum(delta_a) + sum(delta_b) + sum(delta_c) <= 1 - alpha;
            sum(gamma_ab) + sum(gamma_ac) + sum(gamma_bc) <= 1 - alpha;
    cvx_end

    % update Costs
    input_cost(k+1) = U_a'*U_a + U_b'*U_b + U_c'*U_c;
    lambda_sum(k+1) = sum(lambda_ab) + sum(lambda_ac) + sum(lambda_bc);
    total_cost(k+1) = cvx_optval;

    % calculate convergence criteria
    conv_check = abs(input_cost(k+1) - input_cost(k)) + tau(k)*abs((lambda_sum(k+1) - lambda_sum(k)));

    % print statistics
    if ~ quiet
        fprintf('iteration: %d ', k);
        fprintf('\t %5.2f', cvx_optval);
        fprintf('\t %e', conv_check); 
        fprintf('\t %e', lambda_sum(k+1));
        fprintf('\t %4.2f', toc(start_time));
        fprintf('\t %s \n', cvx_status);
    end

    % check for solved status
    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        % check for convergence
        if (conv_check <= epsilon_dc) && (lambda_sum(k+1) <= epsilon_lambda)                 
           break
        end

        % if not converged update previous answer to current answer
        U_a_p = U_a;
        U_b_p = U_b;
        U_c_p = U_c;

    % if results are NaN break before error
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        break
    end

    % update itteration number
     k = k + 1;
end
total_time = toc(start_time);

% make k not less than or equal to kmax
k = min(k, kmax);


% print some useful information
fprintf('%s \n', cvx_status);
fprintf('Computation time (sec): %f \n', total_time);
fprintf('Total Cost: %f \n', total_cost(k+1));
fprintf('Slack Cost: %f \n', lambda_sum(k+1));
fprintf('Input Cost: %f \n', input_cost(k+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify we meet the constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
    verify_dc(1).initial_condition = x_0_a;
    verify_dc(2).initial_condition = x_0_b;
    verify_dc(3).initial_condition = x_0_c;

    verify_dc(1).input = U_a;
    verify_dc(2).input = U_b;
    verify_dc(3).input = U_c;

    verify_dc(1).target_set = target_set_a;
    verify_dc(2).target_set = target_set_b;
    verify_dc(3).target_set = target_set_c;

    verification_dc = verify(10e5, sys, time_horizon, "L2", r, verify_dc, 6);
    verification_dc
else
    return
end