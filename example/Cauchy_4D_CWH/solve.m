%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = tic;
while k <= kmax 

    % update collision avoid probabilities and gradient
    if k ~= 1
        [g_ab, del_g_ab] = update_g(mean_X_a, mean_X_b, Cu, time_horizon, r, 4);
        [g_ac, del_g_ac] = update_g(mean_X_a, mean_X_c, Cu, time_horizon, r, 4);
        [g_bc, del_g_bc] = update_g(mean_X_b, mean_X_c, Cu, time_horizon, r, 4);
    else
        [g_ab, del_g_ab] = update_g(mean_X_a_no_input + Cu * U_a_p, mean_X_b_no_input + Cu * U_b_p, Cu, time_horizon, r, 4);
        [g_ac, del_g_ac] = update_g(mean_X_a_no_input + Cu * U_a_p, mean_X_c_no_input + Cu * U_c_p, Cu, time_horizon, r, 4);
        [g_bc, del_g_bc] = update_g(mean_X_b_no_input + Cu * U_b_p, mean_X_c_no_input + Cu * U_c_p, Cu, time_horizon, r, 4);
    end

    cvx_begin quiet
        variable U_a(2 * time_horizon,1);
        variable U_b(2 * time_horizon,1);
        variable U_c(2 * time_horizon,1);

        variable mean_X_a(4 * time_horizon, 1);
        variable mean_X_b(4 * time_horizon, 1);
        variable mean_X_c(4 * time_horizon, 1);

        variable lambda_ab(time_horizon, 1);
        variable lambda_ac(time_horizon, 1);
        variable lambda_bc(time_horizon, 1);
        
        variable gamma_ab(time_horizon, 1);
        variable gamma_ac(time_horizon, 1);
        variable gamma_bc(time_horizon, 1);
        
        variable norm_c_approx_ab(time_horizon, 1);
        variable norm_c_approx_ac(time_horizon, 1);
        variable norm_c_approx_bc(time_horizon, 1);

        variable delta_a(n_lin_state_a, 1);
        variable delta_b(n_lin_state_b, 1);
        variable delta_c(n_lin_state_c, 1);

        variable cauchy_approx_a(n_lin_state_a, 1);
        variable cauchy_approx_b(n_lin_state_b, 1);
        variable cauchy_approx_c(n_lin_state_c, 1);

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
            % ||x_a - x_b||^2 >= (r + \hat{chi}*min_eig(Sigma_w(k))^1/2)^2 - slack
            % slack variables added for feasibility.
            lambda_ab >= 0;
            lambda_ac >= 0;
            lambda_bc >= 0;
            
            for gamma_indx = 1:time_horizon
                norm_c_approx_ab(gamma_indx) >= q_norm_c_m.* gamma_ab(gamma_indx) + q_norm_c_c;
                norm_c_approx_ac(gamma_indx) >= q_norm_c_m.* gamma_ac(gamma_indx) + q_norm_c_c;
                norm_c_approx_bc(gamma_indx) >= q_norm_c_m.* gamma_bc(gamma_indx) + q_norm_c_c;
            end

            gamma_ab >= lb_approx;
            gamma_ac >= lb_approx;
            gamma_bc >= lb_approx;

            gamma_ab <= ub_approx;
            gamma_ac <= ub_approx;
            gamma_bc <= ub_approx;
            
            g_ab + del_g_ab * [U_a - U_a_p;U_b - U_b_p] - norm_c_approx_ab .* max_sigma ./ dist_scale + lambda_ab >= ...
                 0;
            g_ac + del_g_ac * [U_a - U_a_p;U_c - U_c_p] - norm_c_approx_ac .* max_sigma ./ dist_scale + lambda_ac >= ...
                 0;
            g_bc + del_g_bc * [U_b - U_b_p;U_c - U_c_p] - norm_c_approx_bc .* max_sigma./ dist_scale  + lambda_bc >= ...
                 0;

            %----------------------------
            % terminal state constraint
            %----------------------------

            % approximation of inverse normal in convex region
            for delta_indx = 1:n_lin_state_a
                cauchy_approx_a(delta_indx) >= q_c_m.* delta_a(delta_indx) + q_c_c;
            end
            for delta_indx = 1:n_lin_state_b
                cauchy_approx_b(delta_indx) >= q_c_m.* delta_b(delta_indx) + q_c_c;
            end
            for delta_indx = 1:n_lin_state_c
                cauchy_approx_c(delta_indx) >= q_c_m.* delta_c(delta_indx) + q_c_c;
            end

            % \mu_v in target shrunk by \beta
            target_set_a.A * mean_X_a(end-3:end) + scaled_sigma_a_vec .* cauchy_approx_a ./ dist_scale - target_set_a.b <= ...
                0;
            target_set_b.A * mean_X_b(end-3:end) + scaled_sigma_b_vec .* cauchy_approx_b ./ dist_scale - target_set_b.b <= ...
                0;
            target_set_c.A * mean_X_c(end-3:end) + scaled_sigma_c_vec .* cauchy_approx_c ./ dist_scale - target_set_c.b <= ...
                0;

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
            sum(gamma_ab) + sum(gamma_ac) + sum(gamma_bc) <= 1 - beta;
    cvx_end

    % update Costs
    input_cost(k+1) = U_a'*U_a + U_b'*U_b + U_c'*U_c;
    lambda_sum(k+1) = sum(lambda_ab) + sum(lambda_ac) + sum(lambda_bc);
    total_cost(k+1) = cvx_optval;

    % calculate convergence criteria
    conv_check = abs(input_cost(k+1) - input_cost(k) + tau(k)*(lambda_sum(k+1) - lambda_sum(k)));

    % print statistics
    if ~ quiet
        fprintf('iteration: %d ', k);
        fprintf('\t %f', cvx_optval);
        fprintf('\t %e', conv_check); 
        fprintf('\t %e', lambda_sum(k+1));
        fprintf('\t %s', cvx_status);
        fprintf('\t %s \n', toc(start_time));
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
