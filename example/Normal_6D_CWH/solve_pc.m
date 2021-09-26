%%%%%%%%%%%%%%%%%%%%%%%%
% Blackmore 2011 for goal achievement and collision avoidance
% Primary Coder: Vignesh Sivaramakrishnan
% Modified by: Shawn Priore
%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
% big M arbitrary constant
large_constant = 5000;

% polytoupe defining ||x_i-x_j|| = r

if strcmpi(norm_choice, "L1")
    Avoid_A = [   1, 1, 1;...
                 -1,-1,-1;...
                 -1, 1, 1;...
                  1,-1,-1;
                  1,-1, 1;...
                 -1, 1,-1;...
                  1, 1,-1;...
                 -1,-1, 1];
    Avoid_b = r * sqrt(3) * ones(8,1);
elseif strcmpi(norm_choice, "Linf")
    Avoid_A = [  eye(3);
                -eye(3)];
    Avoid_b = r * ones(6,1);
end

% randomly generate the disturbance vector from the standard normal.
wvec = RandomVector('Gaussian', kron(ones(time_horizon,1), mu), kron(eye(time_horizon), sigma)); 
Wa = wvec.getRealizations(N);
Wb = wvec.getRealizations(N);
Wc = wvec.getRealizations(N);

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.

tic;
cvx_begin quiet
    variable U_a(input_dim * time_horizon,1);
    variable U_b(input_dim * time_horizon,1);
    variable U_c(input_dim * time_horizon,1);

    variable mean_X_a(state_dim * time_horizon, 1);
    variable mean_X_b(state_dim * time_horizon, 1);
    variable mean_X_c(state_dim * time_horizon, 1);    

    variable x_a(state_dim * time_horizon,N);
    variable x_b(state_dim * time_horizon,N);
    variable x_c(state_dim * time_horizon,N);

    variable t_a(N) binary;
    variable t_b(N) binary;
    variable t_c(N) binary;
    variable t_x(N) binary;

    variable c_ab(time_horizon * size(Avoid_A,1), N) binary;
    variable c_ac(time_horizon * size(Avoid_A,1), N) binary;
    variable c_bc(time_horizon * size(Avoid_A,1), N) binary;

    variable sc_ab(time_horizon , N) binary;
    variable sc_ac(time_horizon , N) binary;
    variable sc_bc(time_horizon , N) binary;

    variable ssc_ab(N) binary;
    variable ssc_ac(N) binary;
    variable ssc_bc(N) binary;
    variable ssc_x(N) binary;

    minimize (U_a'*U_a + U_b'*U_b + U_c'*U_c);

    subject to
        mean_X_a == mean_X_a_no_input + Cu * U_a;
        mean_X_b == mean_X_b_no_input + Cu * U_b;
        mean_X_c == mean_X_c_no_input + Cu * U_c;

        x_a(:,1:N) == Cw * Wa + repmat(mean_X_a,1,N);
        x_b(:,1:N) == Cw * Wb + repmat(mean_X_b,1,N);
        x_c(:,1:N) == Cw * Wc + repmat(mean_X_c,1,N);

        input_space_A * U_a <= input_space_b;
        input_space_A * U_b <= input_space_b; 
        input_space_A * U_c <= input_space_b;

        for i = 1:N

            target_set_a.A * x_a(end-5:end,i) - target_set_a.b <= large_constant*t_a(i);
            target_set_b.A * x_b(end-5:end,i) - target_set_b.b <= large_constant*t_b(i);
            target_set_c.A * x_c(end-5:end,i) - target_set_c.b <= large_constant*t_c(i);
            
            t_a(i) + t_b(i) + t_c(i) <= large_constant*t_x(i);

            for t = 1:time_horizon

                Avoid_A * (x_a(6*(t-1) + (1:3),i) - x_b(6*(t-1) + (1:3),i)) - Avoid_b + large_constant * c_ab(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i) >= 0; 
                Avoid_A * (x_a(6*(t-1) + (1:3),i) - x_c(6*(t-1) + (1:3),i)) - Avoid_b + large_constant * c_ac(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i) >= 0; 
                Avoid_A * (x_b(6*(t-1) + (1:3),i) - x_c(6*(t-1) + (1:3),i)) - Avoid_b + large_constant * c_bc(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i) >= 0; 

                sum(c_ab(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i)) - (size(Avoid_A,1) - 1) <= sc_ab(t, i);
                sum(c_ac(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i)) - (size(Avoid_A,1) - 1) <= sc_ac(t, i);
                sum(c_bc(size(Avoid_A,1)*(t-1) + (1:size(Avoid_A,1)) , i)) - (size(Avoid_A,1) - 1) <= sc_bc(t, i);
            end

            sum(sc_ab(:,i)) <= time_horizon * ssc_ab(i);
            sum(sc_ac(:,i)) <= time_horizon * ssc_ac(i);
            sum(sc_bc(:,i)) <= time_horizon * ssc_bc(i);
            
            ssc_ab(i) + ssc_ac(i) + ssc_bc(i) <= large_constant*ssc_x(i);
        end

        1/N * sum(t_x) <= 1-alpha;
        1/N * sum(ssc_x) <= 1-alpha;

cvx_end
total_time = toc;

fprintf('%s \n', cvx_status);
fprintf('Computation time (sec): %f \n', total_time);
fprintf('Input Cost: %f \n', cvx_optval);


%%
if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
    verify_pc(1).initial_condition = x_0_a;
    verify_pc(2).initial_condition = x_0_b;
    verify_pc(3).initial_condition = x_0_c;

    verify_pc(1).input = U_a;
    verify_pc(2).input = U_b;
    verify_pc(3).input = U_c;

    verify_pc(1).target_set = target_set_a;
    verify_pc(2).target_set = target_set_b;
    verify_pc(3).target_set = target_set_c;

    verification_pc = verify(10e5, sys, time_horizon, "L2", r, verify_pc, 6);
    verification_pc
else
    return
end