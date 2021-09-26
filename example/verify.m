function results = verify(samples, sys, time_horizon, norm_choice, r, problem, dim) 
    % define variables
    x_0_a = problem(1).initial_condition;
    x_0_b = problem(2).initial_condition;
    x_0_c = problem(3).initial_condition;

    U_a = problem(1).input;
    U_b = problem(2).input;
    U_c = problem(3).input;
    
    target_set_a = problem(1).target_set;
    target_set_b = problem(2).target_set;
    target_set_c = problem(3).target_set;

    % get mc samples
    [state_realization_a, ~, ~] = ...
        generateMonteCarloSims(samples, sys, x_0_a, ...
        time_horizon, U_a);

    [state_realization_b, ~, ~] = ...
        generateMonteCarloSims(samples, sys, x_0_b, ...
        time_horizon, U_b);

    [state_realization_c, ~, ~] = ...
        generateMonteCarloSims(samples, sys, x_0_c, ...
        time_horizon, U_c);

    % test goal achievement
    P_goal_a = zeros(samples, 1);
    P_goal_b = zeros(samples, 1);
    P_goal_c = zeros(samples, 1);

    for i = 1:samples
       P_goal_a(i) = target_set_a.contains( state_realization_a(end-dim+1:end, i) );
       P_goal_b(i) = target_set_b.contains( state_realization_b(end-dim+1:end, i) );
       P_goal_c(i) = target_set_c.contains( state_realization_c(end-dim+1:end, i) );
    end
    
    % test collision avoidance by time index and sample
    diff_ab = state_realization_a - state_realization_b;    
    diff_ac = state_realization_a - state_realization_c;
    diff_bc = state_realization_b - state_realization_c;
    
    norms_ab = zeros(time_horizon, samples);
    norms_ac = zeros(time_horizon, samples);
    norms_bc = zeros(time_horizon, samples);
    
    if strcmpi(norm_choice, "L2")
        for t = 1:time_horizon
            for i = 1:samples
                norms_ab(t, i) = norm( diff_ab(dim*t + (1:(dim/2)), i) );
                norms_ac(t, i) = norm( diff_ac(dim*t + (1:(dim/2)), i) );
                norms_bc(t, i) = norm( diff_bc(dim*t + (1:(dim/2)), i) );
            end
        end
    elseif strcmpi(norm_choice, "Linf")
        for t = 1:time_horizon
            for i = 1:samples
                norms_ab(t, i) = norm( diff_ab(dim*t + (1:(dim/2)), i) , Inf );
                norms_ac(t, i) = norm( diff_ac(dim*t + (1:(dim/2)), i) , Inf );
                norms_bc(t, i) = norm( diff_bc(dim*t + (1:(dim/2)), i) , Inf );
            end
        end
    end
    
    % test collision avoidance by sample
    norm_ab_sample = (sum( (norms_ab >= r) , 1) == time_horizon)';
    norm_ac_sample = (sum( (norms_ac >= r) , 1) == time_horizon)';
    norm_bc_sample = (sum( (norms_bc >= r) , 1) == time_horizon)';
    
    % test overall safety by sample
    overall_collision_avoidance = norm_ab_sample & norm_ac_sample & norm_bc_sample;
    overall_target_achievement = P_goal_a & P_goal_b & P_goal_c;
    overall = overall_collision_avoidance & overall_target_achievement;

    % compute and return statistics
    results.P_goal_a = sum(P_goal_a) / samples;
    results.P_goal_b = sum(P_goal_b) / samples;
    results.P_goal_c = sum(P_goal_c) / samples;
    
    results.P_collision_ab = sum(norm_ab_sample) / samples;
    results.P_collision_ac = sum(norm_ac_sample) / samples;
    results.P_collision_bc = sum(norm_bc_sample) / samples;
    
    results.overall_collision_avoidance     = sum(overall_collision_avoidance) / samples;
    results.overall_target_achievement      = sum(overall_target_achievement) / samples;
    results.overall                         = sum(overall) / samples;
end