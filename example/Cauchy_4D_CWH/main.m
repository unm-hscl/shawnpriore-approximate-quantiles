%%%%%%%%%%%%%%%%%%%%%%%%
% clear setup
%%%%%%%%%%%%%%%%%%%%%%%%

% clear system and vars
clc; 
clear all;
close all;
cvx_clear;

% add path
addpath '../'
addpath '../../'

%%%%%%%%%%%%%%%%%%%%%%%%
% problem setup
%%%%%%%%%%%%%%%%%%%%%%%%

% time params
time_horizon = 8;
sampling_period = 30; 

% initial states
% format: x, y, z,  x., y., z.
x_0_a = [100; 10.5; 0; 0] ; % satellite A
x_0_b = [106; 00;   0; 0] ; % satellite B
x_0_c = [ 94; 00;   0; 0] ; % satellite C

% target sets
% format: x, y, z, x., y., z.
target_set_a = Polyhedron('lb', [-15;   -10; -0.01; -0.01], ... 
                          'ub', [-10;    -5;  0.01;  0.01]);  
target_set_b = Polyhedron('lb', [- 2.5;  15; -0.01; -0.01], ...
                          'ub', [  2.5;  20;  0.01;  0.01]);    
target_set_c = Polyhedron('lb', [ 10;   -10; -0.01; -0.01], ... 
                          'ub', [ 15;    -5;  0.01;  0.01]);   

 
% collision avoid region radius
r = 15;

% safety threshold
alpha   = 0.9; 
beta    = 0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamics setup
%%%%%%%%%%%%%%%%%%%%%%%%%%

% LTI variables assumed IID
scale = [1e-4; 1e-4; 5e-8; 5e-8];

params = CwhSystemParameters('SamplingPeriod', sampling_period);

sys = getCwhLtiSystem(4, ...
                      Polyhedron('lb', -5*ones(2,1), ...
                                 'ub',  5*ones(2,1)), ...
                      RandomVector('UserDefined', ...
                                   @(N) scale.*trnd(1,[4 N])), ...
                      params);

% polytope representation of \mathcal{U}
[input_space_A, input_space_b] = getConcatInputSpace(sys, time_horizon);

% compute the concatenated matricies
[A, Cu, Cw] = getConcatMats(sys, time_horizon);

% propigate initial conditions through system
mean_X_a_no_input = A*x_0_a;
mean_X_b_no_input = A*x_0_b;
mean_X_c_no_input = A*x_0_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% constraint precomputation
%%%%%%%%%%%%%%%%%%%%%%%%%%

% collision avoidance 
max_sigma = zeros(time_horizon, 1);
rep_scale = kron(ones(time_horizon,1),scale); 
Cw_scaled = Cw * rep_scale;
for i = 1:time_horizon
    index = 4*(i-1) + (1:2); 
    max_sigma(i) =  2 * max(Cw_scaled(index));
end
max_sigma = max_sigma.^2;

% compute the number of polytopic halfspaces to worry about
n_lin_state_a = size(target_set_a.A,1);
n_lin_state_b = size(target_set_b.A,1);
n_lin_state_c = size(target_set_c.A,1);

CW_scaled_last = Cw_scaled(end-3:end);
scaled_sigma_a_vec = abs(target_set_a.A*CW_scaled_last);
scaled_sigma_b_vec = abs(target_set_b.A*CW_scaled_last);
scaled_sigma_c_vec = abs(target_set_c.A*CW_scaled_last);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization setup
%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization input
U_default = 0.5*ones(2 * time_horizon, 1);

% max iterations
kmax = 100;

% convergence perameters
epsilon_dc = 1e-8; % convergence in cost
epsilon_lambda = 1e-8; % convergence of sum of slack variables to zero

% cost of slack variable
tau_max = 100;
gamma = 1.2;
tau = min(tau_max * ones(kmax,1), gamma.^(0:(kmax-1))');

% storage initial cost for convergence check
input_cost = [1e10; zeros(kmax,1)];
lambda_sum = [1e10; zeros(kmax,1)];
total_cost = [1e20; zeros(kmax,1)];

% output
quiet = 0;

% Set defaults for cvx
cvx_solver gurobi
cvx_precision default

% find numerical results
numerical_quantiles

mean_X_a_nq = mean_X_a;
mean_X_b_nq = mean_X_b;
mean_X_c_nq = mean_X_c;
U_a_nq = U_a;
U_b_nq = U_b;
U_c_nq = U_c;

% find analytical results
analytic_quantiles

mean_X_a_aq = mean_X_a;
mean_X_b_aq = mean_X_b;
mean_X_c_aq = mean_X_c;
U_a_aq = U_a;
U_b_aq = U_b;
U_c_aq = U_c;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make some plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

red = [0.6350 0.0780 0.1840];
blue = [0.3010 0.7450 0.9330];
green = [0.4660 0.6740 0.1880];

fig = figure();
fig.Units    = 'inches';
fig.Position = [0.75,-1,13.5,11.5];
hold on
p1 = plot([x_0_a(1); mean_X_a_nq(1:4:end)], [x_0_a(2); mean_X_a_nq(2:4:end)], 'Color', red, 'Marker', 'h', 'LineWidth', 1);
p2 = plot([x_0_b(1); mean_X_b_nq(1:4:end)], [x_0_b(2); mean_X_b_nq(2:4:end)], 'Color', blue, 'Marker', 'p', 'LineWidth', 1);
p3 = plot([x_0_c(1); mean_X_c_nq(1:4:end)], [x_0_c(2); mean_X_c_nq(2:4:end)], 'Color', green, 'Marker', '^', 'LineWidth', 1);
plot([x_0_a(1); mean_X_a_aq(1:4:end)], [x_0_a(2); mean_X_a_aq(2:4:end)], '--', 'Color', red, 'Marker', 'h', 'LineWidth', 1)
plot([x_0_b(1); mean_X_b_aq(1:4:end)], [x_0_b(2); mean_X_b_aq(2:4:end)], '--', 'Color', blue, 'Marker', 'p', 'LineWidth', 1)
plot([x_0_c(1); mean_X_c_aq(1:4:end)], [x_0_c(2); mean_X_c_aq(2:4:end)], '--', 'Color', green, 'Marker', '^', 'LineWidth', 1)
p4 = plot(0,0,'k', 'LineWidth', 2);
p5 = plot(0,0,'k--', 'LineWidth', 2);
p6 = plot( polyshape(Polyhedron(target_set_a.A([1;2;5;6],1:2), target_set_a.b([1;2;5;6])).V),...
    'FaceColor', red, ...
    'FaceAlpha',0.1); 
plot( polyshape(Polyhedron(target_set_b.A([1;2;5;6],1:2), target_set_b.b([1;2;5;6])).V),...
    'FaceColor', blue, ...
    'FaceAlpha',0.1); 
plot( polyshape(Polyhedron(target_set_c.A([1;2;5;6],1:2), target_set_c.b([1;2;5;6])).V),...
    'FaceColor', green, ...
    'FaceAlpha',0.1); 
p7 = plot(x_0_a(1), x_0_a(2), 'Color', 'k', 'Marker', 'h', 'MarkerFaceColor', 'k', 'LineStyle','none');
plot(x_0_b(1), x_0_b(2), 'Color', 'k', 'Marker', 'p', 'MarkerFaceColor', 'k');
plot(x_0_c(1), x_0_c(2), 'Color', 'k', 'Marker', '^', 'MarkerFaceColor', 'k');

xlabel('x (in meters)')
ylabel('y (in meters)')
drawnow()
axis([-20 120 -15 30])
hold off
set(gca, 'OuterPosition', [0.025,0.025,0.975,0.925]);

l = legend([p1,p2,p3,p6,p4,p5,p7], {'Vehicle A', 'Vehicle B', 'Vehicle C', 'Target Set',...
                                       'Numerical', 'Analytical', 'Initial Location' },...
    'Orientation','horizontal', ...
    'Location', 'northoutside', ...
    'NumColumns', 4, ...
    'interpreter', 'latex');

axes('position',[.65 .65 .2 .2])
box on % put box around new pair of axes
hold on
plot(mean_X_a_nq(9:4:17), mean_X_a_nq(10:4:18), 'Color', red, 'Marker', 'h', 'LineWidth', 2);
plot(mean_X_a_aq(9:4:17), mean_X_a_aq(10:4:18), '--', 'Color', red, 'Marker', 'h', 'LineWidth', 2);
plot(mean_X_b_nq(13:4:21), mean_X_b_nq(14:4:22), 'Color', blue, 'Marker', 'p', 'LineWidth', 2);
plot(mean_X_b_aq(13:4:21), mean_X_b_aq(14:4:22), '--', 'Color', blue, 'Marker', 'p', 'LineWidth', 2);
axis([41.705 41.71 4.8155 4.8165])
hold off

set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     22);

%% Create empty vector
dist_ab_nq = zeros(time_horizon,1);
dist_ac_nq = zeros(time_horizon,1);
dist_bc_nq = zeros(time_horizon,1);
dist_ab_aq = zeros(time_horizon,1);
dist_ac_aq = zeros(time_horizon,1);
dist_bc_aq = zeros(time_horizon,1);

% Find distances between vehicles using L_2 norm
for i = 0:time_horizon-1
    dist_ab_nq(i+1) = norm(mean_X_a_nq((1:2)+(i*4))- mean_X_b_nq((1:2)+(i*4)));
    dist_ac_nq(i+1) = norm(mean_X_a_nq((1:2)+(i*4))- mean_X_c_nq((1:2)+(i*4)));
    dist_bc_nq(i+1) = norm(mean_X_b_nq((1:2)+(i*4))- mean_X_c_nq((1:2)+(i*4)));
    dist_ab_aq(i+1) = norm(mean_X_a_aq((1:2)+(i*4))- mean_X_b_aq((1:2)+(i*4)));
    dist_ac_aq(i+1) = norm(mean_X_a_aq((1:2)+(i*4))- mean_X_c_aq((1:2)+(i*4)));
    dist_bc_aq(i+1) = norm(mean_X_b_aq((1:2)+(i*4))- mean_X_c_aq((1:2)+(i*4)));
end
fig = figure();
fig.Units    = 'inches';
fig.Position = [0.75,-1,13.5,11.5];
hold on
p1 = plot((1:time_horizon), dist_ab_nq, 'Color', [0 0.4470 0.7410], 'Marker', 's');
p2 = plot((1:time_horizon), dist_ac_nq, 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
p3 = plot((1:time_horizon), dist_bc_nq, 'Color', [0.4940 0.1840 0.5560], 'Marker', 'd');
plot((1:time_horizon), dist_ab_aq, '--', 'Color', [0 0.4470 0.7410], 'Marker', 's');
plot((1:time_horizon), dist_ac_aq, '--', 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
plot((1:time_horizon), dist_bc_aq, '--', 'Color', [0.4940 0.1840 0.5560], 'Marker', 'd');
p6 = plot((1:time_horizon), r*ones(1,time_horizon), 'k-');
p4 = plot(1,10,'k');
p5 = plot(1,10,'k--');

drawnow()
hold off

xlabel('Time Step, k');
ylabel('Distance (in meters)');

l = legend([p1,p2,p3,p4,p5,p6], {'$||A-B||$','$||A-C||$','$||B-C||$', 'Numerical', 'Analytical', strcat('R=', num2str(r))},...
    'Orientation','horizontal', ...
    'Location', 'northoutside', ...
    'NumColumns', 3, ...
    'interpreter', 'latex');

axes('position',[.15 .65 .2 .2])
box on % put box around new pair of axes
hold on
plot((1:3), dist_ac_nq(1:3), 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
plot((1:3), dist_ac_aq(1:3), '--', 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
hold off
axis([1.9999 2.0001 18.9445 18.945])

set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     22);
