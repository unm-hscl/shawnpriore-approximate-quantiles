%%%%%%%%%%%%%%%%%%%%%%%%
% clear setup
%%%%%%%%%%%%%%%%%%%%%%%%

% clear system and varsC
clc; 
clear all;
close all;
cvx_clear;

% add path
addpath '../'
addpath '../../../'

%%%%%%%%%%%%%%%%%%%%%%%%
% problem setup
%%%%%%%%%%%%%%%%%%%%%%%%

% define the system
input_dim = 3; % size of input vector U_i
state_dim = 6; % size of state vector (position and velocity) X_i

% bounds on input space
umax = 5; % max absolute value for an input

% time steps
time_horizon = 8; % number of steps from initial condition to completion
sampling_period = 30; 

% LTI variables
mu = zeros(state_dim,1); % mean distrubance
sigma = diag([1e-4, 1e-4, 1e-4, 5e-8, 5e-8, 5e-8]); % covariance of disturbance

% initial states
% format: x, y, z,  x., y., z.
x_0_a = [100; 10.5; 1;  0; 0; 0] ; % satellite A
x_0_b = [106; 00;   0;  0; 0; 0] ; % satellite B
x_0_c = [ 94; 00;  -1;  0; 0; 0] ; % satellite C

% target sets
% format: x, y, z, x., y., z.
target_set_a = Polyhedron('lb', [-15;   -10; 0; -0.01; -0.01; -0.01], ... 
                          'ub', [-10;    -5; 5;  0.01;  0.01;  0.01]);  
target_set_b = Polyhedron('lb', [- 2.5;  15; 0; -0.01; -0.01; -0.01], ...
                          'ub', [  2.5;  20; 5;  0.01;  0.01;  0.01]);    
target_set_c = Polyhedron('lb', [ 10;   -10; 0; -0.01; -0.01; -0.01], ... 
                          'ub', [ 15;    -5; 5;  0.01;  0.01;  0.01]);   

% collision avoid region radius
r = 15;

alpha = .9; % safety threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate LTI systems
%%%%%%%%%%%%%%%%%%%%%%%%%%

% LTI variables assumed IID
params = CwhSystemParameters('SamplingPeriod', sampling_period);

sys = getCwhLtiSystem(state_dim, ...
                      Polyhedron('lb', -umax*ones(input_dim,1), ...
                                 'ub',  umax*ones(input_dim,1)), ...
                      RandomVector('Gaussian', ...
                                   mu, ...
                                   sigma), ...
                      params);

sysnoi = LtvSystem('StateMatrix',sys.state_mat, ...
                   'DisturbanceMatrix', sys.dist_mat, ...
                   'Disturbance',sys.dist);

               
%%%%%%%%%%%%%%%%%%%%%%%%%%
% get needed matricies
%%%%%%%%%%%%%%%%%%%%%%%%%%

% polytope representation of \mathcal{U}
[input_space_A, input_space_b] = getConcatInputSpace(sys, time_horizon);

% compute the input concatenated transformations
[A, Cu, Cw] = getConcatMats(sys, time_horizon);

% compute mean trajectories without input
X_a_no_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_a, time_horizon);
X_b_no_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_b, time_horizon);
X_c_no_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_c, time_horizon);

mean_X_a_no_input = X_a_no_input_rv.mean();
mean_X_b_no_input = X_b_no_input_rv.mean();
mean_X_c_no_input = X_c_no_input_rv.mean();

mean_X_a_no_input = mean_X_a_no_input(sysnoi.state_dim+1:end);
mean_X_b_no_input = mean_X_b_no_input(sysnoi.state_dim+1:end);
mean_X_c_no_input = mean_X_c_no_input(sysnoi.state_dim+1:end);

% compute covariance of x 
cov_X_no_input = X_a_no_input_rv.cov();
cov_X_no_input = cov_X_no_input(sysnoi.state_dim+1:end, sysnoi.state_dim+1:end);

% compute the number of polytopic halfspaces to worry about
n_lin_state_a = size(target_set_a.A,1);
n_lin_state_b = size(target_set_b.A,1);
n_lin_state_c = size(target_set_c.A,1);

% compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i} = ||\sqrt\Sigma_X*h_i||
[sqrt_cov_X_no_input, p] = chol(cov_X_no_input(end-5:end,end-5:end));

if p > 0
    sqrt_cov_X_no_input = sqrt(cov_X_no_input);
end

scaled_sigma_a_vec = norms(target_set_a.A * sqrt_cov_X_no_input',2,2);
scaled_sigma_b_vec = norms(target_set_b.A * sqrt_cov_X_no_input',2,2);
scaled_sigma_c_vec = norms(target_set_c.A * sqrt_cov_X_no_input',2,2);

% collision avoidance
sigma_norm_lb = zeros(time_horizon, 1);
for i = 1:time_horizon
    index = 6*(i-1) + (1:3); 
    e = eig(2 * cov_X_no_input(index, index));
    sigma_norm_lb(i) =  sqrt(max(e));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set defaults for cvx
%%%%%%%%%%%%%%%%%%%%%%%%%%

cvx_solver gurobi
cvx_precision default

%%%%%%%%%%%%%%%%%%%%%%%%%%
% run solvers
%%%%%%%%%%%%%%%%%%%%%%%%%%

% find numerical results
quiet = 1;
solve_dc

mean_X_a_dc = mean_X_a;
mean_X_b_dc = mean_X_b;
mean_X_c_dc = mean_X_c;
U_a_dc = U_a;
U_b_dc = U_b;
U_c_dc = U_c;
keyboard
% find particle control results
N = 10;
norm_choice = 'L1';
solve_pc

mean_X_a_pc_1 = mean_X_a;
mean_X_b_pc_1 = mean_X_b;
mean_X_c_pc_1 = mean_X_c;
U_a_pc_1 = U_a;
U_b_pc_1 = U_b;
U_c_pc_1 = U_c;

norm_choice = 'Linf';
solve_pc

mean_X_a_pc_inf = mean_X_a;
mean_X_b_pc_inf = mean_X_b;
mean_X_c_pc_inf = mean_X_c;
U_a_pc_inf = U_a;
U_b_pc_inf = U_b;
U_c_pc_inf = U_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% make some plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

red = [0.6350 0.0780 0.1840];
blue = [0.3010 0.7450 0.9330];
green = [0.4660 0.6740 0.1880];

%%
fig = figure();
fig.Units    = 'inches';
fig.Position = [0.75,-1,13.5,11.5];
hold on
p1 = plot([x_0_a(1); mean_X_a_dc(1:6:end)], [x_0_a(2); mean_X_a_dc(2:6:end)], 'Color', red, 'Marker', 'h', 'LineWidth', 2);
p2 = plot([x_0_b(1); mean_X_b_dc(1:6:end)], [x_0_b(2); mean_X_b_dc(2:6:end)], 'Color', blue, 'Marker', 'p', 'LineWidth', 2);
p3 = plot([x_0_c(1); mean_X_c_dc(1:6:end)], [x_0_c(2); mean_X_c_dc(2:6:end)], 'Color', green, 'Marker', '^', 'LineWidth', 2);
plot([x_0_a(1); mean_X_a_pc_1(1:6:end)], [x_0_a(2); mean_X_a_pc_1(2:6:end)], '--', 'Color', red, 'Marker', 'h', 'LineWidth', 2);
plot([x_0_b(1); mean_X_b_pc_1(1:6:end)], [x_0_b(2); mean_X_b_pc_1(2:6:end)], '--', 'Color', blue, 'Marker', 'p', 'LineWidth', 2);
plot([x_0_c(1); mean_X_c_pc_1(1:6:end)], [x_0_c(2); mean_X_c_pc_1(2:6:end)], '--', 'Color', green, 'Marker', '^', 'LineWidth', 2);
plot([x_0_a(1); mean_X_a_pc_inf(1:6:end)], [x_0_a(2); mean_X_a_pc_inf(2:6:end)], '-.', 'Color', red, 'Marker', 'h', 'LineWidth', 2);
plot([x_0_b(1); mean_X_b_pc_inf(1:6:end)], [x_0_b(2); mean_X_b_pc_inf(2:6:end)], '-.', 'Color', blue, 'Marker', 'p', 'LineWidth', 2);
plot([x_0_c(1); mean_X_c_pc_inf(1:6:end)], [x_0_c(2); mean_X_c_pc_inf(2:6:end)], '-.', 'Color', green, 'Marker', '^', 'LineWidth', 2);
p4 = plot(0,0,'k', 'LineWidth', 2);
p5 = plot(0,0,'k--', 'LineWidth', 2);
p6 = plot(0,0,'k-.', 'LineWidth', 2);
p7 = plot( polyshape(Polyhedron(target_set_a.A([1;2;7;8],1:2), target_set_a.b([1;2;7;8])).V),...
    'FaceColor', red, ...
    'FaceAlpha',0.1); 
plot( polyshape(Polyhedron(target_set_b.A([1;2;7;8],1:2), target_set_b.b([1;2;7;8])).V),...
    'FaceColor', blue, ...
    'FaceAlpha',0.1); 
plot( polyshape(Polyhedron(target_set_c.A([1;2;7;8],1:2), target_set_c.b([1;2;7;8])).V),...
    'FaceColor', green, ...
    'FaceAlpha',0.1); 
p8 = plot(x_0_a(1), x_0_a(2), 'Color', 'k', 'Marker', 'h', 'MarkerFaceColor', 'k', 'LineStyle','none');
plot(x_0_b(1), x_0_b(2), 'Color', 'k', 'Marker', 'p', 'MarkerFaceColor', 'k');
plot(x_0_c(1), x_0_c(2), 'Color', 'k', 'Marker', '^', 'MarkerFaceColor', 'k');

xlabel('x (in meters)')
ylabel('y (in meters)')
drawnow()
axis([-20 120 -15 25])
hold off
set(gca, 'OuterPosition', [0.025,0.025,0.975,0.925]);

l = legend([p1,p2,p3,p7,p4,p5,p6,p8], {'Vehicle A', 'Vehicle B', 'Vehicle C', 'Target Set',...
                                       'Proposed Method', 'Particle Control, $\ell_1$', 'Particle Control, $\ell_\infty$', 'Initial Location' },...
    'Orientation','horizontal', ...
    'Location', 'northoutside', ...
    'NumColumns', 4, ...
    'interpreter', 'latex');

set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     22);


%%
% Create empty vector
dist_ab_dc = zeros(time_horizon,1);
dist_ac_dc = zeros(time_horizon,1);
dist_bc_dc = zeros(time_horizon,1);
dist_ab_pc_1 = zeros(time_horizon,1);
dist_ac_pc_1 = zeros(time_horizon,1);
dist_bc_pc_1 = zeros(time_horizon,1);
dist_ab_pc_inf = zeros(time_horizon,1);
dist_ac_pc_inf = zeros(time_horizon,1);
dist_bc_pc_inf = zeros(time_horizon,1);

% Find distances between vehicles using L_2ty norm
for i = 0:time_horizon-1
    dist_ab_dc(i+1) = norm(mean_X_a_dc((1:3)+(i*6))- mean_X_b_dc((1:3)+(i*6)));
    dist_ac_dc(i+1) = norm(mean_X_a_dc((1:3)+(i*6))- mean_X_c_dc((1:3)+(i*6)));
    dist_bc_dc(i+1) = norm(mean_X_b_dc((1:3)+(i*6))- mean_X_c_dc((1:3)+(i*6)));
    dist_ab_pc_1(i+1) = norm(mean_X_a_pc_1((1:3)+(i*6))- mean_X_b_pc_1((1:3)+(i*6)));
    dist_ac_pc_1(i+1) = norm(mean_X_a_pc_1((1:3)+(i*6))- mean_X_c_pc_1((1:3)+(i*6)));
    dist_bc_pc_1(i+1) = norm(mean_X_b_pc_1((1:3)+(i*6))- mean_X_c_pc_1((1:3)+(i*6)));
    dist_ab_pc_inf(i+1) = norm(mean_X_a_pc_inf((1:3)+(i*6))- mean_X_b_pc_inf((1:3)+(i*6)));
    dist_ac_pc_inf(i+1) = norm(mean_X_a_pc_inf((1:3)+(i*6))- mean_X_c_pc_inf((1:3)+(i*6)));
    dist_bc_pc_inf(i+1) = norm(mean_X_b_pc_inf((1:3)+(i*6))- mean_X_c_pc_inf((1:3)+(i*6)));
end
fig = figure();
fig.Units    = 'inches';
fig.Position = [0.75,-1,13.5,11.5];
hold on

p1 = plot((1:time_horizon), dist_ab_dc, 'Color', [0 0.4470 0.7410], 'Marker', 's');
p2 = plot((1:time_horizon), dist_ac_dc, 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
p3 = plot((1:time_horizon), dist_bc_dc, 'Color', [0.4940 0.1840 0.5560], 'Marker', 'd');
plot((1:time_horizon), dist_ab_pc_1, '--', 'Color', [0 0.4470 0.7410], 'Marker', 's');
plot((1:time_horizon), dist_ac_pc_1, '--', 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
plot((1:time_horizon), dist_bc_pc_1, '--', 'Color', [0.4940 0.1840 0.5560], 'Marker', 'd');
plot((1:time_horizon), dist_ab_pc_inf, '-.', 'Color', [0 0.4470 0.7410], 'Marker', 's');
plot((1:time_horizon), dist_ac_pc_inf, '-.', 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
plot((1:time_horizon), dist_bc_pc_inf, '-.', 'Color', [0.4940 0.1840 0.5560], 'Marker', 'd');
p7 = plot((1:time_horizon), r*ones(1,time_horizon), 'k-');

p4 = plot(1,10,'k');
p5 = plot(1,10,'k--');
p6 = plot(1,10,'k-.');

drawnow()
hold off

xlabel('Time Step, k');
ylabel('Distance (in meters)');

l = legend([p1,p2,p3,p7,p4,p5,p6], {'$||A-B||$','$||A-C||$','$||B-C||$',strcat('R=', num2str(r)), 'Proposed Method', 'Particle Control, $\ell_1$', 'Particle Control, $\ell_\infty$'},...
    'Orientation','horizontal', ...
    'Location', 'northoutside', ...
    'NumColumns', 4, ...
    'interpreter', 'latex');
    
save('full_run.mat');