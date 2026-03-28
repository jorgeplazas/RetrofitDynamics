
%% Entry point: run all retrofit scenarios
% Add paths
addpath(genpath('src'));

% Configure I/O
DATA_IN = fullfile('..','..','data','input');
DATA_OUT = fullfile('..','..','data','output');
if ~exist(DATA_OUT,'dir'); mkdir(DATA_OUT); end
rng(1234,'twister');

% Load baseline parameters and data
params = params_default();
[fleet, retrofits, demand, ef] = demand_module(DATA_IN);

% Solve coupled ODE system (example time grid)
tspan = linspace(0,10,201); % years
x0 = [0;0;0]; % placeholder state
[t, x] = ode45(@(t,x) model_odes(t,x,params,demand,retrofits,ef), tspan, x0);

% Save & plot
save_results(DATA_OUT, t, x, params);
plot_trajectories(t, x, DATA_OUT);
