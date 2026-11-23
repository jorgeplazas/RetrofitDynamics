% Mod1: coherent deployment-feedback model
% Retained Mod3/Mod4_c1 model. Illustrative parameters, not calibration.
% t: year; x: passenger-km/year; co2: kg CO2/year; z: deployment index.
a11 = 0.05;  % activity growth [1/year]
a12 = -0.02; % emissions-to-activity feedback [passenger-km/(kg year)]
a21 = 0.03;  % activity-to-emissions forcing [kg/(passenger-km year)]
a22 = -0.04; % background emissions decay [1/year]
beta = 0.8;  % deployment-dependent emissions decay [1/year]
Nx = 1e12;   % historical scaling constant; only delta/Nx enters the ODE
gamma = 0.1; % deployment exit rate [1/year]
delta = 0.02;% historical coefficient; chi = delta/Nx = 2e-14 [1/kg]
% delta is NOT an obsolescence rate in these retained equations.
y_init = [5e11; 1e9; 0.1];
tspan = linspace(0, 50, 500);
opts = odeset('RelTol', 1e-9, 'AbsTol', [1e-2 1e-3 1e-12], ...
    'MaxStep', 0.1, 'Events', @domain_boundary);
out_dir = fullfile(pwd, 'outputs_matlab');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
eta = 0.6;  % deployment entry rate [1/year], NOT fractional efficiency
sol = ode45(@(t,y) co2_system(t,y,a11,a12,a21,a22,beta,Nx,eta,gamma,delta), ...
    [0 50], y_init, opts);
check_solution(sol);
t = tspan(:);
Y = deval(sol, tspan).';
% Cumulative mass is post-processing of the annual emissions rate.
E = cumtrapz(t, Y(:,2));
fig = figure('Color', 'w');
subplot(3,1,1); plot(t,Y(:,1),'b','LineWidth',2);
ylabel('Activity (passenger-km/year)'); grid on;
subplot(3,1,2); plot(t,Y(:,2)/1e9,'r','LineWidth',2);
ylabel('CO2 rate (Mt/year)'); grid on;
subplot(3,1,3); plot(t,Y(:,3),'g','LineWidth',2);
ylabel('Deployment index'); xlabel('Time (years)'); grid on;
T = array2table([t Y E], 'VariableNames', ...
    {'time_year','x_pkm_per_year','q_kg_per_year','z','E_kg'});
writetable(T, fullfile(out_dir, 'Mod1_trajectory.csv'));
saveas(fig, fullfile(out_dir, 'Mod1_states.png'));

function dydt = co2_system(t, y, a11, a12, a21, a22, beta, Nx, eta, gamma, delta)
    % RHS preserved from Mod3 and Mod4_c1. eta is a rate [1/year].
    x = y(1);
    co2 = y(2);
    z = y(3);
    dxdt = a11 * x + a12 * co2;
    dco2dt = a21 * x + a22 * co2 - beta * co2 * z;
    dzdt = eta * (1 - z) - gamma * z + delta * co2 / Nx;
    dydt = [dxdt; dco2dt; dzdt];
end

function [value, isterminal, direction] = domain_boundary(~, y)
    % Stop on outward boundary crossings; never clip states or alter the RHS.
    value = [y(1)/1e12; y(2)/1e12; y(3); 1-y(3)];
    isterminal = ones(4,1);
    direction = -ones(4,1);
end

function check_solution(sol)
    if sol.x(end) < 50 - 1e-8 || (isfield(sol, 'xe') && ~isempty(sol.xe))
        error('Retrofit:DomainExit', ...
            'Integration stopped at a physical boundary. Revise the scenario; do not discard it silently.');
    end
    if any(~isfinite(sol.y(:))) || any(sol.y(1:2,:) < 0, 'all') || ...
            any(sol.y(3,:) < 0) || any(sol.y(3,:) > 1)
        error('Retrofit:InvalidState', 'Non-finite or out-of-domain state.');
    end
end
