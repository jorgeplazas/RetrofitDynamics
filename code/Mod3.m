% Mod3: coherent deployment-feedback model
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
% eta is deployment entry rate [1/year]; none of these are efficiency fractions.
etas = [0.3 0.6 0.9];
colors = ['b' 'g' 'r'];
legend_entries = cell(1, numel(etas));
line_handles = gobjects(1, numel(etas));
slopes = zeros(numel(etas),1);
final_q = zeros(numel(etas),1);
final_E = zeros(numel(etas),1);
fig = figure('Color','w'); hold on;
for i = 1:numel(etas)
    eta = etas(i);
    sol = ode45(@(t,y) co2_system(t,y,a11,a12,a21,a22,beta,Nx,eta,gamma,delta), ...
        [0 50], y_init, opts);
    check_solution(sol);
    t = tspan(:);
    Y = deval(sol, tspan).';
    co2 = Y(:,2);
    % Evaluate at EXACTLY t=20 using the solver interpolant (not grid index 200).
    Y20 = deval(sol, 20);
    slopes(i) = (Y20(2) - y_init(2))/20; % kg/year^2
    E = cumtrapz(t,co2); % kg; trapezoidal approximation on the output grid
    final_q(i) = co2(end); final_E(i) = E(end);
    line_handles(i) = plot(t, co2/1e9, 'Color',colors(i),'LineWidth',2);
    legend_entries{i} = sprintf('eta = %.1f /yr; 0-20 yr slope %.3f Mt/yr^2', ...
        eta, slopes(i)/1e9);
    T = array2table([t Y E], 'VariableNames', ...
        {'time_year','x_pkm_per_year','q_kg_per_year','z','E_kg'});
    writetable(T,fullfile(out_dir,sprintf('Mod3_eta_%.1f.csv',eta)));
end
% No prespecified stabilization windows or unsupported growth-resumption labels.
title('Annual CO2 emissions by deployment entry rate');
xlabel('Time (years)'); ylabel('CO2 rate (Mt/year)');
legend(line_handles,legend_entries,'Location','northwest','Interpreter','none');
grid on; hold off;
summary = table(etas(:),slopes,final_q,final_E, 'VariableNames', ...
    {'eta_per_year','slope_0_20_kg_per_year2','q50_kg_per_year','E50_kg'});
writetable(summary,fullfile(out_dir,'Mod3_summary.csv'));
saveas(fig,fullfile(out_dir,'Mod3_deployment_scenarios.png'));

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
