
function plot_trajectories(t, x, DATA_OUT)
%PLOT_TRAJECTORIES  Basic figure export
figDir = fullfile(DATA_OUT,'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end
figure('Color','w');
plot(t, x, 'LineWidth', 1.5);
legend({'Adoption','Intensity','Demand'}, 'Location','best');
xlabel('Time (years)'); ylabel('State'); grid on
saveas(gcf, fullfile(figDir,'fig_trajectories.png'));
close gcf
end
