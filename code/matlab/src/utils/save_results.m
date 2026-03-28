
function save_results(DATA_OUT, t, x, params)
%SAVE_RESULTS  Save trajectories and parameters
outDir = fullfile(DATA_OUT,'scenarios');
if ~exist(outDir,'dir'); mkdir(outDir); end
S = struct('t', t, 'x', x, 'params', params);
save(fullfile(outDir,'baseline.mat'), '-struct', 'S');
