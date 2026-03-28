
function params = params_default()
%PARAMS_DEFAULT  Baseline parameter set (stub)
params = struct();
params.alpha = 0.2;      % adoption rate
params.beta  = 0.1;      % coupling to intensity
params.gamma = 0.05;     % demand tracking
params.lead  = 0.05;     % retrofit lead/drag
params.decay = 0.02;     % intensity decay
end
