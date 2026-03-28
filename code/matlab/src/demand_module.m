
function [fleet, retrofits, demand, ef] = demand_module(DATA_IN)
%DEMAND_MODULE  Load baseline fleet/demand/emission factors (stub)
% CSV schemas are described in README.md

fleet = struct(); retrofits = struct(); ef = struct();
% Demand as time series (stub: read CSV if present else linear ramp)
dfile = fullfile(DATA_IN, 'demand_baseline.csv');
if exist(dfile,'file')
    T = readtable(dfile);
    demand.t = T.year - T.year(1);
    demand.y = T.rpk;
else
    demand.t = linspace(0,10,11)';
    demand.y = linspace(1.0,1.2,11)';
end
end
