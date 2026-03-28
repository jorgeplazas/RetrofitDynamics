
function configure_paths()
%CONFIGURE_PATHS  Ensure src is on MATLAB path
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
end
