% Example 1. Benchamark different samplers
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
clearvars,clc
addpath('fxns','.\scripts for paper\test model');

% Initialize COBRA toolbox
try
    initCobraToolbox();
catch
    fprintf('COBRA Toolbox is not in the path or install properly. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
    return;
end

% Load COBRA model
load('iIT341.mat');

% Set up sampling parameters
options.numSamples    = 5e4;    % Number of total samples
options.stepsPerPoint = 1e1;    % Thinning factor
options.loopless      = 0;      % 1 = yes; 0 = no
options.warmUpFlag    = 1;      % 1 = yes; 0 = no
options.parallelFlag  = 0;      % 1 = yes; 0 = no
options.numCores      = 1;      % Number of cores (only useful when running in parallel)
options.diagnostics   = 1;      % 1 = yes; 0 = no

% Test 1: ll-ACHRB
options.algorithm = 'll_ACHRB';
sample_ll_ACHRB = looplessFluxSampler(model,options);
save('sample_ll_ACHRB.mat','sample_ll_ACHRB');

% Test 2: ADSB
options.algorithm = 'ADSB';
sample_ADSB = looplessFluxSampler(model,options);
save('sample_ADSB.mat','sample_ADSB');

% Performance comparison
display(['Average time per effective sample ll-ACHRB: ',num2str(sample_ll_ACHRB.samplingTime/mean(sample_ll_ACHRB.Neff(isfinite(sample_ll_ACHRB.Neff))))]);
display(['Average time per effective sample ADSB: ',num2str(sample_ADSB.samplingTime/mean(sample_ADSB.Neff(isfinite(sample_ADSB.Neff))))]);