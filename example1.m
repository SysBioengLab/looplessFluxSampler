% Example1. Benchamark of ADSB against ll-ACHRB
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
clearvars,clc
addpath('fxns','test model');

% Initialize COBRA toolbox
try
    initCobraToolbox();
catch
    fprintf('COBRA Toolbox is not in the path or install properly. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
    return;
end

% Load model
load('e_coli_core.mat');

% Set up sampling parameters
options.numSamples    = 1e5;
options.stepsPerPoint = 2e1;

% ADSB (By default this algorithm is employed)
sample_ADSB = looplessFluxSampler(model,options);

% ll-ACHRB
options.algorithm = 'll_ACHRB';
sample_ll_ACHRB = looplessFluxSampler(model,options);

% Performance comparison
display(['Average time per effective sample ADSB: ',num2str(sample_ADSB.samplingTime/mean(sample_ADSB.Neff(isfinite(sample_ADSB.Neff))))]);
display(['Average time per effective sample ll-ACHRB: ',num2str(sample_ll_ACHRB.samplingTime/mean(sample_ll_ACHRB.Neff(isfinite(sample_ll_ACHRB.Neff))))]);
