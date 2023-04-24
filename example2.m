% Example 2. Use samples from a previous run of ADSB for constructing the initial current set 
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
clearvars,clc
addpath('fxns','test model');

% Initialize COBRA toolbox
% try
%     initCobraToolbox();
% catch
%     fprintf('COBRA Toolbox is not in the path or install properly. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
%     return;
% end

% Load model
load('e_coli_core.mat');

% Set up sampling parameters
options.numSamples    = 1e4;
options.stepsPerPoint = 1e2;

% ADSB (By default this algorithm is employed)
sample_ADSB1 = looplessFluxSampler(model,options);

% Re-use the previous points to construct the initial current set
model.points = sample_ADSB1.points;
sample_ADSB2 = looplessFluxSampler(model,options);

% Performance comparison
display(['Average time per effective sample ADSB (run 1): ',num2str(sample_ADSB1.samplingTime/mean(sample_ADSB1.Neff(isfinite(sample_ADSB1.Neff))))]);
display(['Average time per effective sample ADSB (run 2): ',num2str(sample_ADSB2.samplingTime/mean(sample_ADSB2.Neff(isfinite(sample_ADSB2.Neff))))]);
