% Example 1. Benchamark of ll-ACHRB and EDHRB sampling algorithms
% -------------------- Copyright (C) 2018 Pedro A. Saa --------------------
clearvars,clc
addpath('looplessFluxSampler','test models');

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
options.numSamples    = 2e5;
options.stepsPerPoint = 1e1;

% EDHRB (By default this algorithm is employed)
sample_EDHRB = looplessFluxSampler(e_coli_core,options);

% ll-ACHRB
options.algorithm = 'll_ACHRB';
sample_ll_ACHRB = looplessFluxSampler(e_coli_core,options);

% Performance comparison
display(['Average time per effective sample EDHRB: ',num2str(sample_EDHRB.samplingTime/mean(sample_EDHRB.Neff(isfinite(sample_EDHRB.Neff))))]);
display(['Average time per effective sample ll-ACHRB: ',num2str(sample_ll_ACHRB.samplingTime/mean(sample_ll_ACHRB.Neff(isfinite(sample_ll_ACHRB.Neff))))]);