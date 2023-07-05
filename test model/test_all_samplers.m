% Example 1. Benchamark of ADSB against ll-ACHRB
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
clearvars,clc
addpath('fxns','test model');

% % Initialize COBRA toolbox
% try
%     initCobraToolbox();
% catch
%     fprintf('COBRA Toolbox is not in the path or install properly. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
%     return;
% end

% Load model
load('e_coli_core.mat');
model = changeRxnBounds(model,'EX_o2_e',-.1,'u');

% Set up sampling parameters
options.numSamples    = 1e5;
options.stepsPerPoint = 2e1;
options.loopless      = 0;
options.warmUpFlag    = 0;
options.parallelFlag  = 1;
options.numCores      = 2;
options.diagnostics   = 1;

% Test1
options.algorithm = 'HRB';
sample = looplessFluxSampler(model,options);
save('sample_HRB.mat','sample');

% Test2
options.algorithm = 'll_ACHRB';
sample = looplessFluxSampler(model,options);
save('sample_ll_ACHRB.mat','sample');

% Test3
options.algorithm = 'ADSB';
sample = looplessFluxSampler(model,options);
save('sample_ADSB.mat','sample');
return

% Performance comparison
display(['Average time per effective sample ADSB: ',num2str(sample_ADSB.samplingTime/mean(sample_ADSB.Neff(isfinite(sample_ADSB.Neff))))]);
display(['Average time per effective sample ll-ACHRB: ',num2str(sample_ll_ACHRB.samplingTime/mean(sample_ll_ACHRB.Neff(isfinite(sample_ll_ACHRB.Neff))))]);
