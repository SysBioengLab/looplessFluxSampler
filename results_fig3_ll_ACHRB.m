% Generate results for ll-ACHRB sampler
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
clearvars,clc
addpath('fxns','models');

% Initialize COBRA toolbox
try
    initCobraToolbox();
catch
    fprintf('COBRA Toolbox is not in the path. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
    return;
end

% Set up sampling parameters
options.numSamples    = 2e5;
options.stepsPerPoint = 1e2;
options.numDiscarded  = 0.1*options.numSamples;
options.algorithm     = 'll_ACHRB';
options.loopless      = 1;
options.warmUpFlag    = 0;
options.parallelFlag  = 1;
options.diagnostics   = 1;

% Load E coli core
load('e_coli_core.mat');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_ecoli_core_ll_ACHRB.mat sample
clearvars -except options

% Load iIT341
load('iIT341.mat');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_iIT341_ll_ACHRB.mat sample
clearvars -except options

% Load iYO844
load('iYO844.mat');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_iYO844_ll_ACHRB.mat sample
clearvars -except options

% Load iMM904
load('iMM904.mat');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_iMM904_ll_ACHRB.mat sample
clearvars -except options