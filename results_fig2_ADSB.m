% Generate results for ADSB sampler
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
options.algorithm     = 'ADSB';
options.loopless      = 1;
options.warmUpFlag    = 0;
options.parallelFlag  = 1;
options.diagnostics   = 1;

% Load model R0
load('e_coli_core.mat');
model = changeRxnBounds(model,'FRD7',0,'b');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_R0_ADSB.mat sample
clearvars -except options

% Run model R2
load('e_coli_core.mat');
model = changeRxnBounds(model,'FRD7',0,'b');
model = changeRxnBounds(model,{'AKGDH', 'GLUN', 'NADTRHD', 'PPC'},-1e3,'l');
model = changeRxnBounds(model,{'AKGDH', 'GLUN', 'NADTRHD', 'PPC'},1e3,'u');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_R2_ADSB.mat sample
clearvars -except options

% Run model R10
load('e_coli_core.mat');
model = changeRxnBounds(model,'FRD7',0,'b');
model = changeRxnBounds(model,{'CS', 'SUCCt2_2', 'SUCCt3', 'SUCDi'},-1e3,'l');
model = changeRxnBounds(model,{'CS', 'SUCCt2_2', 'SUCCt3', 'SUCDi'},1e3,'u');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_R10_ADSB.mat sample
clearvars -except options

% Run model R12
load('e_coli_core.mat');
model = changeRxnBounds(model,'FRD7',0,'b');
model = changeRxnBounds(model,{'GND', 'NADH16', 'NADTRHD', 'PGL'},-1e3,'l');
model = changeRxnBounds(model,{'GND', 'NADH16', 'NADTRHD', 'PGL'},1e3,'u');
sample = looplessFluxSampler(model,options);
save -v7.3 sample_R12_ADSB.mat sample
clearvars -except options