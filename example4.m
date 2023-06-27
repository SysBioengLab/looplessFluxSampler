% Example 4. Comparison between ADSB against ll-ACHRB
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
%% Part 1: Initialization

clearvars,clc
addpath('fxns','test models');

Initialize COBRA toolbox
try
    initCobraToolbox();
catch
    fprintf('COBRA Toolbox is not in the path or install properly. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
    return;
end

%% Part 2: Sampling
% Set up seed for reproducible results
rng('default')

% Set up sampling parameters
options.numSamples    = 2e3; % Temporary: little number for tests
options.stepsPerPoint = 1e1;

% E coli
load('e_coli_core.mat');
options.algorithm = 'ADSB';
m1_ADSB = looplessFluxSampler(model,options);
save("samples\m1_ADSB","m1_ADSB")
clear m1_ADSB
options.algorithm = 'll_ACHRB';
m1_llACHRB = looplessFluxSampler(model,options);
save("samples\m1_llACHRB","m1_llACHRB")
clear m1_llACHRB model

% iIT341
load('iIT341.mat')
options.algorithm = 'ADSB';
m2_ADSB = looplessFluxSampler(model,options);
save("samples\m2_ADSB","m2_ADSB")
clear m2_ADSB
options.algorithm = 'll_ACHRB';
m2_llACHRB = looplessFluxSampler(model,options);
save("samples\m2_llACHRB","m2_llACHRB")
clear m2_llACHRB model

% iYO844
load('iYO844.mat')
options.algorithm = 'ADSB';
m3_ADSB = looplessFluxSampler(model,options);
save("samples\m3_ADSB","m3_ADSB")
clear m3_ADSB
options.algorithm = 'll_ACHRB';
m3_llACHRB = looplessFluxSampler(model,options);
save("samples\m3_llACHRB","m3_llACHRB")
clear m3_llACHRB model

% iMM904
load('iMM904.mat')
options.algorithm = 'ADSB';
m4_ADSB = looplessFluxSampler(model,options);
save("samples\m4_ADSB","m4_ADSB")
clear m4_ADSB
options.algorithm = 'll_ACHRB';
m4_llACHRB = looplessFluxSampler(model,options);
save("samples\m4_llACHRB","m4_llACHRB")
clear m4_llACHRB model
