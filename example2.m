% Example 2. Estimation of `optimal` thining factor for EDHRB
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
load('iIT341.mat');

% Set up sampling parameters
options.numSamples    = 2e4;
options.stepsPerPoint = 1e0;

% EDHRB
sample = looplessFluxSampler(iIT341,options);
tau1   = mean(sample.tau(isfinite(sample.tau)));

% Get the maximum thining factor and use it to produce the final sample
sample.stepsPerPoint = max(sample.thin(isfinite(sample.tau)));
[sample.points,sample.samplingTime] = EDHRB(sample,1);
sample.points = reshape(sample.points,sample.numRxns,sample.numChains*sample.pointsPerChain);

% Split samples resulting chain in 10 segments for convergence analysis
sample.points = reshape(sample.points,sample.numRxns,fix(sample.numSamples/10),10);
sample.points = permute(sample.points,[2,1,3]);

% Calculate potential scale reduction statistics and recover original chain
[sample.R,sample.Rint,sample.Neff,sample.tau,sample.thin] = psrf(sample.points);
sample.points = permute(sample.points,[2,3,1]);
sample.points = reshape(sample.points,sample.numRxns,sample.numSamples);
disp(['Average autocorrelation time (no thining): ',num2str(tau1)]);
disp(['Average autocorrelation time (`optimal` thining): ',num2str(mean(sample.tau(isfinite(sample.tau))))]);