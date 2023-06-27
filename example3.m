% Example 3. Benchmark of scalability of ADSB against HRB
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
clearvars,clc
addpath('fxns','test models');

% Initialize COBRA toolbox
% try
%     initCobraToolbox();
% catch
%     fprintf('COBRA Toolbox is not in the path or install properly. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
%     return;
% end

%% Part 1: Model construction
% Load model
load('e_coli_core.mat');
% model = e_coli_core;
model = changeRxnBounds(model,'EX_o2_e',-0.1,'u'); % Aerobic regime

% Remove blocked reaction to maintain full-dimensional space
model = changeRxnBounds(model,'FRD7',0,'b');

% Model without loops
ZeroLoopModel = changeRxnBounds(model,'EX_o2_e',0,'b'); % Without o2 consumption there are no loops

% Model 2 loops
model2 = changeRxnBounds(model,'AKGDH',-1000,'l');
model2 = changeRxnBounds(model2,'AKGDH',1000,'u');
model2 = changeRxnBounds(model2,'GLUN',-1000,'l');
model2 = changeRxnBounds(model2,'GLUN',1000,'u');
model2 = changeRxnBounds(model2,'NADTRHD',-1000,'l');
model2 = changeRxnBounds(model2,'NADTRHD',1000,'u');
model2 = changeRxnBounds(model2,'PPC',-1000,'l');
model2 = changeRxnBounds(model2,'PPC',1000,'u');
% Model 10
model10 = changeRxnBounds(model,'CS',-1000,'l');
model10 = changeRxnBounds(model10,'CS',1000,'u');
model10 = changeRxnBounds(model10,'SUCCt2_2',-1000,'l');
model10 = changeRxnBounds(model10,'SUCCt2_2',1000,'u');
model10 = changeRxnBounds(model10,'SUCCt3',-1000,'l');
model10 = changeRxnBounds(model10,'SUCCt3',1000,'u');
model10 = changeRxnBounds(model10,'SUCDi',-1000,'l');
model10 = changeRxnBounds(model10,'SUCDi',1000,'u');
% Model 12
model12 = changeRxnBounds(model,'GND',-1000,'l');
model12 = changeRxnBounds(model12,'GND',1000,'u');
model12 = changeRxnBounds(model12,'NADH16',-1000,'l');
model12 = changeRxnBounds(model12,'NADH16',1000,'u');
model12 = changeRxnBounds(model12,'NADTRHD',-1000,'l');
model12 = changeRxnBounds(model12,'NADTRHD',1000,'u');
model12 = changeRxnBounds(model12,'PGL',-1000,'l');
model12 = changeRxnBounds(model12,'PGL',1000,'u');

%% Part 2: Sampling
% Set up seed for reproducible results
rng('default')

% Set up sampling parameters
options.numSamples    = 1e6;
options.stepsPerPoint = 1e1;

% Run with ADSB
R0 = looplessFluxSampler(ZeroLoopModel,options);
save('samples\R0.mat',"R0")
clear R0
R2 = looplessFluxSampler(model2,options);
save('samples\R2.mat',"R2")
clear R2
R10 = looplessFluxSampler(model10,options);
save('samples\R10.mat',"R10")
clear R10
R12 = looplessFluxSampler(model12,options);
save('samples\R12.mat',"R12")
clear R12

% Run with HRB
options.algorithm = 'HRB';
options.diagnostics = 0; % Without MCMC Diagnosis
R0_HRB = looplessFluxSampler(ZeroLoopModel,options);
save('samples\R0_HRB.mat',"R0_HRB")
clear R0_HRB
R2_HRB = looplessFluxSampler(model2,options);
save('samples\R2_HRB.mat',"R2_HRB")
clear R2_HRB
R10_HRB = looplessFluxSampler(model10,options);
save('samples\R10_HRB.mat',"R10_HRB")
clear R10_HRB
R12_HRB = looplessFluxSampler(model12,options);
save('samples\R12_HRB.mat',"R12_HRB")
clear R12_HRB
