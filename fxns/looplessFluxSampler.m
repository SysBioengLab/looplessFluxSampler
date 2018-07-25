function sample = looplessFluxSampler(model,options)
% Loopless Flux Sampler
%
% Performs random sampling of the loopless space of metabolic models 
%
% USAGE:
%              sample = looplessFluxSampler(model, options)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x 1` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * rxns - `n x 1` rxn identifiers (cell array)
%              options (structure):  (the following fields are required - others can be supplied)
%                                    * numSamples - number of points (double)
%
% OPTIONAL INPUTS:
%              options (structure):  (the following fields are optional)
%                                    * numDiscarded - Burn-in (double)
%                                    * stepsPerPoint - Thining (double)
%                                    * algorithm - 'EDHRB' (default) or 'll_ACHRB'
%                                    * loopless - Loop removal option, true (default) or false
%                                    * stepFxn - Step function option,'linear' (default) or 'uniform'
%                                    * vTol - Numerical flux tolerance  (default 1e-8)
%                                    * parallelFlag - Parallel sampling option, true or false (default)
%
% OUTPUT:
%              sample (structure):   sampling structure containing #numSamples random points
%
% -------------------- Copyright (C) 2018 Pedro A. Saa --------------------

% Check inputs
if (nargin<2)
    fprintf('Not enough input arguments.');
    return;
else
    % Check if COBRA toolbox is in the path
    if ~exist('initCobraToolbox','file')
        fprintf('COBRA Toolbox is not in the path. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
        return;
    end
    
    % Number of samples
    if isfield(options,'numSamples'); sample.numSamples = options.numSamples;
    else fprintf('Field numSamples has not been defined.'); return; end;
    
    % Samples discarded or burn-in (default 10% of samples)
    if isfield(options,'numDiscarded'); sample.numDiscarded = options.numDiscarded;
    else sample.numDiscarded = fix(.1*sample.numSamples); end;
    
    % Steps per point or thinning (default 1e2)
    if isfield(options,'stepsPerPoint'); sample.stepsPerPoint = options.stepsPerPoint;
    else sample.stepsPerPoint = 1e2; end;
    
    % Sampling algorithm 'EDHRB' (default) or 'll_ACHRB'
    if isfield(options,'algorithm'); sample.algorithm = options.algorithm;
    else sample.algorithm = 'EDHRB'; end;
    
    % Sampling mode (default loopless true or 1)
    if isfield(options,'loopless'); sample.loopless = options.loopless;
    else sample.loopless = 1; end;
    
    % Step-calculation fxn 'linear' (default) or 'uniform'
    if isfield(options,'stepFxn'); sample.stepFxn = options.stepFxn;
    else sample.stepFxn = 'linear'; end;
    
    % Numerical flux tolerance  (default 1e-8)
    if isfield(options,'vTol'); sample.vTol = options.vTol;
    else sample.vTol = 1e-8; end;
    
    % Parallel sampling option (default false or 0)
    sample.parallelFlag  = 0;
    if isfield(options,'parallelFlag'); sample.parallelFlag = options.parallelFlag;
    else sample.parallelFlag = 0; end;
end

% Define solver parameters
changeCobraSolverParams('MILP','intTol',1e-9);
changeCobraSolverParams('MILP','relMipGapTol',1e-12);
changeCobraSolverParams('LP','optTol',1e-9);
changeCobraSolverParams('LP','feasTol',1e-9);

% Initialize temporal model structure
tempModel = model;
[tempModel,rxnList] = parseInternalRxns(tempModel,model.rxns);
fprintf('Model loaded: initial model contains %d rxns and %d mets.\n',numel(model.rxns),numel(model.mets));

%% I. Model reduction
% Time model pre-processing and sampling preparation
t0 = cputime;

% Reduce model by removing singletons. First, build LP structure
LPproblem = buildLinearProblem(tempModel);

% Run FVA to remove singletons
[tempModel.lb,tempModel.ub] = generalFVA(LPproblem,sample.vTol);
[tempModel,zeroRxns] = removeBlockedSets(tempModel);
rxnList(zeroRxns)    = [];
[tempModel,rxnList]  = parseInternalRxns(tempModel,rxnList);

% Run ll-FVA to remove rxns involved in infeasible loops active
if (sample.loopless)
    
    % Get a sparse null-space matrix of internal rxns
    Nint = fastSNP(tempModel.S(:,tempModel.intRxns),-1e2*(tempModel.lb(tempModel.intRxns)<0),...
        1e2*(tempModel.ub(tempModel.intRxns)>0),sample.vTol);
    if ~isempty(Nint)
        looplessProblem = buildLooplessProblem(tempModel,normMatrixEntries(Nint)');
        [tempModel.lb,tempModel.ub] = generalFVA(looplessProblem,sample.vTol,'loopless');
        [tempModel,zeroRxns] = removeBlockedSets(tempModel);
        rxnList(zeroRxns)    = [];
        [tempModel,rxnList]  = parseInternalRxns(tempModel,rxnList);
    else sample.loopless = 0; end;
    
    % Get an updated sparse null-space matrix of internal rxns (if relevant)
    if sample.loopless
        Nint = fastSNP(tempModel.S(:,tempModel.intRxns),-1e2*(tempModel.lb(tempModel.intRxns)<0),...
            1e2*(tempModel.ub(tempModel.intRxns)>0),sample.vTol);
        if ~isempty(Nint); Nint = normMatrixEntries(Nint);
        else sample.loopless = 0; end;
    end
    if sample.loopless; fprintf('The model has %d infeasible loop law(s) potentially active.\n',size(Nint,2));
    else fprintf('The model does not contain infeasible loops active.\n'); end;
end
fprintf('Model reduced: pre-processed model has %d rxns.\n',numel(rxnList));

%% II. Preparation phase
% Define fields related to the model in sample structure
sample.rxns     = rxnList;
sample.numRxns  = numel(tempModel.lb);
sample.S        = tempModel.S;
sample.lb       = tempModel.lb;
sample.ub       = tempModel.ub;
sample.intRxns  = tempModel.intRxns;
sample.exchRxns = tempModel.exchRxns;
sample.uTol     = 1e-10;                                                   % Direction tolerance
sample.bTol     = 1e-10;                                                   % Minimum allowed distance to closest constraint
if sample.loopless
    sample.Nint       = Nint;
    sample.nnzEntries = sum(Nint~=0);
    sample.loopMatrix = sparse(sign(Nint));
    
    % Define thermodynamic feasibility oracle
    sample.isFeasible = @(points) looplessCheck(points,sample);
end

% Define fxn to ensure samples are kept within the bounds
sample.keepWithinBounds = @(points) bringToBoundary(points,sample.lb,sample.ub);

% Generate initial warmup points
sample.warmUpPoints = warmupLooplessACHRB(sample);

% Calculate centroid
centroid = mean(sample.warmUpPoints,2);
centroid(abs(centroid)<sample.vTol) = 0;

% Check feasibility of the initial point. If centroid infeasible, find closest feasible point
if sample.loopless && ~sample.isFeasible(centroid)
    sample.centroid = findNearestFeasiblePoint(sample,centroid,Nint',sample.vTol);
else
    sample.centroid = centroid;
end

% Calculate centroid and seeds for EDHRB 
if strcmp(sample.algorithm,'EDHRB')    
    
    % Run ll-ACHRB and generate initial seedpoints (cPoints)
    numSamples = 2e4;
    cWeight    = size(sample.warmUpPoints,2);
    cCentroid  = [];
    cPoints    = [];
    while true
        [vPoints,~,sample.centroid,cWeight] = ll_ACHRB(sample,numSamples,0,1,1,cWeight);
        if isempty(vPoints); return; end;
        
        % Calculate centroid trace evolution
        cPoints   = [cPoints,vPoints];
        cCentroid = [cCentroid,bsxfun(@rdivide,cumsum(vPoints,2),size(cCentroid,2)+1:size(cCentroid,2)+size(vPoints,2))];
        Rfactor   = psrf(cCentroid');
        
        % Check the convergence of the reduction factor
        if (abs(median(Rfactor(isfinite(Rfactor)))-1) < .1); break;
        else numSamples = 2*numSamples; end;
    end
    
    % Calculate directions norms and select those spanning elongated
    % coordinates
    sample.udir     = bsxfun(@minus,cPoints,sample.centroid);
    sample.udir(:,~sample.isFeasible(sample.udir)) = [];                                                   % Drop infeasible directions
    numDir          = min([rank(sample.warmUpPoints)*sample.numRxns,size(sample.udir,2)]);                 % Define number of directions to employ
    [dLength,ixDir] = sort(sqrt(sum(sample.udir.^2)),'descend');
    sample.udir     = bsxfun(@rdivide,sample.udir(:,ixDir(1:numDir)),dLength(ixDir(1:numDir)));
    
    % Finally, define number of chains, points per chain and initial points
    maxNumChains          = 1e2;
    sample.numChains      = min([max([rank(sample.warmUpPoints),maxNumChains]),maxNumChains]);             % Run 1e2 independent parallel chains
    sample.pointsPerChain = ceil(sample.numSamples/sample.numChains);                                      % Define number of points per chain
    
    % Sample an equally spaced initial seed (this provides polydispersed starting points)
    for ix = 1:sample.numChains
        sample.x0(:,ix)   = cPoints(:,fix(size(cPoints,2)/sample.numChains)*(ix-1) + randi(fix(size(cPoints,2)/sample.numChains)));
    end
end
if strcmp(sample.algorithm,'EDHRB')
    fprintf('%d coordinates spanning the feasible space have been defined.\n',size(sample.udir,2));
else
    fprintf('%d coordinates spanning the feasible space have been defined.\n',size(sample.warmUpPoints,2));
end
sample.prepTime = (cputime-t0)/60;

%% III. Sampling
% Run appropriate sampler in either single core or parallel mode
if ~sample.parallelFlag
    workerIdx = 1;
    if strcmp(sample.algorithm,'EDHRB')
        fprintf('Sampling in progress (single core)...\n---------------------------\n');
        [sample.points,sample.samplingTime] = EDHRB(sample,workerIdx);
    elseif strcmp(sample.algorithm,'ll_ACHRB')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        [sample.points,sample.samplingTime,sample.centroid] = ll_ACHRB(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint,workerIdx);
        if isempty(sample.points); return; end;
    end
else
    % Delete active workers
    if ~isempty(gcp('nocreate')); delete(gcp('nocreate')); end;
    
    % Define number of workers to use. First, try to use the number of workers requested.
    % If not posible, use the number of workers allocated to the 'local' profile.
    if isfield(sample,'numCores')
        try
            parpool(sample.numCores);
        catch
            sample = rmfield(sample,'numCores');
        end
    end
    if ~isfield(sample,'numCores')
        parWorkers = parpool('local');
        sample.numCores = parWorkers.NumWorkers;
    end
    
    % Replicate sample structure for parallel sampling
    numCores  = sample.numCores;
    numChains = sample.numChains;
    samples{numCores}      = [];
    numParticles(numCores) = 0;
    if strcmp(sample.algorithm,'EDHRB')
        ixPrev = 1;
        for ix = 1:sample.numCores
            samples{ix}           = sample;
            numParticles(ix)      = round(numChains/numCores);
            samples{ix}.numChains = numParticles(ix);
            samples{ix}.x0        = sample.x0(:,ixPrev:sum(numParticles));
            numChains             = numChains-numParticles(ix);
            numCores              = numCores-1;
            ixPrev                = sum(numParticles)+1;
        end
        
        % Perform parallel sampling (EDHRB)
        fprintf('Sampling in progress (multiple cores)...\n---------------------------\n');
        points{sample.numCores}       = [];
        samplingTime{sample.numCores} = 0;
        parfor workerIdx = 1:sample.numCores
            rng('shuffle');
            [points{workerIdx},samplingTime{workerIdx}] = EDHRB(samples{workerIdx},workerIdx);
        end
        delete(gcp('nocreate'));
        
        % Build definitive structure
        sample.points       = zeros(sample.pointsPerChain,sample.numRxns,sample.numChains);
        sample.samplingTime = 0;
        countIdx = 1;
        for ix = 1:sample.numCores
            for jx = 1:samples{ix}.numChains
                sample.points(:,:,countIdx) = reshape(points{ix}(:,jx,:),sample.numRxns,sample.pointsPerChain)';
                countIdx = countIdx + 1;
            end
            sample.samplingTime = max([sample.samplingTime,samplingTime{ix}]);
        end
        
    elseif strcmp(sample.algorithm,'ll_ACHRB')
        for ix = 1:sample.numCores
            samples{ix}.points       = [];
            samples{ix}.samplingTime = 0;
            numParticles(ix)         = round(numSamples/numCores);
            numSamples               = numSamples-numParticles(ix);
            numCores                 = numCores-1;
        end
        
        % Perform parallel sampling (ll-ACHRB)
        fprintf('Sampling in progress (multiple cores)...\n--------------------------------------------\n');
        samplingTime{sample.numCores} = 0;
        centroid{sample.numCores}     = [];
        points{sample.numCores}       = [];
        parfor workerIdx = 1:sample.numCores
            rng('shuffle');
            [points{workerIdx},samplingTime{workerIdx},centroid{workerIdx}] = ll_ACHRB(sample,max(numParticles),sample.numDiscarded,sample.stepsPerPoint,1,workerIdx);
        end
        delete(gcp('nocreate'));
        
        % Build definitive structure
        sample.samplingTime = 0;
        sample.points       = [];
        sample.centroid     = [];
        for ixCore = 1:sample.numCores
            sample.points       = [sample.points,points{ixCore}];
            sample.samplingTime = max([sample.samplingTime,samplingTime{ixCore}]);
            sample.centroid     = [sample.centroid,centroid{ixCore}];
        end
    end
end

%%  IV. Run sampling and MCMC diagnostics
% Run convergence and mixing diagnostics
fprintf('Starting MCMC diagnostics...\n');
if ~sample.parallelFlag
    
    % Append independent chains (EDHRB)
    if strcmp(sample.algorithm,'EDHRB'); sample.points = reshape(sample.points,sample.numRxns,sample.numChains*sample.pointsPerChain); end;
    
    % Split samples resulting chain in 10 segments for convergence analysis
    sample.points = reshape(sample.points,sample.numRxns,fix(sample.numSamples/10),10);
    sample.points = permute(sample.points,[2,1,3]);
end

% Calculate potential scale reduction statistics and recover original chain
[sample.R,sample.Rint,sample.Neff,sample.tau,sample.thin] = psrf(sample.points);
sample.points = permute(sample.points,[2,3,1]);
sample.points = reshape(sample.points,sample.numRxns,sample.numSamples);

% Calculate chains statistics
sample.mu    = mean(sample.points,2);
sample.sigma = std(sample.points')';

% Check if there are exclusive topological modes
if (sample.loopless)
    
    % Determine unique flux patterns and continue sampling (if necessary)
    fluxPattern = (~sample.sigma)&(~sample.mu);
    if any(fluxPattern)
        uniquePatterns = unique(fluxPattern','rows');
        fprintf('%d blocked reaction(s) detected in the final sample.\n',size(uniquePatterns,1));
    end
end
