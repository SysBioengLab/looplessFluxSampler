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
%                                    * numDiscarded - Burn-in (double) (only used in ll-ACHRB) (default 0)
%                                    * stepsPerPoint - Thinning or number of steps per effective point (double). In ADSB, it refers to the expected number of times a point is moved 
%                                    * algorithm - 'ADSB' (default) or 'll_ACHRB'
%                                    * loopless - Loop removal flag, true (default) or false
%                                    * vTol - Numerical flux tolerance (default 1e-8)
%                                    * parallelFlag - Parallel sampling option, false (default) or true
%                                    * numCores - Number of cores for parallel sampling (double), empty (default)                               
%                                    * diagnostics - MCMC diagnostics flag, true (default) or false
%                                    * populationScale - Size of current set relative to the feasible space dimension (integer => 1), 3 (default)
%                                    * points - Points from a previous sampling, empty (default)
%
% OUTPUT:
%              sample (structure):   sampling structure containing #numSamples random points
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

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
    else fprintf('Field numSamples has not been defined.'); return; end

    % Samples discarded or burn-in (default 0)
    if isfield(options,'numDiscarded'); sample.numDiscarded = options.numDiscarded;
    else sample.numDiscarded = 0; end

    % Steps per point (thinning in ACHRB) or number of expected moves (ADSB) (default 1e2)
    if isfield(options,'stepsPerPoint'); sample.stepsPerPoint = options.stepsPerPoint;
    else sample.stepsPerPoint = 1e2; end

    % Sampling algorithm 'ADSB' (default) or 'll_ACHRB'
    if isfield(options,'algorithm'); sample.algorithm = options.algorithm;
    else sample.algorithm = 'ADSB'; end

    % Sampling mode (default loopless true or 1)
    if isfield(options,'loopless'); sample.loopless = options.loopless;
    else sample.loopless = 1; end

    % Numerical flux tolerance  (default 1e-8)
    if isfield(options,'vTol'); sample.vTol = options.vTol;
    else sample.vTol = 1e-8; end
    
    % Parallel sampling option (default false or 0)
    if isfield(options,'parallelFlag'); sample.parallelFlag = options.parallelFlag;
    else sample.parallelFlag = 0; end
    
    % Number of cores for parallel sampling (default empty)
    if isfield(options,'numCores'); sample.numCores = options.numCores;
    else sample.numCores = []; end
    
    % Run MCMC diagnostics
    if isfield(options,'diagnostics'); sample.diagnostics = options.diagnostics;
    else sample.diagnostics = 1; end
    
    % Fold number of particles in the population relative to dim(Omega) (only for ADSB, default 1) 
    if strcmp(sample.algorithm,'ADSB') && isfield(options,'populationScale'); sample.populationScale = options.populationScale;
    else sample.populationScale = 3; end
    
    % Restart sampler from previous results
    if ~isfield(model,'points'); sample.points = [];
    else sample.points = model.points; end
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

        % If there are no possible active loops, then turn-off loopless option and
        % run simple FVA, not ll-FVA
    else
        sample.loopless = 0;
        LPproblem       = buildLinearProblem(tempModel);
        [tempModel.lb,tempModel.ub] = generalFVA(LPproblem,sample.vTol);
    end
    [tempModel,zeroRxns] = removeBlockedSets(tempModel);
    rxnList(zeroRxns)    = [];
    [tempModel,rxnList]  = parseInternalRxns(tempModel,rxnList);

    % Get an updated sparse null-space matrix of internal rxns (if relevant)
    if sample.loopless
        Nint = fastSNP(tempModel.S(:,tempModel.intRxns),-1e2*(tempModel.lb(tempModel.intRxns)<0),...
            1e2*(tempModel.ub(tempModel.intRxns)>0),sample.vTol);
        if ~isempty(Nint); Nint = normMatrixEntries(Nint);
        else sample.loopless = 0; end
    end
    if sample.loopless; fprintf('The model has %d infeasible loop law(s) potentially active.\n',size(Nint,2));
    else fprintf('The model does not contain infeasible loops active.\n'); end
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

% Generate seeds for ADSB
if strcmp(sample.algorithm,'ADSB')
    
    % Generate seed points for the current set: Run ll-ACHRB
    if isempty(sample.points)
        nWarmups   = 2e4;
        maxPoints  = 2e6;
        cWeight    = size(sample.warmUpPoints,2);
        cCentroid  = [];
        prevPoints = [];
        while (nWarmups < maxPoints)
            [vPoints,~,sample.centroid,cWeight] = ll_ACHRB(sample,nWarmups,0,1,1,cWeight);
            if isempty(vPoints); return; end
            
            % Calculate centroid trace evolution
            prevPoints = [prevPoints,vPoints];
            cCentroid  = [cCentroid,bsxfun(@rdivide,cumsum(vPoints,2),size(cCentroid,2)+1:size(cCentroid,2)+size(vPoints,2))];
            Rfactor    = psrf(cCentroid');
            
            % Check the convergence of the reduction factor
            if (abs(median(Rfactor(isfinite(Rfactor)))-1) <= .1); break;
            else nWarmups = 2*nWarmups; end
        end
    else
        % Use previous points in the current set
        prevPoints = sample.points;
    end
    
    % Determine size of the feasible space
    omegaSize = rank(prevPoints(:,randsample(size(prevPoints,2),sample.numRxns,false)));
    if (sample.populationScale == 1); nDim = sample.populationScale*omegaSize+1;
    else nDim = sample.populationScale*omegaSize; end
    
    % Define number of chains and points per chain (at least three points are required)
    minPoints             = 3;
    sample.pointsPerChain = max([nDim,minPoints]);
    sample.numChains      = ceil(sample.numSamples/sample.pointsPerChain);
    
    % Define current sets
    if isempty(sample.points)
        sample.points = zeros(sample.numRxns,nDim,sample.numChains);
        for ix = 1:sample.numChains
            sample.points(:,:,ix) = prevPoints(:,randsample(size(prevPoints,2),nDim,'false'));
        end
    else        
        extraPoints = sample.pointsPerChain*sample.numChains-sample.numSamples;
        if (extraPoints ~= 0)
            extraPoints = prevPoints(:,randsample(size(prevPoints,2),extraPoints,'false'));
            prevPoints   = [prevPoints,extraPoints];
        end
        prevPoints = reshape(prevPoints,sample.numRxns,nDim,sample.numChains);
        sample.points = prevPoints;
    end
end

if strcmp(sample.algorithm,'ADSB')
    fprintf('%d points spanning the feasible space have been incorporated to the current set.\n',sample.pointsPerChain);
else
    fprintf('%d warmup points spanning the feasible space have been defined.\n',size(sample.warmUpPoints,2));
end
sample.prepTime = (cputime-t0)/60;
clearvars -except sample

%% III. Sampling
% Run appropriate sampler in either single core or parallel mode
if ~sample.parallelFlag
    if strcmp(sample.algorithm,'ADSB')
        fprintf('Sampling in progress (single core)...\n---------------------------\n');
        [sample.points,sample.samplingTime] = ADSB(sample);
    elseif strcmp(sample.algorithm,'ll_ACHRB')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        [sample.points,sample.samplingTime,sample.centroid] = ll_ACHRB(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint);
        if isempty(sample.points); return; end
    end
else
    % Delete active workers
    if ~isempty(gcp('nocreate')); delete(gcp('nocreate')); end

    % Define number of workers to use. First, try to use the number of workers requested.
    % If not posible, use the number of workers allocated to the 'local' profile
    if ~isempty(sample.numCores)
        try
            parpool(sample.numCores);
        catch
            parWorkers = parpool('local');
            sample.numCores = parWorkers.NumWorkers;
        end
    else
        parWorkers = parpool('local');
        sample.numCores = parWorkers.NumWorkers;
    end

    % Choose appropriate sampling algorithm
    if strcmp(sample.algorithm,'ADSB')
        
        % Replicate sample structure for parallel sampling
        for ix = 1: sample.numChains
           samples{ix}           = sample;
           samples{ix}.numChains = 1;
           samples{ix}.points    = sample.points(:,:,ix);
        end

        % Perform parallel sampling using ADSB
        fprintf('Sampling in progress (multiple cores)...\n');
        points{sample.numChains}       = [];
        samplingTime{sample.numChains} = 0;
        parfor workerIdx = 1:sample.numChains
            rng('shuffle');
            [points{workerIdx},samplingTime{workerIdx}] = ADSB(samples{workerIdx},0);
        end
        fprintf('Sampling finalized.\n');
        delete(gcp('nocreate'));
        
        % Build definitive structure
        sample.points       = zeros(sample.numRxns,sample.pointsPerChain,sample.numChains);
        sample.samplingTime = 0;
        for ix = 1:sample.numChains            
            sample.points(:,:,ix) = points{ix}(:,:,1);
            sample.samplingTime   = max([sample.samplingTime,samplingTime{ix}]);
        end

    elseif strcmp(sample.algorithm,'ll_ACHRB')
        
        % Replicate sample structure for parallel sampling
        numCores   = sample.numCores;
        numSamples = sample.numSamples;
        samples{numCores}      = [];
        numParticles(numCores) = 0;
        for ix = 1:sample.numCores
            samples{ix}.points       = [];
            samples{ix}.samplingTime = 0;
            numParticles(ix)         = round(numSamples/numCores);
            numSamples               = numSamples-numParticles(ix);
            numCores                 = numCores-1;
        end

        % Perform parallel sampling using ll-ACHRB
        fprintf('Sampling in progress (multiple cores)...\n--------------------------------------------\n');
        samplingTime{sample.numCores} = 0;
        centroid{sample.numCores}     = [];
        points{sample.numCores}       = [];
        parfor workerIdx = 1:sample.numCores
            rng('shuffle');
            [points{workerIdx},samplingTime{workerIdx},centroid{workerIdx}] = ll_ACHRB(sample,max(numParticles),sample.numDiscarded,sample.stepsPerPoint,1,workerIdx);
        end
        delete(gcp('nocreate'));

        % Build final sample structure
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

% Remove auxiliary fields
sample = rmfield(sample,'warmUpPoints');

%%  IV. Run sampling and MCMC diagnostics
if sample.diagnostics
    fprintf('Starting MCMC diagnostics...\n');
    if ~strcmp(sample.algorithm,'ADSB')
        
        % Split samples resulting from the chain in 10 segments for
        % convergence analysis of ll-ACHRB (this is not necessary for ADSB)
        sample.numChains = 10;
        sample.pointsPerChain = fix(sample.numSamples/sample.numChains);
        sample.points = reshape(sample.points,sample.numRxns,sample.pointsPerChain,sample.numChains);
    end
    
    % Calculate potential scale reduction statistics and recover original chain
    sample.points = permute(sample.points,[2,1,3]);
    [sample.R,sample.Rint,sample.Neff,sample.tau,sample.thin] = psrf(sample.points);
    sample.points = permute(sample.points,[2,3,1]);
    sample.points = reshape(sample.points,sample.numRxns,sample.numChains*sample.pointsPerChain);
    sample.points = sample.points(:,1:sample.numSamples);
    
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
else
    if strcmp(sample.algorithm,'ADSB')
        sample.points = reshape(sample.points,sample.numRxns,sample.numChains*sample.pointsPerChain);
    end
    sample.points = sample.points(:,1:sample.numSamples);
end
