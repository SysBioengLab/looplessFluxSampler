function [points,samplingTime] = EDHRB(sample,workerIdx)
% Elongated Directions Hit and Run on a Box
%
% Uses EDHRB to generate a uniform random sample from the loopless flux
% solution space
%
% USAGE:
%              points = EDHRB(sample)
%
% INPUTS:
%              sample (structure):    (the following fields are required - others can be supplied)
%                                    * pointsPerChain - Number of points for the chains (double)
%                                    * numDiscarded - Burn-in (double)
%                                    * stepsPerPoint - Thining (double)
%                                    * loopless - Loop removal option, true (default) or false
%                                    * stepFxn - Step function option,'linear' (default) or 'uniform'
%                                    * numChains - Number of Markov chains (double) 
%                                    * x0  - `n x M` Matrix of initial points for #numChains Markov chains
%                                    * udir - `n x K` Matrix with feasible directions of movement
%                                    * keepWithinBounds - Fxn handle for bringToBoundary.m fxn
%                                    * isFeasible - Fxn handle for looplessCheck.m fxn
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%
% OPTIONAL INPUTS:
%              workerIdx:  id of the working processor, 1 (default)
%
% OUTPUT:
%              points:   Matrix with n x #numChains x #pointsPerChains (loopless) flux solutions
%
% OPTIONAL OUTPUT:
%              samplingTime:   Runtime of EDHRB
%
% -------------------- Copyright (C) 2018 Pedro A. Saa --------------------

% Check input arguments
if nargin<2; workerIdx = 1; end;

% Define mapping fxn for new points
mappingFxn = @(x_0,udir,step) x_0 + bsxfun(@times,udir,step);

% Allocating storage for particles from (numChains) iterations
points = repmat(sample.x0,[1,1,sample.pointsPerChain]);

% Build step-calculating function
if strcmp(sample.stepFxn,'linear')
    stepSizeFxn = @(x_0,L,R) L.*(1 + (-2*x_0./L - 1) - sqrt( (1 + (-2*x_0./L - 1)).^2 - 4*(-2*x_0./L - 1).*R ))./(2*(-2*x_0./L - 1)) + x_0;
elseif strcmp(sample.stepFxn,'uniform')
    stepSizeFxn = @(x_0,L,R) x_0 + L.*R;
end

% Start clock and sampler
currPoint  = sample.x0;
nextPoint  = zeros(size(currPoint));
iter       = 0;
iterSample = 0;
t0         = cputime;                                         % Start clock
if (workerIdx==1); fprintf('%%Prog \t Time \t Time left\n---------------------------\n'); end;
while (iterSample<sample.pointsPerChain)
    
    % Print progress information every 1e4 steps
    if (iterSample>0) && ~mod(iter-1,1e4) && (workerIdx==1)
        timeElapsed = (cputime-t0)/60;
        timePerStep = timeElapsed/iterSample;
        fprintf('%d\t%8.2f\t%8.2f\n',round(1e2*iterSample/sample.pointsPerChain),timeElapsed,(sample.pointsPerChain-iterSample)*timePerStep);        
    end
           
    % Figure out maximum distance to boundaries
    distUb = bsxfun(@plus,-currPoint,sample.ub);
    distLb = bsxfun(@minus,currPoint,sample.lb);
        
    % Pre-allocate variables
    flags  = true(1,sample.numChains);
    
    % Check distance to the boundaries based on chosen direction and correct if necessary
    while any(flags)
        
        % Get random direction
        idxUdir(flags) = randi(size(sample.udir,2),1,sum(flags));
        udir(:,flags)  = sample.udir(:,idxUdir(flags));
        
        % Figure out positive and negative directions (only for directions )
        posDir(:,flags)      = bsxfun(@gt,udir(:,flags),sample.uTol);
        negDir(:,flags)      = bsxfun(@lt,udir(:,flags),-sample.uTol);
        posStepTemp(:,flags) = bsxfun(@rdivide,distUb(:,flags),udir(:,flags));
        negStepTemp(:,flags) = bsxfun(@rdivide,-distLb(:,flags),udir(:,flags));

        % Figure out the true max & min step sizes
        for ix = find(flags)
            posStep(ix) = min([posStepTemp(posDir(:,ix),ix);negStepTemp(negDir(:,ix),ix)]);
            negStep(ix) = max([negStepTemp(posDir(:,ix),ix);posStepTemp(negDir(:,ix),ix)]);
        end

        % Draw new direction if too close to the boundaries
        flags(flags) = (posStep(flags)-negStep(flags) < sample.bTol) | (posStep(flags) < 0) | (negStep(flags) > 0);          
    end

    % Sample from the shrinking interval as in slice sampling
    flags = true(1,sample.numChains);
    Lcord = posStep-negStep;             % Determine length of maximum cord
    while any(flags)
        
        % Perform a random step
        step(flags)        = stepSizeFxn(negStep(flags),Lcord(flags),rand(1,sum(flags)));
        nextPoint(:,flags) = mappingFxn(currPoint(:,flags),udir(:,flags),step(flags));
        nextPoint(:,flags) = sample.keepWithinBounds(nextPoint(:,flags));
        if sample.loopless
            flags(flags)   = ~sample.isFeasible(nextPoint(:,flags));
            
            % Shrink hypercube if necessary
            posSign        = (step(flags)>0);
            posStep(flags) = posStep(flags).*(~posSign) + step(flags).*(posSign);
            negStep(flags) = negStep(flags).*(posSign) + step(flags).*(~posSign);
            
            % Keep shrinking only if the bounding box is not too small
            Lcord(flags) = posStep(flags)-negStep(flags);
            flags(flags) = bsxfun(@gt,Lcord(flags),sample.bTol);
        else
            break;
        end
    end

    % Reset position of the next point to the current point for zero-length steps
    fixed_point = (Lcord<sample.bTol);
    nextPoint(:,fixed_point) = currPoint(:,fixed_point);

    % Save new point every stepsPerPoint samples
    if (iter>sample.numDiscarded) && ~mod(iter,sample.stepsPerPoint)
        points(:,:,iterSample+1) = nextPoint;
        iterSample = iterSample+1;
    end
    
    % Update general counter and point
    iter = iter+1;
    currPoint = nextPoint;
end
samplingTime = (cputime-t0)/60;
if (workerIdx==1); fprintf('---------------------------\nSampling finalized.\n'); end;