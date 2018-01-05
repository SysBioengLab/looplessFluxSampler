function points = bringToBoundary(points,lb,ub)
% Bring points to boundary
%
% Bring points ouside the boundaries back to the feasible region
%
% USAGE:
%              points = bringToBoundary(points, lb, ub)
%
% INPUTS:
%              points:  n x k Matrix of points
%              lb:      n x 1 Lower bounds
%              ub:      n x 1 Upper bounds
%
% OUTPUT:
%              points:  n x p Matrix of points inside the feasible region
%
% -------------------- Copyright (C) 2018 Pedro A. Saa --------------------

points = bsxfun(@times,bsxfun(@gt,points,ub),ub) + bsxfun(@times,~bsxfun(@gt,points,ub),points);
points = bsxfun(@times,bsxfun(@lt,points,lb),lb) + bsxfun(@times,~bsxfun(@lt,points,lb),points);