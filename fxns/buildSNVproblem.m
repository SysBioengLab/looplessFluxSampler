function SNVproblem = buildSNVproblem(S,lb,ub,x,epsilon)
% Builds Sparse Vector Pursuit (SNV) problem
%
% Structure has the following form
%
%    min      ~& ||v||_1 \\
%    s.t.     ~& S v = 0 \\
%             ~& x v \geq \epsilon \\
%             ~& or x v \leq -\epsilon \\
%             ~& lb \leq v \leq ub \\
%
% USAGE:
%              SNVproblem = buildSNVproblem(S,lb,ub,x,epsilon)
%
% INPUTS:
%              S:  m x 1 Stoichiometric matrix
%              lb: n x 1 Lower bounds
%              ub: n x 1 Upper bounds
%              x:  1 x n Projection matrix onto the null space generated by the previous sparse vectors
%              epsilon:  Small positive constant
%
% OUTPUT:
%              SNVproblem (structure):   SVP problem structure
%
% -------------------- Copyright (C) 2019 Pedro A. Saa --------------------

% Parameter initialization
[m,n] = size(S);

% Build the appropriate LP structure (minimal and maximal)
SNVproblem.A = sparse([S,zeros(m,n);...
    -eye(n),eye(n);...
    eye(n),eye(n);...
    x,zeros(1,n)]);

% Formulate rhs
SNVproblem.b      = zeros(m+2*n+1,1);
SNVproblem.b(end) = epsilon;                                               % This ensures that the solution is linearly independent from the previous solutions

% Define constraint sense
SNVproblem.csense = blanks(m+2*n+1);
for ix = 1:m+2*n+1
    if ix<=m
        SNVproblem.csense(ix) = 'E';
    else
        SNVproblem.csense(ix) = 'G';
    end
end

% Define objective, osense and bounds
SNVproblem.c      = [zeros(n,1);ones(n,1)];
SNVproblem.osense = 1;                                                     % Objective sense: minimization (required to eliminate infeasible loops)
SNVproblem.lb     = [lb;zeros(n,1)];
SNVproblem.ub     = [ub;1e2*ones(n,1)];