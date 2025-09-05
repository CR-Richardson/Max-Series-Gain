%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% CR Richardson and MC Turner
% ECS
% University of Southampton
% UK
%
% Date: 5/09/25
%
% Purpose: 
% Compute the maximum series gain (alpha) when using the Circle-like 
% criterion as defined by Theorem 1.
%
% Parameters:
% syst:      Structure containing the system matrices of an example.
% eps:       Loop termination accuracy.
% alpha_up:  Largest series gain to be tested.
% alpha_low: Lowest series gain to be tested.
%
% Returns:
% alpha_low: Maximum series gain.
% data:      Structure containing solutions of the LMI parametrised by alpha.
% dec:       # of decision variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha_low,data,dec]=Circle_Like(syst,eps,alpha_up,alpha_low)

%% Parameters
A     = syst.a;
B     = syst.b;
C     = syst.c;
D     = syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

%% Initialising alpha
alpha = alpha_up;

%% Determine max alpha by repeatedly solving LMI for largest feasible alpha 

while ((alpha_up-alpha_low)/alpha_up) > eps

% define variables
P = sdpvar(n,n,'symmetric');
V = sdpvar(m,m,'full');
Q11 = sdpvar(m,m,'full');

% define constraints
epsilon = 1e-3;
In = eye(n);
Im = eye(m);
c1 = P >= epsilon*In; % PD
c2 = V(~Im) <= 0; % Z-matrix
c3 = Q11(:) >= 0; % Non-negative matrix
c4 = [A'*P + P*A, ...
      alpha*P*B + C'*V';
      alpha*B'*P + V*C, ...
      Q11 + Q11' - V - V'] <= -epsilon*eye(n+m); %ND

Constraints = [c1, c2, c3, c4];

% solve LMI
Objective = [];
options = sdpsettings('solver', 'mosek', 'verbose', 0);
result = optimize(Constraints, Objective, options);

% Update alpha upper/lower bound plus new test value
if result.problem == 0 % if LMIs are feasible
    alpha_low = alpha;
else 
    alpha_up = alpha; % if LMIs are infeasible
end
  
alpha = (alpha_up + alpha_low)/2;
end

%% output

% Get variable indices
vars = [P(:); V(:); Q11(:)];
dec = length(unique(getvariables(vars)));

% solutions
data.P = value(P);
data.V = value(V);
data.Q11 = value(Q11);

end