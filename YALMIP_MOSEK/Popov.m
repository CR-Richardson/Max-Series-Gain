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
% Compute the maximum series gain (alpha) when using the Popov criterion
% as defined by Theorem 2 and Remark 4.
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

function [alpha_low,data,dec]=Popov(syst,eps,alpha_up,alpha_low)

%% Parameters
A     = syst.a;
B     = syst.b;
C     = syst.c;
D     = syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

% Raise error if D != 0.
if ~isequal(D, zeros(m))
    error('D is not a zero matrix of size %d x %d.', m, m);
end

%% Initialising alpha
alpha = alpha_up;

%% Determine max alpha by repeatedly solving LMI for largest feasible alpha 

while ((alpha_up-alpha_low)/alpha_up) > eps

% define variables
P = sdpvar(n,n,'symmetric');
L = sdpvar(m,m,'diagonal');
V = sdpvar(m,m,'diagonal');

% define constraints
epsilon = 1e-3;
In = eye(n);
Im = eye(m);
c1 = P >= epsilon*In; % PD
c2 = L >= epsilon*Im; % PD
c3 = V >= epsilon*Im; % PD

X11 = P*A + A'*P;
X12 = alpha*P*B + C'*V + A'*C'*L;
X22 = alpha*L*C*B + alpha*B'*C'*L - 2*V;
c4 = [X11 , X12;
      X12', X22] <= -epsilon*eye(n+m); % ND


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
vars = [P(:); L(:); V(:)];
dec = length(unique(getvariables(vars)));

% solutions
data.P  =  value(P);
data.L  =  value(L);
data.V  =  value(V);

end