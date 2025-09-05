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
% Compute the maximum series gain (alpha) when using the Popov-like criterion
% as defined by Theorem 2 and Corollary 1.
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

function [alpha_low,data,dec]=Popov_Like1(syst,eps,alpha_up,alpha_low)

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

while ((alpha_up - alpha_low)/alpha_up) > eps

eta = 10000;
W = eta*eye(m);

% define variables
P    = sdpvar(n,n,'symmetric');
V    = sdpvar(m,m,'symmetric'); % full
H    = sdpvar(m,m,'symmetric'); % full
Q11  = sdpvar(m,m,'symmetric'); % full
Q11_ = sdpvar(m,m,'symmetric'); % full


% define constraints
epsilon = 1e-3;
In = eye(n);
Im = eye(m);
c1 = P >= epsilon*In; % PD
c2 = V(~Im) <= 0; % Z-matrix
c3 = Q11(:) >= 0; % Non-negative matrix
c4 = Q11_(:) >= 0; % Non-negative matrix

X11 = A'*P + P*A;
X12 = alpha*P*B + C'*V' + A'*C'*H';
X13 = A'*C'*H' + C'*(H' - Im)*W;
X22 = (Q11_ + Q11 - V + alpha*H*C*B) + (Q11_ + Q11 - V + alpha*H*C*B)';
X23 = alpha*B'*C'*H' + Q11_;
X33 = -2*W;
c5 = [X11 , X12 , X13;
      X12', X22 , X23;
      X13', X23', X33] <= -epsilon*eye(n+2*m); % ND
 
Constraints = [c1, c2, c3, c4, c5];
 
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
vars = [P(:); V(:); H(:); Q11(:); Q11_(:)];
dec = length(unique(getvariables(vars)));
 
% solutions
data.P = value(P);
data.V = value(V);
data.H = value(H);
data.Q11 = value(Q11);
data.Q11_ = value(Q11_);

end