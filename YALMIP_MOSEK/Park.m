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
% Compute the maximum series gain (alpha) when using the Park criterion
% as defined by Theorem 2 from his 2002 IEEETAC paper.
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

function [alpha_low,data,dec]=Park(syst,eps,alpha_up,alpha_low)

%% Parameters
A     = syst.a;
B     = -syst.b;
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

% define variables
X  = sdpvar(n,n,'symmetric');
Y  = sdpvar(n,m,'full');
Z  = sdpvar(m,m,'symmetric');
M  = sdpvar(m,m,'diagonal');
N1 = sdpvar(m,m,'diagonal');
N2 = sdpvar(m,m,'diagonal');

% define constraints
epsilon = 1e-3;
c1 = [X,  Y;
      Y', Z] >= 0; % PSD
c2 = N1 >= 0; % PSD
c3 = N2 >= 0; % PSD

X11 = X*A + A'*X;
X12 = -alpha*X*B + A'*Y + A'*C'*M + C'*N1;
X13 = Y + A'*C'*N2;
X22 = -(Y' + M*C)*alpha*B - alpha*B'*(Y + C'*M) - 2*N1;
X23 = Z - alpha*B'*C'*N2;
X33 = -2*N2;

c4 = [X11 , X12 , X13;
      X12', X22 , X23;
      X13', X23', X33] <= -epsilon*eye(n+2*m); % ND

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
vars = [X(:); Y(:); Z(:); M(:); N1(:); N2(:)];
dec = length(unique(getvariables(vars)));

% solutions
data.X  =  value(X);
data.Y  =  value(Y);
data.Z  =  value(Z);
data.M  =  value(M);
data.N1 =  value(N1);
data.N2 =  value(N2);

end