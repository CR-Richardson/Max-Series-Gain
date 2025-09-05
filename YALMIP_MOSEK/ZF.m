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
% Compute the maximum series gain (alpha) when using the Zames-Falb criterion
% as defined in Reference 30.
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

function [alpha_low,data,dec]=ZF(syst,eps,alpha_up,alpha_low)

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

% define variables
P11  = sdpvar(n,n,'symmetric');
S11  = sdpvar(n,n,'symmetric');
N    = sdpvar(n,n,'symmetric');
H0   = sdpvar(m,m,'diagonal');
etaP = sdpvar(m,m,'diagonal');
VP   = sdpvar(m,m,'diagonal');
Acb = sdpvar(n,n,'symmetric');
Aab = sdpvar(n,n,'symmetric');
Bcb_vec = sdpvar(n,1);
Bab_vec = sdpvar(n,1);
Bcb = repmat(Bcb_vec, 1, m); % nxm matrix with identical columns.
Bab = repmat(Bab_vec, 1, m); % nxm matrix with identical columns.

% define constraints
epsilon = 1e-3;
In = eye(n);
Im = eye(m);
WB = Im;
WC = Im;

% define constraints: positive definite conditions
c1 = S11 >= epsilon*In;
c2 = N >= epsilon*In;
c3 = P11 - S11 - N >= epsilon*In;
c4 = VP >= epsilon*Im;

% define constraints: L1 bound LMI
e1   = Im(:,1);
X11  = -1;
X12  = -e1'*Bcb'; % 1xn
X13  = -e1'*Bab'; % 1xn
X22  = -Acb;
X23  = zeros(n,n);
X33  = -Aab;
c5   = [X11 , X12 , X13; 
        X12', X22 , X23;
        X13', X23', X33] <= -epsilon*eye(1+2*n);

% define constraints: LMI
Y11  = S11*A + A'*S11 - Bab*WC*C - C'*WC'*Bab';
Y12  = S11*A + A'*P11 + Acb' + Aab' + C'*WB*Bcb' - Bab*WC*C;
Y13  = A'*N + C'*WC*Bab' + Aab';
Y14  = alpha*S11*B + Bcb*WC + Bab*WC + C'*H0' + A'*C'*etaP + C'*VP';
Y22  = P11*A + A'*P11 + Bcb*WB*C + C'*WB'*Bcb';
Y23  = A'*N + C'*WC*Bab' - Aab;
Y24  = alpha*P11*B - Bcb*WB - Bab*WB + C'*H0' + A'*C'*etaP + C'*VP';
Y33  = -2*Aab;
Y34  = alpha*N*B - Bab*(WB + WC);
Y44  = -2*H0 + alpha*etaP*C*B + alpha*B'*C'*etaP' - 2*VP;
c6   = [Y11 , Y12 , Y13 , Y14;
        Y12', Y22 , Y23 , Y24;
        Y13', Y23', Y33 , Y34;
        Y14', Y24', Y34', Y44] <= -epsilon*eye(m+3*n); % ND

Constraints = [c1, c2, c3, c4, c5, c6];

% define constraints: L1 bound LMIs
for i=1:m
    Constraints = [Constraints, trace(WB)*Im(i,:)*WC*Im(:,i) - Im(i,:)*H0*Im(:,i) <= 0];
end

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
vars = [P11(:); S11(:); N(:); H0(:); Acb(:); Aab(:); Bcb(:); Bab(:); etaP(:); VP(:)];
dec = length(unique(getvariables(vars)));
      
% solutions
data.P11 =  value(P11);
data.S11 =  value(S11);
data.N   =  value(N);
data.H0  =  value(H0);
data.Acb =  value(Acb);
data.Aab =  value(Aab);
data.Bcb =  value(Bcb);
data.Bab =  value(Bab);
data.etaP =  value(etaP);
data.VP   =  value(VP);

end