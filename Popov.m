function [alpha,data]=Popov(syst)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% MC Turner and CR Richardson 
% ECS
% University of Southampton
% UK
%
% Date: 03/03/23
%
% Purpose: 
% Compute the maximum series gain (alpha) when using the Popov
% criterion as defined by Theorem 2 and Remark 6.
%
% Parameters:
% syst: Structure containing the system matrices of an example.
%
% Returns:
% alpha: Maximum series gain (float)
% data: Structure containing solutions of the LMI parametrised by alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Parameters

A=syst.a;
B=syst.b;
C=syst.c;
D=syst.d;

[n,m]=size(B); % n = dimension of state, m = dimension of output

%%
% Initialising alpha

if m == 1
   Gm=margin(ss(A,B,-C,-D));
   if Gm > 10000
      Gm = 10000;
   end
else
    Gm = 1000;
end

% Determine initial upper/lower bound and inital test value
alpha_up=Gm*0.999;
alpha_low=0; % We know alpha = 0 is always feasible as system's are stable
alpha=alpha_up;

%%
% Determine alpha by repeatedly solving LMI until the largest alpha is 
% found where LMI is feasible

while ((alpha_up-alpha_low)/alpha_up) > 0.0001

setlmis([]);

P = lmivar(1,[n,1]);
V = lmivar(1,kron([1,0],ones(m,1)));
L = lmivar(1,kron([1,0],ones(m,1)));

% LMI
lmiterm([1,1,1,P],A',1,'s');

lmiterm([1,1,2,P],1*alpha,B); 
lmiterm([1,1,2,L],A'*C',1);
lmiterm([1,1,2,-V],C',1);

lmiterm([1,2,2,L],1,C*B*alpha,'s');
lmiterm([1,2,2,V],-1,1,'s');

% P > 0
lmiterm([2,1,1,P],-1,1);

% V > 0
lmiterm([3,1,1,V],-1,1);

% L > 0
lmiterm([4,1,1,L],-1,1);

LMISYS=getlmis;
[tmin,xfeas] = feasp(LMISYS,[1e-20 5000 -0.1 1000 1]);

% Update alpha upper/lower bound plus new test value
 if tmin < 0  % if LMIs are feasible
    alpha_low = alpha;
 else 
    alpha_up = alpha; % if LMIs are infeasible
 end
  
alpha=(alpha_up+alpha_low)/2;
end

%%
% Return solutions
data.P  =  dec2mat(LMISYS,xfeas,P);
data.V  =  dec2mat(LMISYS,xfeas,V);
data.L  =  dec2mat(LMISYS,xfeas,L);

end

