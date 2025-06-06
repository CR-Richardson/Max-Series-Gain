function [alpha,data,dec]=Popov_Like1(syst,eta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% MC Turner and CR Richardson 
% ECS
% University of Southampton
% UK
%
% Date: 15/05/23
%
% Purpose: 
% Compute the maximum series gain (alpha) when using the relaxed Popov-like
% criterion as defined by Theorem 2 and Corollary 1.
%
% Parameters:
% syst: Structure containing the system matrices of an example.
% eta:  Chosen value of relaxation where W = eta*I
%
% Returns:
% alpha: Maximum series gain (float)
% data:  Structure containing solutions of the LMI parametrised by alpha
% dec:   # number of decision variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
A     = syst.a;
B     = syst.b;
C     = syst.c;
D     = syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

%% Initialising alpha

if m == 1
   Gm = margin(ss(A,B,-C,-D));
   if Gm > 10000
      Gm = 10000;
   end
else
    Gm = 100000;
end

% Determine initial upper/lower bound and initial test value
alpha_up  = Gm*0.999;
alpha_low = 0; % We know alpha = 0 is always feasible as system's are stable
alpha     = alpha_up;

eps = 1e-6;

%%
% Determine alpha by repeatedly solving LMI until the largest alpha is 
% found where LMI is feasible 

while ((alpha_up - alpha_low)/alpha_up) > eps

setlmis([]);

W = eta*eye(m);

P    = lmivar(1,[n,1]);
V    = lmivar(2,[m,m]);
H    = lmivar(2,[m,m]);
Q11_ = lmivar(2,[m,m]);
Q11  = lmivar(2,[m,m]);

% LMI
lmiterm([1,1,1,P],A',1,'s');

lmiterm([1,1,2,P],1*alpha,B);
lmiterm([1,1,2,-V],C',1);
lmiterm([1,1,2,-H],A'*C',1);

lmiterm([1,1,3,-H],A'*C',1);
lmiterm([1,1,3,-H],C',W);
lmiterm([1,1,3,0],-C'*W);

lmiterm([1,2,2,V],-1,1,'s');
lmiterm([1,2,2,H],1,C*B*alpha,'s');
lmiterm([1,2,2,Q11_],1,1,'s');
lmiterm([1,2,2,Q11],1,1,'s');

lmiterm([1,2,3,-H],alpha*B'*C',1);
lmiterm([1,2,3,Q11_],1,1);

lmiterm([1,3,3,0],-2*W);

% P > 0
lmiterm([2,1,1,P],-1,1);

% V: Z-matrix conditions
 JJ = eye(m);
 count = 3;       
 for i = 1:m
     for j = 1:m
         e1 = JJ(i,:); e2 = JJ(:,j);
         if i ~= j
            lmiterm([count,1,1,V],e1,e2,'s');
            count = count+1;
         end
     end
 end
        
% Q11_: Positivity matrix conditions
 for i = 1:m
     for j = 1:m
         e1 = JJ(i,:); e2 = JJ(:,j);
         lmiterm([count,1,1,Q11_],-e1,e2,'s');
         count = count+1;
     end
 end

 % Q11: Positivity matrix conditions
for i = 1:m
    for j = 1:m
        e1 = JJ(i,:); e2 = JJ(:,j);
        lmiterm([count,1,1,Q11],-e1,e2,'s');
        count = count+1;
     end
end 

LMISYS = getlmis;
[tmin,xfeas] = feasp(LMISYS,[1e-20 5000 -0.1 1000 1]);

% Update alpha upper/lower bound plus new test value
 if tmin < 0  % if LMIs are feasible
    alpha_low = alpha;
  else 
    alpha_up = alpha; % if LMIs are infeasible
  end
  
alpha = (alpha_up + alpha_low)/2;
end

%% Return solutions
dec        =  decnbr(LMISYS); % returns number of decision varibles
data.P     =  dec2mat(LMISYS,xfeas,P);
data.V     =  dec2mat(LMISYS,xfeas,V);
data.H     =  dec2mat(LMISYS,xfeas,H);
data.Q11_  =  dec2mat(LMISYS,xfeas,Q11_);
data.Q11   =  dec2mat(LMISYS,xfeas,Q11);

if D ~= zeros(m)
    disp('D not equal to zero. Popov-like criterion may not be applied!');
    alpha      = nan;
    dec        =  nan;
    data.P     =  nan;
    data.V     =  nan;
    data.H     =  nan;
    data.Q11_  =  nan;
    data.Q11   =  nan;
end

end

