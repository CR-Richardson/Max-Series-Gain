function [alpha,data,dec]=Park(syst)

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
% Compute the maximum series gain (alpha) when using the Park criterion as 
% defined by Theorem 2 from his 2002 IEEETAC paper
%
% Note: Park uses negative feedback
%
% Parameters:
% syst: Structure containing the system matrices of an example.
%
% Returns:
% alpha: Maximum series gain (float)
% data:  Structure containing solutions of the LMI parametrised by alpha
% dec:   # number of decision variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
A     =  syst.a;
B     =  syst.b;
C     = -syst.c;
D     = -syst.d;
[n,m] = size(B); % n = dimension of state, m = dimension of output

%% Initialising alpha

if m == 1
   Gm = margin(ss(A,B,C,D));			   
   if Gm > 10000
      Gm = 10000;
   end
else
    Gm = 1000;
end

% Determine initial upper/lower bound and initial test value
alpha_up  = Gm*0.999;
alpha_low = 0; % We know alpha = 0 is always feasible as system's are stable
alpha     = alpha_up;

%%
% Determine alpha by repeatedly solving LMI until the largest alpha is 
% found where LMI is feasible

while ((alpha_up - alpha_low)/alpha_up) > 0.0001

setlmis([]);

X = lmivar(1,[n,1]);
Y = lmivar(2,[n,m]);
Z = lmivar(1,[m,1]);
clear structM1 structM2;
for i=1:m
    structM1(i) = 1;
    structM2(i) = 0;
end
[M,~,~]  = lmivar(1,[structM1' structM2']);
[N1,~,~] = lmivar(1,[structM1' structM2']);
[N2,~,~] = lmivar(1,[structM1' structM2']);

% LMI
lmiterm([1,1,1,X],1,A,'s');

lmiterm([1,1,2,X],-1*alpha,B);
lmiterm([1,1,2,Y],A',1);
lmiterm([1,1,2,M],A'*C',1);
lmiterm([1,1,2,N1],C',1);

lmiterm([1,1,3,Y],1,1);
lmiterm([1,1,3,N2],A'*C',1);

lmiterm([1,2,2,-Y],-1*alpha,B,'s');
lmiterm([1,2,2,M],-1*alpha,C*B,'s');
lmiterm([1,2,2,N1],-1,1,'s');

lmiterm([1,2,3,Z],1,1);
lmiterm([1,2,3,N2],B'*C'*alpha,-1);

lmiterm([1,3,3,N2],-1,1,'s');

% N1 > 0
lmiterm([2,1,1,N1],-1,1);

% N2 > 0
lmiterm([3,1,1,N2],-1,1);

% P > 0 (where P is a block matrix of X,Y,Z)
lmiterm([4,1,1,X],-1,0.5,'s');
lmiterm([4,1,2,Y],-1,1);
lmiterm([4,2,2,Z],-1,0.5,'s');

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
dec     =  decnbr(LMISYS); % returns number of decision varibles
data.X  =  dec2mat(LMISYS,xfeas,X);
data.Y  =  dec2mat(LMISYS,xfeas,Y);
data.Z  =  dec2mat(LMISYS,xfeas,Z);
data.M  =  dec2mat(LMISYS,xfeas,M);
data.N1 =  dec2mat(LMISYS,xfeas,N1);
data.N2 =  dec2mat(LMISYS,xfeas,N2);

if D ~= zeros(m)
    disp('D not equal to zero. Park criterion may not be applied!'); 
    alpha   = nan;
    dec     =  nan;
    data.X  =  nan;
    data.Y  =  nan;
    data.Z  =  nan;
    data.M  =  nan;
    data.N1 =  nan;
    data.N2 =  nan;
end
end
