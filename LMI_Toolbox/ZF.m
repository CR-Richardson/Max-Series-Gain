function [alpha,ZFmult,data,dec]=ZF(syst,WB,WC,Pflag,alpha_low)

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
% Compute the maximum series gain (alpha) when using the ZF criterion as 
% defined in Reference 30. 
%
% Parameters:
% syst:      Structure containing the system matrices of an example.
% WB=WC:     User defined diagonal matrices
% Pflag:     1 = Popov multiplier used; 0 Popov multiplier not used
% alpha_low: lower bound of series gain (e.g. from Circle Crit.)

% Returns:
% alpha:  Maximum series gain (float)
% ZFmult: Associated Zames-Falb multiplier which proves stability
% data:   Structure containing solutions of the LMI parametrised by alpha
% dec:    # number of decision variables
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
    Gm = 1000;
end

% Determine initial upper/lower bound and initial test value
alpha_up = Gm*0.999;
if ((alpha_up - alpha_low)/alpha_up) < 0.0001
   alpha  = -alpha_low;
   ZFmult = tf(0,1);
   Vm     = 0;
   etam   = 0;
else
   alpha  = alpha_up;
   ZFmult = [];
end
 
eps = 1e-6;

%%
% Determine alpha by repeatedly solving LMI until the largest alpha is 
% found where LMI is feasible

while ((alpha_up - alpha_low)/alpha_up) > eps
  
setlmis([]);

P11 = lmivar(1,[n,1]);
S11 = lmivar(1,[n,1]);
N   = lmivar(1,[n,1]);
H0  = lmivar(1,kron([1,0],ones(m,1)));
if Pflag == 1
   etaP = lmivar(1,kron([1,0],ones(m,1))); % Popov mutiplier - diagonal  
   VP   = lmivar(1,kron([1,0],ones(m,1))); % Circle multiplier - diagonal
end 
Acb = lmivar(1,[n,1]);
Aab = lmivar(1,[n,1]);
[~,~,sBc] = lmivar(2,[n,1]);
[~,~,sBa] = lmivar(2,[n,1]);
Bcb = lmivar(3,kron(sBc,ones(1,m))); % columns identical
Bab = lmivar(3,kron(sBa,ones(1,m)));
        
% LMI
lmiterm([1,1,1,S11],1,A,'s');
lmiterm([1,1,1,Bab],-1,WC*C,'s');

lmiterm([1,1,2,S11],1,A);
lmiterm([1,1,2,P11],A',1);
lmiterm([1,1,2,-Acb],1,1);
lmiterm([1,1,2,-Aab],1,1);
lmiterm([1,1,2,-Bcb],C'*WB,1);
lmiterm([1,1,2,Bab],-1,WC*C);

lmiterm([1,1,3,N],A',1);
lmiterm([1,1,3,-Bab],C'*WC,1);
lmiterm([1,1,3,-Aab],1,1);

lmiterm([1,1,4,S11],1*alpha,B);
lmiterm([1,1,4,Bcb],1,WC);
lmiterm([1,1,4,Bab],1,WC);
lmiterm([1,1,4,-H0],C',1);
if Pflag == 1
   lmiterm([1,1,4,etaP],A'*C',1);   % Popov
   lmiterm([1,1,4,-VP],C',1);       % Circle
end

lmiterm([1,2,2,P11],1,A,'s');
lmiterm([1,2,2,Bcb],1,WB*C,'s');
 
lmiterm([1,2,3,N],A',1);
lmiterm([1,2,3,-Bab],C'*WC,1);
lmiterm([1,2,3,Aab],-1,1);

lmiterm([1,2,4,P11],1*alpha,B);
lmiterm([1,2,4,Bcb],-1,WB);
lmiterm([1,2,4,Bab],-1,WB);
lmiterm([1,2,4,-H0],C',1);
if Pflag == 1
   lmiterm([1,2,4,etaP],A'*C',1); % Popov
   lmiterm([1,2,4,-VP],C',1);     % Circle
end

lmiterm([1,3,3,Aab],-1,1,'s');

lmiterm([1,3,4,N],1*alpha,B);
lmiterm([1,3,4,Bab],-1,WB+WC);

lmiterm([1,4,4,H0],-1,1,'s');
if Pflag == 1
   lmiterm([1,4,4,etaP],1*alpha,C*B,'s'); % Popov
   lmiterm([1,4,4,VP],1,-1,'s');          % Circle 
end
        
% L1 bound LMI
E  = eye(m);
e1 = E(:,1);
lmiterm([2,1,1,0],-1);
lmiterm([2,1,2,-Bcb],-e1',1);
lmiterm([2,1,3,-Bab],-e1',1);
lmiterm([2,2,2,Acb],-1,0.5,'s');
lmiterm([2,2,3,0],zeros(n,n));
lmiterm([2,3,3,Aab],-1,0.5,'s');
        
% L1 bound LMIs - Assumes H0 is diagonal
for i = 1:m
    lmiterm([2+i,1,1,0],trace(WB)*E(i,:)*WC*E(:,i));
    lmiterm([2+i,1,1,H0],-E(i,:),E(:,i));
end
lmiterm([2+m+1,1,1,S11],-1,1,'s'); % Positive definite conditions
lmiterm([2+m+2,1,1,N],-1,1,'s');
lmiterm([2+m+3,1,1,P11],-1,1,'s');
lmiterm([2+m+3,1,1,S11],1,1,'s');
lmiterm([2+m+3,1,1,N],1,1,'s');
if Pflag == 1
   lmiterm([2+m+4,1,1,VP],1,-1,'s');
end 

LMISYS = getlmis;
[tmin,xfeas] = feasp(LMISYS,[1e-20 5000 -0.1 1000 1]);

% Update alpha upper/lower bound plus new test value 
if tmin < 0  
   alpha_low = alpha; % if LMIs are feasible
else 
   alpha_up = alpha; % if LMIs are infeasible
end

alpha = (alpha_up + alpha_low)/2;
end
       
%% Return solutions
dec      =  decnbr(LMISYS); % returns number of decision varibles
data.P11 =  dec2mat(LMISYS,xfeas,P11);
data.S11 =  dec2mat(LMISYS,xfeas,S11);
data.N   =  dec2mat(LMISYS,xfeas,N);
data.H0  =  dec2mat(LMISYS,xfeas,H0);
data.Acb =  dec2mat(LMISYS,xfeas,Acb);
data.Aab =  dec2mat(LMISYS,xfeas,Aab);
data.Bcb =  dec2mat(LMISYS,xfeas,Bcb);
data.Bab =  dec2mat(LMISYS,xfeas,Bab);
if Pflag == 1
   data.etaP =  dec2mat(LMISYS,xfeas,etaP); % Popov multiplier
   data.VP   =  dec2mat(LMISYS,xfeas,VP);   % Circle multiplier
end

if D ~= zeros(m)
    disp('D not equal to zero. ZF criterion may not be applied!'); 
    alpha     = nan;
    dec       = nan;
    data.P11  = nan;
    data.S11  = nan;
    data.N    = nan;
    data.H0   = nan;
    data.Acb  = nan;
    data.Aab  = nan;
    data.Bcb  = nan;
    data.Bab  = nan;
    data.etaP = nan;
    data.VP   = nan;
end

end
