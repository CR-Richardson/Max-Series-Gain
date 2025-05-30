function [N] = var_count(n,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
% CR Richardson 
% ECS
% University of Southampton
% UK
%
% Date: 24/05/23
%
% Purpose: 
% Compute the number of decision variables for each criteria.
%
% Note: ZF is incorrect when checked against decnbr(LMISYS)!

% Parameters:
% n: n = dim(x)
% m: m = dim(u) = dim(y)
%
% Returns:
% N:   vector containing the number of decision variables for each criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Formulas used for counting the number of variables of different types of square nxn matrices
% Full matrix: n^2
% Diagonal matrix: n
% Symmetric matrix: 0.5*(n^2 - n) + n = 0.5*n*(n + 1)
%%
x   = 7; % number of criteria being compared
N   = zeros(1,x);

%% Circle
N(1)   = 0.5*n*(n+1) + m;

%% Circle-like
N(2)   = 0.5*n*(n+1) + 2*(m^2);

%% Popov
N(3)   = 0.5*n*(n+1) + 2*m;

%% Popov-like 1
N(4)   = 0.5*n*(n+1) + 4*(m^2);

%% Popov-like 2
N(5)   = 0.5*n*(n+1) + 2*(m^2) + m;

%% Park
N(6)   = 0.5*(n+m)*(n+m+1) + 10*m;

%% ZF (includes Circle and Popov multipliers - subtract 4*m to exclude them)
N(7)   = 5*0.5*n*(n+1) + 2*n + 3*m;

%%
% Display results

disp(' ');
title_str=['     ', '         Circle', '    Circle-Like', '          Popov', ...
            '   Popov-Like 1', '   Popov-Like 2', '           Park', ...
            '             ZF'];
fprintf('%5s %15s %15s %15s %15s %15s %15s %15s',title_str);
disp(' ');
fprintf('    N'); fprintf('%15d %14d %14d %14d %14d %14d %14d\n', N);

end