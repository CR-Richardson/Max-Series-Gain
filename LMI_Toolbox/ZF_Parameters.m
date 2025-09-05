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
% Script containing the Zames-Falb parameters for the max. series gain
% experiments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Syst{1}
W{1}  = diag([1 0.1 0.001]);
WB{1} = W{1};
WC{1} = W{1};

%% Syst{2}
W{2}  = eye(3);
WB{2} = W{2};
WC{2} = W{2};

%% Syst{3}
W{3}  = eye(4);
WB{3} = W{3};
WC{3} = W{3};

%% Syst{4}
W{4}  = eye(4);
WB{4} = W{4};
WC{4} = W{4};

%% Syst{5}
W{5}  = eye(4);
WB{5} = W{5};
WC{5} = W{5};

%% Syst{6}
W{6}  = eye(4);
WB{6} = W{6};
WC{6} = W{6};

%% Syst{7}
W{7}  = eye(4);
WB{7} = W{7};
WC{7} = W{7};

%% Syst{8}
W{8}  = eye(5);
WB{8} = W{8};
WC{8} = W{8};
