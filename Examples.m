
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
% This script builds various linear systems and computes the maximum series
% gain (alpha) according to various criteria for which the Lurie system is
% stable when the repeated ReLU is placed in the feebdack path. 
% 
% Note: As stated in the paper - the D matrix is assumed to be zero.
%
% Functions:
% Circle: Circle criterion - See Theorem 1 and Remark 4.
% Circle_Like: Circle-Like criterion - See Theorem 1.
% Popov: Popov criterion - See Theorem 2 and Remark 6.
% Popov_Like1: Relaxed Popov-Like criterion - See Theorem 2 and Corollary 1.
% Popov_Like2: Relaxed Popov-Like criterion - See Theorem 2 and Corollary 2.
% Park: Park's criterion - See Reference 13 Theorem 2.
% ZF: Zames-Falb criterion - See References 24 and 25.
%
% Variables:
% Ex_array: The set of example systems to compute alpha for.
% eta: Parameter of Popov_Like1
% Pflag: Use Popov multiplier with ZF multiplier 1/0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all; close all;
Ex_array = 1:8;
eta      = 10000;
Pflag    = 1;

%%
% Example systems

%--------------------------------------------------------------------------
% Park - IEEETAC 2002 - Example 3
A        = -diag([1 4 6 2 9 8 3 10 12]);
B        = -[1 0 0 0 1 0 0 1 0;
             0 1 0 1 0 0 0 0 1;
             0 0 1 0 0 1 1 0 0]';
C        =  [1 1 1 0 0 0 0 0 0;
             0 0 0 1 1 1 0 0 0;
             0 0 0 0 0 0 1 1 1];
D        =  zeros(3,3);
Syst{1}  = ss(A,B,C,D);
clear A B C D;
%--------------------------------------------------------------------------
% Fetzer and Scherer - IJRNC 2017 - Example 5.6
delt    = 0.00001;
P11     = tf([0 -1 0 +1],conv([1 delt 1],[1 10]));
P12     = tf([0 1],[1 1]);
P22     = tf([0 0 0 0 -40],conv(conv([1 delt],[1 1]),[1 0.8 16]));
Pbig    = [P11 P12;P12 P22];
Pbig    = ss(Pbig);
Pbig    = minreal(Pbig);
Syst{2} = Pbig;
clear P11 P12 P22 Pbig;
%--------------------------------------------------------------------------
% Fetzer and Scherer - IJRNC 2017 - Example 5.7 
delt    = 0.00001;
P22     = tf([0 0 0 0 -40],conv(conv([1 delt],[1 1]),[1 0.8 16]));
Pbig    = P22;
Pbig    = ss(Pbig);
Syst{3} = Pbig;
clear P22 Pbig;
%--------------------------------------------------------------------------
% Fetzer and Scherer - IJRNC 2017 - Example 5.6
delt    = 0.00001;
P11     = tf([0 -1 0 -1],conv([1 delt 1],[1 10]));
P12     = tf([0 1],[1 1]);
P22     = tf([0 0 0 0 -40],conv(conv([1 delt],[1 1]),[1 0.8 16]));
Pbig    = [P11 P12;P12 P22];
Pbig    = ss(Pbig);
Syst{4} = Pbig;
clear P11 P12 P22 Pbig;
%--------------------------------------------------------------------------
% Drummond, Guiver, and Turner - IEEETAC 2022 - Example 2
A = diag([-1 -10 -30 -60 -100])+0.1*ones(5,5);
B = eye(5);
C = 1*ones(5,5)-eye(5);
D = zeros(5,5);
Syst{5} = ss(A,B,C,D);
clear A B C D;
%--------------------------------------------------------------------------
% Turner, Kerr, and Sofrony - IJRNC 2015 - Example 17
A = [-5.4966   -1.7013   -1.7006   -1.0871   -0.3885   -0.2560
      8.0000         0         0         0         0         0
           0    2.0000         0         0         0         0
           0         0    2.0000         0         0         0
           0         0         0    2.0000         0         0
           0         0         0         0    0.5000         0];
B = [0.5751    0.3840    0.0158    0.6315
     0.4514    0.6831    0.0164    0.7176
     0.0439    0.0928    0.1901    0.6927
     0.0272    0.0353    0.5869    0.0841
     0.3127    0.6124    0.0576    0.4544
     0.0129    0.6085    0.3676    0.4418];
C = [0.3533    0.7275    0.4508    0.2548    0.9084    0.0784
     0.1536    0.4784    0.7159    0.8656    0.2319    0.6408
     0.6756    0.5548    0.8928    0.2324    0.2393    0.1909
     0.6992    0.1210    0.2731    0.8049    0.0498    0.8439];
D = zeros(4,4);
Syst{6}=ss(A,B,C,D);
clear A B C D;
%--------------------------------------------------------------------------
% Turner, Kerr, and Sofrony - IJRNC 2015 - Example 18
A = [-5.4966   -1.7013   -1.7006   -1.0871   -0.3885   -0.2560
      8.0000         0         0         0         0         0
           0    2.0000         0         0         0         0
           0         0    2.0000         0         0         0
           0         0         0    2.0000         0         0
           0         0         0         0    0.5000         0];
B = [-0.1747    0.1242   -0.1618   -0.4348
      0.1069   -0.4609    0.3498    0.2012
     -0.2069   -0.2954   -0.1913    0.5389
     -0.3018    0.4648    0.6416    0.3872
     -0.3921   -0.0244    0.4655   -0.5162
     -0.8676   -0.6227    0.1863    0.2499];
C = [-0.0874   -0.4380    0.2355   -0.1751   -0.3313   -0.0954
     -0.5379    0.0904   -0.0862    0.0422    0.0770   -0.2257
     -0.3988    0.1378   -0.1092    0.2767    0.1642   -0.2763
      0.0327   -0.4871   -0.4322    0.2185   -0.8372    0.2476];
D = zeros(4,4);
Syst{7}=ss(A,B,C,D);
clear A B C D;
%--------------------------------------------------------------------------
% Turner, Kerr, and Sofrony - IJRNC 2015 - Example 22
A = [-5.8097   -1.7258   -1.3496   -0.8882   -0.3991   -0.2843   -0.3009   -0.3524
      8.0000         0         0         0         0         0         0         0
           0    2.0000         0         0         0         0         0         0
           0         0    2.0000         0         0         0         0         0
           0         0         0    2.0000         0         0         0         0
           0         0         0         0    1.0000         0         0         0
           0         0         0         0         0    0.5000         0         0
           0         0         0         0         0         0    0.2500         0];
B = [-0.0697    0.1982   -0.4746   -0.3289
     -0.4567    0.3788   -0.2270    0.5888
      0.3737   -0.6601    0.6961   -0.3923
     -0.7929   -0.6439   -0.0461   -0.2640
      0.0148    0.4926    0.5402    0.7878
     -0.5594   -0.0041    0.2616    0.4186
      0.0385   -0.5487   -0.5681    0.5668
     -0.1095    0.0133    0.6951    0.7202];
C = 0.01*[0.5336   -0.2416    0.4056    0.0167    0.0115   -0.6213    0.4604    0.1259
          0.1910    0.1106    0.3421   -0.1956    0.2760    0.9206    0.4323   -0.7717
          0.2772    0.2116   -0.2094    0.4855    0.0615   -0.3314   -0.7378   -0.6914
         -0.6986   -0.0247   -0.0576    0.2415    0.3123    0.2021   -0.0033   -0.1274];
D = zeros(4,4);
Syst{8}=ss(A,B,C,D);
clear A B C D;
%%
% Defining parameters needed for Zames-Falb method

%--------------------------------------------------------------------------
% Syst{1}
W{1}  = diag([1 0.1 0.001]);
WB{1} = W{1};
WC{1} = W{1};
%--------------------------------------------------------------------------
% Syst{2}
W{2} = diag([1 0.001]);
WB{2} = W{2};
WC{2} = W{2};
%--------------------------------------------------------------------------
% Syst{3} (SISO)
W{3} = 1;
WB{3} = W{3};
WC{3} = W{3};
%--------------------------------------------------------------------------
% Syst{4}
W{4} = diag([0.0001  1]);
WB{4} = W{4};
WC{4} = W{4};
%--------------------------------------------------------------------------
% Syst{5}
W{5} = eye(5);
WB{5} = W{5};
WC{5} = W{5};
%--------------------------------------------------------------------------
% Syst{6}
W{6} = eye(4);
WB{6} = W{6};
WC{6} = W{6};
%--------------------------------------------------------------------------
% Syst{7}
W{7} = eye(4);
WB{7} = W{7};
WC{7} = W{7};
%--------------------------------------------------------------------------
% Syst{8}
W{8} = eye(4);
WB{8} = W{8};
WC{8} = W{8};
%--------------------------------------------------------------------------
%% 
% Calculate maximum series gain using various criteria

All_ex = 1:8; % All example systems

% Arrays for storing the maximum series gain (alpha) for each example.
Alpha_Circle      = zeros(All_ex);
Alpha_Circle_Like = zeros(All_ex);
Alpha_Popov       = zeros(All_ex);
Alpha_Popov_Like1 = zeros(All_ex);
Alpha_Popov_Like2 = zeros(All_ex);
Alpha_Park        = zeros(All_ex);
Alpha_ZF          = zeros(All_ex);

for i=Ex_array
    disp(['Example ',num2str(i),' ']);
    
    disp('Circle calculations...'); 
    [Alpha_Circle(i), data1(i)] = Circle(Syst{i});
    
    disp('Circle-like calculations...'); 
    [Alpha_Circle_Like(i), data2(i)] = Circle_Like(Syst{i});

    disp('Popov calculations...'); 
    [Alpha_Popov(i), data3(i)] = Popov(Syst{i});
    
    disp('Popov-like 1 calculations...'); 
    [Alpha_Popov_Like1(i), data4(i)] = Popov_Like1(Syst{i},eta);

    disp('Popov-like 2 calculations...'); 
    [Alpha_Popov_Like2(i), data5(i)] = Popov_Like2(Syst{i});

    disp('Park calculations...'); 
    [Alpha_Park(i), data6(i)] = Park(Syst{i});

    disp('Zames-Falb calculations...');
    Alpha_ZF(i) = ZF(Syst{i},WB{i},WC{i},Pflag,0.9*Alpha_Circle(i));
end

%%
% Display results

disp(' ');
title_str=['        Example', '         Circle', '    Circle-Like', '          Popov', ...
           '   Popov-Like 1', '   Popov-Like 2', '           Park', '             ZF'];
mat_data =[Ex_array' Alpha_Circle(Ex_array)' Alpha_Circle_Like(Ex_array)' Alpha_Popov(Ex_array)' ...
           Alpha_Popov_Like1(Ex_array)' Alpha_Popov_Like2(Ex_array)' Alpha_Park(Ex_array)' Alpha_ZF(Ex_array)'];
fprintf('%15s %15s %15s %15s %15s %15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f\n',mat_data');
