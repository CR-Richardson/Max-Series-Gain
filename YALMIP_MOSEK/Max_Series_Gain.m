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
% This script builds various linear systems and computes the maximum series
% gain (alpha) according to various criteria for which the Lurie system is
% stable when the repeated ReLU is placed in the feebdack path. 
% 
%
% Scripts
% Examples:    Contains example linear systems.
%
% Functions
% Circle:      Circle criterion - See Theorem 1 and Remark 3.
% Circle_Like: Circle-Like criterion - See Theorem 1.
% Popov:       Popov criterion - See Theorem 2 and Remark 4.
% Popov_Like1: Relaxed Popov-Like criterion - See Theorem 2 and Corollary 1.
% Popov_Like2: Relaxed Popov-Like criterion - See Theorem 2 and Corollary 2.
% Park:        Park's criterion - See Reference 10 Theorem 2.
% ZF:          Zames-Falb criterion - See Reference 30.
%
% Variables
% Total_Ex:    Total number of examples.
% Ex_array:    The set of example systems to compute alpha for.
% eps:         Loop termination accuracy - try 1e-6 for Examples 1-8 and 1e-2 for Examples 9-12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex   = 12;
Ex_array   = 1:8;
eps        = 1e-2;

%% Makes example systems accessible to script
Examples;

%% Nyquist gains of each example
Alpha_up = [100000.0, 89.9, 0.6983, 0.0020, 0.0869, 0.8202, ...
            0.2002, 2.0221, 2.1, 2.8, 2.3, 2.8];

Alpha_up = Alpha_up + 0.1; % ensures alpha_low != alpha_up for each example.

%% Arrays for storing the maximum series gain (alpha) and # of decision variables.

All_ex = 1:Total_Ex; % All example systems

Alpha_Circle      = zeros(All_ex);
Alpha_Circle_Like = zeros(All_ex);
Alpha_Popov       = zeros(All_ex);
Alpha_Popov_Like1 = zeros(All_ex);
Alpha_Popov_Like2 = zeros(All_ex);
Alpha_Park        = zeros(All_ex);
Alpha_ZF          = zeros(All_ex);

decs_Circle      = zeros(All_ex);
decs_Circle_Like = zeros(All_ex);
decs_Popov       = zeros(All_ex);
decs_Popov_Like1 = zeros(All_ex);
decs_Popov_Like2 = zeros(All_ex);
decs_Park        = zeros(All_ex);
decs_ZF          = zeros(All_ex);

%% Calculate maximum series gain using various criteria

for i=Ex_array
    disp(['Example ',num2str(i),' ']);
    
    disp('Circle calculations...');
    tic;
    [Alpha_Circle(i), data_Circle(i), decs_Circle(i)] = Circle(Syst{i}, eps, Alpha_up(i));
    toc;

    disp('Popov calculations...'); 
    tic;
    [Alpha_Popov(i), data_Popov(i), decs_Popov(i)] = Popov(Syst{i}, eps, Alpha_up(i), Alpha_Circle(i));
    toc;

    disp('Circle-like calculations...'); 
    tic;
    [Alpha_Circle_Like(i), data_Circle_Like(i), decs_Circle_Like(i)] = Circle_Like(Syst{i}, eps, Alpha_up(i), Alpha_Circle(i));
    toc;

    disp('Popov-like 1 calculations...'); 
    tic;
    [Alpha_Popov_Like1(i), data_Popov_Like1(i), decs_Popov_Like1(i)] = Popov_Like1(Syst{i}, eps, Alpha_up(i), Alpha_Circle(i));
    toc;

    disp('Popov-like 2 calculations...'); 
    tic;
    [Alpha_Popov_Like2(i), data_Popov_Like2(i), decs_Popov_Like2(i)] = Popov_Like2(Syst{i}, eps, Alpha_up(i), Alpha_Popov(i));
    toc;

    disp('Park calculations...'); 
    tic;
    [Alpha_Park(i), data_Park(i), decs_Park(i)] = Park(Syst{i}, eps, Alpha_up(i), Alpha_Popov(i));
    toc;

    disp('Zames-Falb calculations...');
    tic;
    [Alpha_ZF(i), data_ZF(i), decs_ZF(i)] = ZF(Syst{i}, eps, Alpha_up(i), Alpha_Popov(i));
    toc;

end

%% Display max series gain
disp(' ');
disp('Max. series gain');
title_str=['        Example', '         Circle', '    Circle-Like', '          Popov', ...
           '   Popov-Like 1', '   Popov-Like 2', '           Park', '             ZF'];
mat_data =[Ex_array' Alpha_Circle(Ex_array)' Alpha_Circle_Like(Ex_array)' Alpha_Popov(Ex_array)' ...
           Alpha_Popov_Like1(Ex_array)' Alpha_Popov_Like2(Ex_array)' Alpha_Park(Ex_array)' Alpha_ZF(Ex_array)'];
fprintf('%15s %15s %15s %15s %15s %15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f\n',mat_data');

%% Display # of decision variables
disp(' ');
disp('# of decision variables');
title_str=['        Example', '         Circle', '    Circle-Like', '          Popov', ...
           '   Popov-Like 1', '   Popov-Like 2', '           Park', '             ZF'];
mat_data =[Ex_array' decs_Circle(Ex_array)' decs_Circle_Like(Ex_array)' decs_Popov(Ex_array)' ...
           decs_Popov_Like1(Ex_array)' decs_Popov_Like2(Ex_array)' decs_Park(Ex_array)' decs_ZF(Ex_array)'];
fprintf('%15s %15s %15s %15s %15s %15s %15s %15s\n',title_str);
disp(' ');
fprintf('%15d %14.0f %14.0f %14.0f %14.0f %14.0f %14.0f %14.0f\n',mat_data');
