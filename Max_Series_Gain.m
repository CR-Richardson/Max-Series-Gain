
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
% This script builds various linear systems and computes the maximum series
% gain (alpha) according to various criteria for which the Lurie system is
% stable when the repeated ReLU is placed in the feebdack path. 
% 
% Note: Only Circle and Circle-like criteria can handle D~=0. All other 
% criteria will return nan. 
%
% Scripts
% Examples:     Contains example linear systems
% ZF_Paramters: Contains Zames-Falb parameters for each example
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
% Total_Ex: Total number of examples
% Ex_array: The set of example systems to compute alpha for
% eta:      Parameter of Popov_Like1
% Pflag:    Use Popov multiplier with ZF multiplier 1/0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script variables
clear all; close all;
Total_Ex   = 8;
Ex_array   = 1:8;
eta        = 10000;
Pflag      = 1;

%% Makes example systems and ZF parameters accessible to script
Examples;
ZF_Parameters;

%% Calculate maximum series gain using various criteria

All_ex = 1:Total_Ex; % All example systems

% Arrays for storing the maximum series gain (alpha) for each example.
Alpha_Circle      = zeros(All_ex);
Alpha_Circle_Like = zeros(All_ex);
Alpha_Popov       = zeros(All_ex);
Alpha_Popov_Like1 = zeros(All_ex);
Alpha_Popov_Like2 = zeros(All_ex);
Alpha_Park        = zeros(All_ex);
Alpha_ZF          = zeros(All_ex);

% Arrays for storing the # of decision variables for each example.
decs_Circle      = zeros(All_ex);
decs_Circle_Like = zeros(All_ex);
decs_Popov       = zeros(All_ex);
decs_Popov_Like1 = zeros(All_ex);
decs_Popov_Like2 = zeros(All_ex);
decs_Park        = zeros(All_ex);
decs_ZF          = zeros(All_ex);

for i=Ex_array
    disp(['Example ',num2str(i),' ']);
    
    disp('Circle calculations...'); 
    [Alpha_Circle(i), data1(i), decs_Circle(i)] = Circle(Syst{i});
    
    disp('Circle-like calculations...'); 
    [Alpha_Circle_Like(i), data2(i), decs_Circle_Like(i)] = Circle_Like(Syst{i});

    disp('Popov calculations...'); 
    [Alpha_Popov(i), data3(i), decs_Popov(i)] = Popov(Syst{i});
    
    disp('Popov-like 1 calculations...'); 
    [Alpha_Popov_Like1(i), data4(i), decs_Popov_Like1(i)] = Popov_Like1(Syst{i},eta);

    disp('Popov-like 2 calculations...'); 
    [Alpha_Popov_Like2(i), data5(i), decs_Popov_Like2(i)] = Popov_Like2(Syst{i});

    disp('Park calculations...'); 
    [Alpha_Park(i), data6(i), decs_Park(i)] = Park(Syst{i});

    disp('Zames-Falb calculations...');
    [Alpha_ZF(i),~, data7(i), decs_ZF(i)] = ZF(Syst{i},WB{i},WC{i},Pflag,0.9*Alpha_Circle(i));
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
