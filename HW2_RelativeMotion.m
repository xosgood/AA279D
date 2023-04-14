% Paxton Scott and Chris Osgood
% 4/14/23
% AA279D HW1
clc; clear; close all;
addpath(genpath("Functions/"));

%% a) initial orbits
% chief
a_0 = 7000; % km
e_0 = 0.001;
i_0 = deg2rad(98);
RAAN_0 = 0;
omega_0 = deg2rad(90);
nu_0 = 0;

% deputy
a_1 = a_0;
e_1 = 0.001;
i_1 = deg2rad(98.25);
RAAN_1 = deg2rad(0.25);
omega_1 = deg2rad(90.25);
nu_1 = deg2rad(0.5);

%% b) relative orbit sim


%% c) absolute orbit sim of chief and deputy



%% d) error between relative sim and absolute sim


