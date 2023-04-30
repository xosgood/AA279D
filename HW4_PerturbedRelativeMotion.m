% Paxton Scott and Chris Osgood
% 4/30/23
% AA279D HW4
clc; clear; close all;
addpath(genpath("Functions/"));

%% 1) initial chief orbit and setup
% constants
mu = 3.986e5;

% chief OE
a_c = 7000; % km
e_c = 0.001;
i_c = deg2rad(98);
RAAN_c = 0;
omega_c = deg2rad(90);
nu_c = 0;
oe_c = [a_c; e_c; i_c; RAAN_c; omega_c; nu_c];

%% 2) initial deputy roe
% deputy ROE (quasi non-singular) [da, dlambda, dex, dey, dix, diy]^T
roe_d = [0; 0.100; 0.050; 0.100; 0.030; 0.200] / a_c ; % [m]

oe_d = ROE2OE(oe_c, roe_d);

%% 3) full non-linear simulation

