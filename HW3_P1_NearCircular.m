% Paxton Scott and Chris Osgood
% 4/21/23
% AA279D HW3 P1
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
oe_c = [a_0; e_0; i_0; RAAN_0; omega_0; nu_0];

% deputy
a_1 = a_0;
e_1 = 0.001;
i_1 = deg2rad(98.01);
RAAN_1 = deg2rad(0.05);
omega_1 = deg2rad(90.05);
nu_1 = deg2rad(0.01);
oe_d = [a_1; e_1; i_1; RAAN_1; omega_1; nu_1];

%% b) convert initial conditions to useful coord frames
[r_ECI_0, v_ECI_0] = OE2ECI(a_0, e_0, i_0, RAAN_0, omega_0, nu_0);
[r_ECI_1, v_ECI_1] = OE2ECI(a_1, e_1, i_1, RAAN_1, omega_1, nu_1);

[r_RTN_d, v_RTN_d] = ECI2RTN(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

% roe = [da, dlambda, dex, dey, dix, diy]^T
roe = OE2ROE(oe_c, oe_d);

v_ECI_c_expressedInRTN = rECI2RTN([r_ECI_0; v_ECI_0]) * v_ECI_0;

rho = norm(r_RTN_d);
r0 = norm(r_ECI_0);
r0_dot = v_ECI_c_expressedInRTN(1);

assert(rho / r0 <= 0.001);

theta = nu_0 + omega_0;
theta_dot = norm(cross(r_ECI_0, v_ECI_0))/r0^2; % v0 / r0
state0 = [r_RTN_d; theta; r0; v_RTN_d; theta_dot; r0_dot];

%% c)
