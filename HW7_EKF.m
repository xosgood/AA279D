% Chris Osgood
% 5/19/23
% AA279D HW7
% Continous Control
clc; clear; close all;
addpath(genpath("Functions/"));

%% constants
mu = 3.986e5;

%% initialize chief
% chief mean OE
a_c = 7000; % km
e_c = 0.001;
i_c = deg2rad(98);
RAAN_c = 0;
omega_c = deg2rad(90);
nu_c = deg2rad(-15);
M_c = TrueToMeanAnomaly(nu_c, e_c);
oe_c_mean = [a_c; e_c; i_c; RAAN_c; omega_c; nu_c];

n_c = sqrt(mu / a_c^3);

% chief osculating OE
oe_c_osc = mean2osc(oe_c_mean, 1);

% get chief initial state for ode113
[r_c_osc_0_ECI, v_c_osc_0_ECI] = OE2ECI(oe_c_osc);
x_c_osc_0 = [r_c_osc_0_ECI; v_c_osc_0_ECI];

%% initialize deputy
% deputy mean ROE
roe_d_mean = [0; 100; 0; 0; 0; 0] / a_c;

% deputy mean OE
oe_d = ROE2OE(oe_c_mean, roe_d_mean); % deputy mean oe
a_d = oe_d(1);
e_d = oe_d(2);
i_d = oe_d(3);
RAAN_d = oe_d(4);
omega_d = oe_d(5);
nu_d = oe_d(6);
M_d = TrueToMeanAnomaly(nu_d, e_d);

% deputy osculating OE
oe_d_mean = oe_d;
oe_d_osc = mean2osc(oe_d_mean, 1);
roe_d_osc = OE2ROE(oe_c_osc, oe_d_osc);

% get deputy initial state for ode113
[r_d_osc_0_ECI, v_d_osc_0_ECI] = OE2ECI(oe_d_osc);
x_d_osc_0 = [r_d_osc_0_ECI; v_d_osc_0_ECI];

%% EKF



