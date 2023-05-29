% Chris Osgood & Paxton Scott
% 5/29/23
% AA279D HW7
% Implementing the Unscented Khalman filter
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

%% simulation parameters
n_orbits = 300;
n_steps_per_orbit = 100;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
dt = tspan(2) - tspan(1);
orbit_span = (1:n_iter)/n_steps_per_orbit;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

%% chief absolute orbit full nonlinear propagator 
[~, x_c_ECI] = ode113(@AbsoluteOrbitWithJ2DiffEq, tspan, x_c_osc_0, options);
x_c_ECI = x_c_ECI';
r_c_ECI = x_c_ECI(1:3,:);
v_c_ECI = x_c_ECI(4:6,:);

%% deputy absolute orbit full nonlinear propagator (for ground truth)
[~, x_d_ECI] = ode113(@AbsoluteOrbitWithJ2DiffEq, tspan, x_d_osc_0, options);
x_d_ECI = x_d_ECI';
r_d_ECI = x_d_ECI(1:3,:);
v_d_ECI = x_d_ECI(4:6,:);
