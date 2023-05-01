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
[a_d; e_d; i_d; RAAN_d; omega_d; nu_d]  = oe_d;


%% 3) full non-linear simulation

% sim parameters
n_orbits = 15;
n_steps_per_orbit = 30;
n_iter = n_steps_per_orbit * n_orbits;
mu = 3.986e5; % (km^3 / s^2)
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter); % [0, t_f];

% set times
t_0_MJD = mjuliandate(2023,07,04,12,0,0); % 07/04/2023 converted to MJD
sec_per_day = 86400;
t_f_MJD = t_0_MJD + n_orbits * T / sec_per_day;
geod_station = [0; 0; 0]; % fake value since we don't care about this

M_c = TrueToMeanAnomaly(nu_c, e_c);
M_d = TrueToMeanAnomaly(nu_d, e_d);

abs_data_c = SimulateOrbitFromOE(a_c, e_c, i_c, RAAN_c, omega_c, M_c, geod_station, t_0_MJD, t_f_MJD, n_iter);
abs_data_d = SimulateOrbitFromOE(a_d, e_d, i_d, RAAN_d, omega_d, M_d, geod_station, t_0_MJD, t_f_MJD, n_iter);

r_ECI_c = abs_data_c.r_ECI_vec;
v_ECI_c = abs_data_c.v_ECI_vec;
r_ECI_d = abs_data_d.r_ECI_vec;
v_ECI_d = abs_data_d.v_ECI_vec;

% turn ECI to RTN for deputy state
[r_RTN_abssim, v_RTN_abssim] = ECI2RTN_Vectorized(r_ECI_c, v_ECI_c, r_ECI_d, v_ECI_d);

% plot orbits in ECI
figure;
PlotEarth();
plot3(r_ECI_c(1,:), r_ECI_c(2,:), r_ECI_c(3,:), 'r.');
plot3(r_ECI_d(1,:), r_ECI_d(2,:), r_ECI_d(3,:), 'g.');
title("Absolute ECI Orbits");
legend("Chief", "Deputy");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

