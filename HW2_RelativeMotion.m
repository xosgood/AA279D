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

%% sim parameters
n_orbits = 5;
n_iter = 30 * n_orbits;
mu = 3.986e5; % (km^3 / s^2)
T = 2 * pi * sqrt(a_0^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter); % [0, t_f];

%% b) relative orbit sim
[r_ECI_0, v_ECI_0] = OE2ECI(a_0, e_0, i_0, RAAN_0, omega_0, nu_0);
[r_ECI_1, v_ECI_1] = OE2ECI(a_1, e_1, i_1, RAAN_1, omega_1, nu_1);

[r_RTN, v_RTN] = ECI2RTN(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

r0 = norm(r_ECI_0);
v0 = norm(v_ECI_0);

theta = nu_0 - omega_0;
theta_dot = v0/r0;
state0 = [r_RTN; theta; r0; v_RTN; theta_dot; v0];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, state_out] = ode113(@RelativeMotionDifEqRTN, tspan, state0, options);

%% c) absolute orbit sim of chief and deputy
% set times
t_0_MJD = mjuliandate(2023,07,04,12,0,0); % 07/04/2023 converted to MJD
sec_per_day = 86400;
t_f_MJD = t_0_MJD + n_orbits * T / sec_per_day;
geod_station = [0; 0; 0]; % fake value since we don't care about this

M_0 = TrueToMeanAnomaly(nu_0, e_0);
M_1 = TrueToMeanAnomaly(nu_1, e_1);

abs_data_0 = SimulateOrbitFromOE(a_0, e_0, i_0, RAAN_0, omega_0, M_0, geod_station, t_0_MJD, t_f_MJD, n_iter);
abs_data_1 = SimulateOrbitFromOE(a_1, e_1, i_1, RAAN_1, omega_1, M_1, geod_station, t_0_MJD, t_f_MJD, n_iter);

r_ECI_0 = abs_data_0.r_ECI_vec;
v_ECI_0 = abs_data_0.v_ECI_vec;
r_ECI_1 = abs_data_1.r_ECI_vec;
v_ECI_1 = abs_data_1.v_ECI_vec;

% turn ECI to RTN for deputy state
[r_RTN_1, v_RTN_1] = ECI2RTN_Vectorized(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

% plot orbits in ECI
figure;
PlotEarth();
plot3(r_ECI_0(1,:), r_ECI_0(2,:), r_ECI_0(3,:), 'r.');
plot3(r_ECI_1(1,:), r_ECI_1(2,:), r_ECI_1(3,:), 'g.');
title("Absoulte ECI Orbits");
legend("Chief", "Deputy");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

% plot orbits in RTN
figure;
sgtitle("Relative Orbit (RTN) from Absolute Orbit Simulation")
PlotRTN(t, r_RTN_1, v_RTN_1);
% TODO




%% d) error between relative sim and absolute sim


