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
e_1 = 0.01;
i_1 = deg2rad(98.5);
RAAN_1 = deg2rad(0.5);
omega_1 = deg2rad(95);
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
r0_dot = v_RTN(1);

theta = nu_0 + omega_0;
theta_dot = norm(cross(r_ECI_0, v_ECI_0))/r0^2; % v0 / r0
state0 = [r_RTN; theta; r0; v_RTN; theta_dot; r0_dot];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, state_out] = ode113(@RelativeMotionDifEqRTN, tspan, state0, options);

% pull out RTN position and velocity
r_RTN_relsim = state_out(:,1:3)';
v_RTN_relsim = state_out(:,6:8)';

% plot RTN orbit
figure;
sgtitle("Relative Orbit (RTN) from Relative Motion Simulation")
PlotRTN(t, r_RTN_relsim, v_RTN_relsim);
% plot RTN orbit, non-dimensionalized by semi-major axis
figure;
sgtitle("Non-Dimensionalized Relative Orbit (RTN) from Relative Motion Simulation")
PlotRTN_NonDim(t, r_RTN_relsim, v_RTN_relsim, a_0);

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
[r_RTN_abssim, v_RTN_abssim] = ECI2RTN_Vectorized(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

% plot orbits in ECI
figure;
PlotEarth();
plot3(r_ECI_0(1,:), r_ECI_0(2,:), r_ECI_0(3,:), 'r.');
plot3(r_ECI_1(1,:), r_ECI_1(2,:), r_ECI_1(3,:), 'g.');
title("Absolute ECI Orbits");
legend("Chief", "Deputy");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

% plot RTN orbit
figure;
sgtitle("Relative Orbit (RTN) from Absolute Orbit Simulation")
PlotRTN(t, r_RTN_abssim, v_RTN_abssim);
% plot RTN orbit, non-dimensionalized by semi-major axis
figure;
sgtitle("Non-Dimensionalized Relative Orbit (RTN) from Absolute Orbit Simulation")
PlotRTN_NonDim(t, r_RTN_abssim, v_RTN_abssim, a_0);
% plot RTN in 3D
% figure;
% Plot3RTN(r_RTN_abssim, 0);
% title("Relative Motion in RTN from Absolute Orbit Simulation");
% figure;
% Plot3RTN(r_RTN_abssim, a_0);
% title("Non-Dimensionalized Relative Motion in RTN from Absolute Orbit Simulation");

%% d) error between relative sim and absolute sim
r_RTN_err = abs(r_RTN_relsim - r_RTN_abssim);
v_RTN_err = abs(v_RTN_relsim - v_RTN_abssim);
figure;
PlotRTN(t, r_RTN_err, v_RTN_err);
sgtitle("RTN Error Between Relative Motion and Absolute Orbit Simulations");

% Repeat parts b-c with a difference in semi-major axis
% deputy
da = 100;
a_1 = a_0 + da;

% repeating part c
[r_ECI_0, v_ECI_0] = OE2ECI(a_0, e_0, i_0, RAAN_0, omega_0, nu_0);
[r_ECI_1, v_ECI_1] = OE2ECI(a_1, e_1, i_1, RAAN_1, omega_1, nu_1);

[r_RTN_abssim, v_RTN_abssim] = ECI2RTN(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

r0 = norm(r_ECI_0);
v0 = norm(v_ECI_0);

theta = nu_0 + omega_0;
theta_dot = norm(cross(r_ECI_0, v_ECI_0))/r0^2;
state0 = [r_RTN_abssim; theta; r0; v_RTN_abssim; theta_dot; v0];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, state_out] = ode113(@RelativeMotionDifEqRTN, tspan, state0, options);

% pull out RTN position and velocity
r_RTN_relsim = state_out(:,1:3)';
v_RTN_relsim = state_out(:,6:8)';

% plot RTN orbit
figure;
sgtitle("Relative Orbit (RTN) from Relative Motion Simulation, with da")
PlotRTN(t, r_RTN_relsim, v_RTN_relsim);
% plot RTN orbit, non-dimensionalized by semi-major axis
figure;
sgtitle("Non-Dimensionalized Relative Orbit (RTN) from Relative Motion Simulation, with da")
PlotRTN_NonDim(t, r_RTN_relsim, v_RTN_relsim, a_0);

% repeating part c
abs_data_0 = SimulateOrbitFromOE(a_0, e_0, i_0, RAAN_0, omega_0, M_0, geod_station, t_0_MJD, t_f_MJD, n_iter);
abs_data_1 = SimulateOrbitFromOE(a_1, e_1, i_1, RAAN_1, omega_1, M_1, geod_station, t_0_MJD, t_f_MJD, n_iter);

r_ECI_0 = abs_data_0.r_ECI_vec;
v_ECI_0 = abs_data_0.v_ECI_vec;
r_ECI_1 = abs_data_1.r_ECI_vec;
v_ECI_1 = abs_data_1.v_ECI_vec;

% turn ECI to RTN for deputy state
[r_RTN_abssim, v_RTN_abssim] = ECI2RTN_Vectorized(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

% plot orbits in ECI
figure;
PlotEarth();
plot3(r_ECI_0(1,:), r_ECI_0(2,:), r_ECI_0(3,:), 'r.');
plot3(r_ECI_1(1,:), r_ECI_1(2,:), r_ECI_1(3,:), 'g.');
title("Absolute ECI Orbits");
legend("Chief", "Deputy");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

% plot RTN orbit
figure;
sgtitle("Relative Orbit (RTN) from Absolute Orbit Simulation, with da")
PlotRTN(t, r_RTN_abssim, v_RTN_abssim);
% plot RTN orbit, non-dimensionalized by semi-major axis
figure;
sgtitle("Non-Dimensionalized Relative Orbit (RTN) from Absolute Orbit Simulation, with da")
PlotRTN_NonDim(t, r_RTN_abssim, v_RTN_abssim, a_0);
% plot RTN in 3D
% figure;
% Plot3RTN(r_RTN_abssim, 0);
% title("Relative Motion in RTN from Absolute Orbit Simulation, with da");
% figure;
% Plot3RTN(r_RTN_abssim, a_0);
% title("Non-Dimensionalized Relative Motion in RTN from Absolute Orbit Simulation, with da");

% plot errors for this new sim with da
r_RTN_err = abs(r_RTN_relsim - r_RTN_abssim);
v_RTN_err = abs(v_RTN_relsim - v_RTN_abssim);
figure;
PlotRTN(t, r_RTN_err, v_RTN_err);
sgtitle("RTN Error Between Relative Motion and Absolute Orbit Simulations with da");


%% e+f) apply manuever

% use T of deputy
T = 2 * pi * sqrt(a_1^3 / mu);

n_orbits_pre_dv = 2 - nu_1/(2*pi); % after 2 orbits (minus the initial true anomaly), it is back at periapsis
n_iter = 30 * n_orbits_pre_dv;
t_f = n_orbits_pre_dv * T;
tspan = linspace(0, t_f, n_iter); % [0, t_f];

[t_pre_dv, state_out_pre_dv] = ode113(@RelativeMotionDifEqRTN, tspan, state0, options);

n_orbits_post_dv = 5;
n_iter = 30 * (n_orbits_post_dv - n_orbits_pre_dv);
t_f_2 = n_orbits_post_dv * T;
tspan = linspace(t_f, t_f_2, n_iter); % [0, t_f];

% compute dv
n = sqrt(mu / a_1^3);
r = norm([state_out_pre_dv(end,1) + state_out_pre_dv(end,5), state_out_pre_dv(end,2:3)]);
dv_T = -da * n * r / (2 * a_1 * sqrt(1 - e_1^2));

% apply dv
state0 = state_out_pre_dv(end,:);
state0(7) = state0(7) + dv_T;
[t_post_dv, state_out_post_dv] = ode113(@RelativeMotionDifEqRTN, tspan, state0, options);

r_RTN_pre_dv = state_out_pre_dv(:,1:3)';
v_RTN_pre_dv = state_out_pre_dv(:,6:8)';
r_RTN_post_dv = state_out_post_dv(:,1:3)';
v_RTN_post_dv = state_out_post_dv(:,6:8)';

t_tot_maneuver = [t_pre_dv; t_post_dv];
r_RTN_tot_maneuver = [r_RTN_pre_dv, r_RTN_post_dv];
v_RTN_tot_maneuver = [v_RTN_pre_dv, v_RTN_post_dv];

figure;
PlotRTN(t_tot_maneuver, r_RTN_tot_maneuver, v_RTN_tot_maneuver);
sgtitle("Relative Motion with Along-Track Maneuver to Cancel da");
