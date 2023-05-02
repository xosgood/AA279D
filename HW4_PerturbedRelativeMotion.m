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
a_d = oe_d(1);
e_d = oe_d(2);
i_d = oe_d(3);
RAAN_d = oe_d(4);
omega_d = oe_d(5);
nu_d = oe_d(6);


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
figure(1);
PlotEarth();
plot3(r_ECI_c(1,:), r_ECI_c(2,:), r_ECI_c(3,:), 'r');
plot3(r_ECI_d(1,:), r_ECI_d(2,:), r_ECI_d(3,:), 'g');
title("Absolute ECI Orbits, without J2");
legend("Earth", "Chief", "Deputy", "Location", "best");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

abs_data_c_j2 = SimulateOrbitFromOE_WithJ2(a_c, e_c, i_c, RAAN_c, omega_c, M_c, geod_station, t_0_MJD, t_f_MJD, n_iter);
abs_data_d_j2 = SimulateOrbitFromOE_WithJ2(a_d, e_d, i_d, RAAN_d, omega_d, M_d, geod_station, t_0_MJD, t_f_MJD, n_iter);

r_ECI_c_j2 = abs_data_c_j2.r_ECI_vec;
v_ECI_c_j2 = abs_data_c_j2.v_ECI_vec;
r_ECI_d_j2 = abs_data_d_j2.r_ECI_vec;
v_ECI_d_j2 = abs_data_d_j2.v_ECI_vec;

% turn ECI to RTN for deputy state.
[r_RTN_abssim_j2, v_RTN_abssim_j2] = ECI2RTN_Vectorized(r_ECI_c_j2, v_ECI_c_j2, r_ECI_d_j2, v_ECI_d_j2);

% plot orbits in ECI with the j2 effect. 
figure(2);
PlotEarth();
plot3(r_ECI_c_j2(1,:), r_ECI_c_j2(2,:), r_ECI_c_j2(3,:), 'r');
plot3(r_ECI_d_j2(1,:), r_ECI_d_j2(2,:), r_ECI_d_j2(3,:), 'g');
title("Absolute ECI Orbits, with J2");
legend("Earth", "Chief", "Deputy", "Location", "best");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

% convert r, v in time series into orbital elements.
oe_c_series = zeros(6, n_iter);
oe_d_series = zeros(6, n_iter);
oe_c_j2_series = zeros(6, n_iter);
oe_d_j2_series = zeros(6, n_iter);

QNS_oe_c_series = zeros(6, n_iter);
QNS_oe_d_series = zeros(6, n_iter);
QNS_oe_c_j2_series = zeros(6, n_iter);
QNS_oe_d_j2_series = zeros(6, n_iter);

QNS_roe_series = zeros(6, n_iter);
QNS_roe_j2_series = zeros(6, n_iter);

oe_c_mean_series = zeros(6, n_iter);
oe_d_mean_series = zeros(6, n_iter);
oe_c_mean_j2_series = zeros(6, n_iter);
oe_d_mean_j2_series = zeros(6, n_iter);

QNS_oe_c_mean_series = zeros(6, n_iter);
QNS_oe_d_mean_series = zeros(6, n_iter);
QNS_oe_c_mean_j2_series = zeros(6, n_iter);
QNS_oe_d_mean_j2_series = zeros(6, n_iter);

QNS_roe_mean_series = zeros(6, n_iter);
QNS_roe_mean_j2_series = zeros(6, n_iter);

for iter = 1:size(r_ECI_c,2)
    oe_c_series(:, iter) = ECI2OE(r_ECI_c(:, iter), v_ECI_c(:, iter))';
    oe_d_series(:, iter) = ECI2OE(r_ECI_d(:, iter), v_ECI_d(:, iter))';
    oe_c_j2_series(:, iter) = ECI2OE(r_ECI_c_j2(:, iter), v_ECI_c_j2(:, iter))';
    oe_d_j2_series(:, iter) = ECI2OE(r_ECI_d_j2(:, iter), v_ECI_d_j2(:, iter))';

    % a. Osculating quasi-non-singular orbital elements
    QNS_oe_c_series(:, iter) = OE2QNS_OE(oe_c_series(:, iter));
    QNS_oe_d_series(:, iter) = OE2QNS_OE(oe_d_series(:, iter));
    QNS_oe_c_j2_series(:, iter) = OE2QNS_OE(oe_c_j2_series(:, iter));
    QNS_oe_d_j2_series(:, iter) = OE2QNS_OE(oe_d_j2_series(:, iter));

    % b. Osculating relative quasi-non-singular orbital elements
    QNS_roe_series(:, iter) = OE2ROE(oe_c_series(:, iter), oe_d_series(:, iter));
    QNS_roe_j2_series(:, iter) = OE2ROE(oe_c_j2_series(:, iter), oe_d_j2_series(:, iter));
    
    % c. Mean quasi-non-singular orbital elements
    % convert [a [km]; e; i; RAAN; omega; nu] to [a [km]; e; i; RAAN; w; M]
    osc_elem_c = oe_c_series(:, iter);
    osc_elem_c(1) = 1000 * osc_elem_c(1);
    osc_elem_c(6) = TrueToMeanAnomaly(osc_elem_c(6), osc_elem_c(2));
    osc_elem_d = oe_d_series(:, iter);
    osc_elem_d(1) = 1000 * osc_elem_d(1);
    osc_elem_d(6) = TrueToMeanAnomaly(osc_elem_d(6), osc_elem_d(2));
    osc_elem_c_j2 = oe_c_j2_series(:, iter);
    osc_elem_c_j2(1) = 1000 * osc_elem_c_j2(1);
    osc_elem_c_j2(6) = TrueToMeanAnomaly(osc_elem_c_j2(6), osc_elem_c_j2(2));
    osc_elem_d_j2 = oe_d_j2_series(:, iter);
    osc_elem_d_j2(1) = 1000 * osc_elem_d_j2(1);
    osc_elem_d_j2(6) = TrueToMeanAnomaly(osc_elem_d_j2(6), osc_elem_d_j2(2));
    
    oe_c_mean_series(:, iter) = osc2mean(osc_elem_c, 0);
    oe_d_mean_series(:, iter) = osc2mean(osc_elem_d, 0);
    oe_c_mean_j2_series(:, iter) = osc2mean(osc_elem_c_j2, 1);
    oe_d_mean_j2_series(:, iter) = osc2mean(osc_elem_d_j2, 1);
    
    oe_c_mean_series(1, iter) = oe_c_mean_series(1, iter) / 1000;
    oe_d_mean_series(1, iter) = oe_d_mean_series(1, iter) / 1000;
    oe_c_mean_j2_series(1, iter) = oe_c_mean_j2_series(1, iter) / 1000;
    oe_d_mean_j2_series(1, iter) = oe_d_mean_j2_series(1, iter) / 1000;
    
    oe_c_mean_series(6, iter) = MeanToTrueAnomaly(oe_c_mean_series(6, iter), oe_c_mean_series(2, iter));
    oe_d_mean_series(6, iter) = MeanToTrueAnomaly(oe_d_mean_series(6, iter), oe_d_mean_series(2, iter));
    oe_c_mean_j2_series(6, iter) = MeanToTrueAnomaly(oe_c_mean_j2_series(6, iter), oe_c_mean_j2_series(2, iter));
    oe_d_mean_j2_series(6, iter) = MeanToTrueAnomaly(oe_d_mean_j2_series(6, iter), oe_d_mean_j2_series(2, iter));

    QNS_oe_c_mean_series(:, iter) = OE2QNS_OE(oe_c_mean_series(:, iter));
    QNS_oe_d_mean_series(:, iter) = OE2QNS_OE(oe_d_mean_series(:, iter));
    QNS_oe_c_mean_j2_series(:, iter) = OE2QNS_OE(oe_c_mean_j2_series(:, iter));
    QNS_oe_d_mean_j2_series(:, iter) = OE2QNS_OE(oe_d_mean_j2_series(:, iter));
    
    % d. Mean relative quasi-non-singular orbital elements
    QNS_roe_mean_series(:, iter) = OE2ROE(oe_c_mean_series(:, iter), oe_d_mean_series(:, iter));
    QNS_roe_mean_j2_series(:, iter) = OE2ROE(oe_c_mean_j2_series(:, iter), oe_d_mean_j2_series(:, iter));
end

% plotting qns oe
ylabels = ["a [km]", "u [rad]", "e_x", "e_y", "i [rad]", "RAAN [rad]"];
figure(3);
PlotOEvsTime(tspan, QNS_oe_c_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_c_mean_series, ylabels);
legend("Osculating", "Mean");
sgtitle("Chief quasi non-singular mean and osculating orbital elements vs time, without J2");
figure(4);
PlotOEvsTime(tspan, QNS_oe_d_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_d_mean_series, ylabels);
legend("Osculating", "Mean");
sgtitle("Deputy quasi non-singular mean and osculating orbital elements vs time, without J2");

figure(5);
PlotOEvsTime(tspan, QNS_oe_c_j2_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_c_mean_j2_series, ylabels);
legend("Osculating", "Mean");
sgtitle("Chief quasi non-singular mean and osculating orbital elements vs time, with J2");
figure(6);
PlotOEvsTime(tspan, QNS_oe_d_j2_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_d_mean_j2_series, ylabels);
legend("Osculating", "Mean");
sgtitle("Deputy quasi non-singular mean and osculating orbital elements vs time, with J2");

% plotting qns roe
ylabels = ["\delta a", "\delta \lambda", "\delta e_x", "\delta e_y", "\delta i_x", "\delta i_y"];
figure(7);
PlotOEvsTime(tspan, QNS_roe_series, ylabels);
PlotOEvsTime(tspan, QNS_roe_mean_series, ylabels);
legend("Osculating", "Mean");
sgtitle("Relative quasi non-singular mean and osculating relative orbital elements vs time, without J2");

figure(8);
PlotOEvsTime(tspan, QNS_roe_j2_series, ylabels);
PlotOEvsTime(tspan, QNS_roe_mean_j2_series, ylabels);
legend("Osculating", "Mean");
sgtitle("Relative quasi non-singular mean and osculating relative orbital elements vs time, with J2");

%% 4) RTN plots
% figures 9-12
PlotRTNSpace_meters(1000 * [r_RTN_abssim; v_RTN_abssim]');
figure(9);
sgtitle("3D Relative position in RTN, without J2");
figure(10);
sgtitle("Planar Relative position in RTN, without J2");

PlotRTNSpace_meters(1000 * [r_RTN_abssim_j2; v_RTN_abssim_j2]');
figure(11);
sgtitle("3D Relative position in RTN, with J2");
figure(12);
sgtitle("Planar Relative position in RTN, with J2");

%% 5) QNS plots
% figures 13-15
figure(13);
PlotQNSROE_meters(QNS_roe_series, a_c*1000);
PlotQNSROE_meters(QNS_roe_mean_series, a_c*1000);
sgtitle("Relative Motion, without J2");

figure(14);
PlotQNSROE_meters(QNS_roe_j2_series, a_c*1000);
PlotQNSROE_meters(QNS_roe_mean_j2_series, a_c*1000);
sgtitle("Relative Motion, with J2");


%% 6) 
