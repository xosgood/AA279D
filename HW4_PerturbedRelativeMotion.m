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

[r_c_0_ECI, v_c_0_ECI] = OE2ECI(a_c, e_c, i_c, RAAN_c, omega_c, nu_c);

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

[r_d_0_ECI, v_d_0_ECI] = OE2ECI(a_d, e_d, i_d, RAAN_d, omega_d, nu_d);

%% 3) full non-linear simulation
%%% simulate
% sim parameters
n_orbits = 15;
n_steps_per_orbit = 100;
n_iter = n_steps_per_orbit * n_orbits;
mu = 3.986e5; % (km^3 / s^2)
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter) ; % [0, t_f];

x_c_0 = [r_c_0_ECI; v_c_0_ECI];
x_d_0 = [r_d_0_ECI; v_d_0_ECI];

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
[~, x_c_ECI] = ode113(@func, tspan, x_c_0, options);
[~, x_d_ECI] = ode113(@func, tspan, x_d_0, options);
[~, x_c_j2_ECI] = ode113(@func_J2, tspan, x_c_0, options);
[t, x_d_j2_ECI] = ode113(@func_J2, tspan, x_d_0, options);

x_c_ECI = x_c_ECI';
x_d_ECI = x_d_ECI';
x_c_j2_ECI = x_c_j2_ECI';
x_d_j2_ECI = x_d_j2_ECI';

r_c_ECI = x_c_ECI(1:3,:);
v_c_ECI = x_c_ECI(4:6,:);
r_d_ECI = x_d_ECI(1:3,:);
v_d_ECI = x_d_ECI(4:6,:);
r_c_j2_ECI = x_c_j2_ECI(1:3,:);
v_c_j2_ECI = x_c_j2_ECI(4:6,:);
r_d_j2_ECI = x_d_j2_ECI(1:3,:);
v_d_j2_ECI = x_d_j2_ECI(4:6,:);

% plot orbits in ECI
figure(1);
PlotEarth();
plot3(r_c_ECI(1,:), r_c_ECI(2,:), r_c_ECI(3,:), 'r');
plot3(r_d_ECI(1,:), r_d_ECI(2,:), r_d_ECI(3,:), 'g');
title("Absolute ECI Orbits, without J2");
legend("Earth", "Chief", "Deputy", "Location", "best");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

% plot orbits in ECI with the J2 effect
figure(2);
PlotEarth();
plot3(r_c_j2_ECI(1,:), r_c_j2_ECI(2,:), r_c_j2_ECI(3,:), 'r');
plot3(r_d_j2_ECI(1,:), r_d_j2_ECI(2,:), r_d_j2_ECI(3,:), 'g');
title("Absolute ECI Orbits, with J2");
legend("Earth", "Chief", "Deputy", "Location", "best");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');


%%% Conversions
% convert to RTN
[r_RTN, v_RTN] = ECI2RTN_Vectorized(r_c_ECI, v_c_ECI, r_d_ECI, v_d_ECI);
[r_j2_RTN, v_j2_RTN] = ECI2RTN_Vectorized(r_c_j2_ECI, v_c_j2_ECI, r_d_j2_ECI, v_d_j2_ECI);

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

for iter = 1:size(r_c_ECI,2)
    oe_c_series(:, iter) = ECI2OE(r_c_ECI(:, iter), v_c_ECI(:, iter))';
    oe_d_series(:, iter) = ECI2OE(r_d_ECI(:, iter), v_d_ECI(:, iter))';
    oe_c_j2_series(:, iter) = ECI2OE(r_c_j2_ECI(:, iter), v_c_j2_ECI(:, iter))';
    oe_d_j2_series(:, iter) = ECI2OE(r_d_j2_ECI(:, iter), v_d_j2_ECI(:, iter))';

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

%%% Plotting
% plotting qns oe
ylabels = ["a [km]", "u [rad]", "e_x", "e_y", "i [rad]", "RAAN [rad]"];
figure(3);
PlotOEvsTime(tspan, QNS_oe_c_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_c_mean_series, ylabels);
subplot(6,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Chief quasi non-singular mean and osculating orbital elements vs time, without J2");
figure(4);
PlotOEvsTime(tspan, QNS_oe_d_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_d_mean_series, ylabels);
subplot(6,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Deputy quasi non-singular mean and osculating orbital elements vs time, without J2");

figure(5);
PlotOEvsTime(tspan, QNS_oe_c_j2_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_c_mean_j2_series, ylabels);
subplot(6,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Chief quasi non-singular mean and osculating orbital elements vs time, with J2");
figure(6);
PlotOEvsTime(tspan, QNS_oe_d_j2_series, ylabels);
PlotOEvsTime(tspan, QNS_oe_d_mean_j2_series, ylabels);
subplot(6,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Deputy quasi non-singular mean and osculating orbital elements vs time, with J2");

% plotting qns roe
ylabels = ["a \delta a [m]", "a \delta \lambda [m]", "a \delta e_x [m]", "a \delta e_y [m]", "\delta i_x [m]", "\delta i_y [m]"];
figure(7);
PlotOEvsTime(tspan, 1000 * oe_c_series(1,:) .* QNS_roe_series, ylabels);
PlotOEvsTime(tspan, 1000 * oe_c_mean_series(1,:) .* QNS_roe_mean_series, ylabels);
subplot(6,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Relative quasi non-singular mean and osculating relative orbital elements vs time, without J2");

figure(8);
PlotOEvsTime(tspan, 1000 * oe_c_j2_series(1,:) .* QNS_roe_j2_series, ylabels);
PlotOEvsTime(tspan, 1000 * oe_c_mean_j2_series(1,:) .* QNS_roe_mean_j2_series, ylabels);
subplot(6,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Relative quasi non-singular mean and osculating relative orbital elements vs time, with J2");

%% 4) RTN plots
% figures 9-12
PlotRTNSpace_meters(1000 * [r_RTN; v_RTN]');
figure(9);
sgtitle("3D Relative position in RTN, without J2");
figure(10);
sgtitle("Planar Relative position in RTN, without J2");

PlotRTNSpace_meters(1000 * [r_j2_RTN; v_j2_RTN]');
figure(11);
sgtitle("3D Relative position in RTN, with J2");
figure(12);
sgtitle("Planar Relative position in RTN, with J2");

%% 5) QNS plots
% figures 13-15
figure(13);
PlotQNSROE_meters(QNS_roe_series, a_c*1000);
PlotQNSROE_meters(QNS_roe_mean_series, a_c*1000);
subplot(3,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Relative Motion, without J2");

figure(14);
PlotQNSROE_meters(QNS_roe_j2_series, a_c*1000);
PlotQNSROE_meters(QNS_roe_mean_j2_series, a_c*1000);
subplot(3,1,1);
legend("Osculating", "Mean", "Location", "best");
sgtitle("Relative Motion, with J2");

%% 6) 

%% 7)
% new deputy ROE (quasi non-singular) [da, dlambda, dex, dey, dix, diy]^T
roe_d = [0; 0.100; 0.050; 0.100; 0.000; 0.200] / a_c ; % [m]

oe_d = ROE2OE(oe_c, roe_d);
a_d = oe_d(1);
e_d = oe_d(2);
i_d = oe_d(3);
RAAN_d = oe_d(4);
omega_d = oe_d(5);
nu_d = oe_d(6);

[r_d_0_ECI, v_d_0_ECI] = OE2ECI(a_d, e_d, i_d, RAAN_d, omega_d, nu_d);


%% ODE Functions
function statedot = func(t, state)
    % State vector is [rx ry rz vx vy vz]’.
    % Although required by form, input value t will go unused here.
    mu = 3.986e5; % (km^3 / s^2)
    r = state(1:3);
    rdot = state(4:6);
    vdot = CentralBodyAccel(mu, r);
    statedot = [rdot;
                vdot];
end

function statedot = func_J2(t, state)
    % State vector is [rx ry rz vx vy vz]’.
    % Although required by form, input value t will go unused here.
    mu = 3.986e5; % (km^3 / s^2) for earth
    R_E = 6378.1; % equatorial radius of earth in km
    J2 = 0.108263e-2; % for earth
    r = state(1:3);
    rdot = state(4:6);
    vdot = CentralBodyAccel(mu, r) + J2AccelECI(J2, mu, R_E, r(1), r(2), r(3));
    statedot = [rdot;
                vdot];
end