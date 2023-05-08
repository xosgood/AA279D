% Paxton Scott and Chris Osgood
% 5/07/2023
% AA279D HW5

clc; clear; close all;
addpath(genpath("Functions/"));

% Implement and verify the models for your mission in open loop.
% Models:
%   1. Absolute non-linear, non-circular absolute model for chief. 
%   2. Relative linear (first order), circular ROE STM with J2. 
%   3. Relative non-linear ROE with dif. eqs.
%
% constants
mu = 3.986e5;


%% 1) Absolute orbit propagator.
% chief mean OE
a_c = 7000; % km
e_c = 0.001;
i_c = deg2rad(98);
RAAN_c = 0;
omega_c = deg2rad(90);
nu_c = 0;
M_c = TrueToMeanAnomaly(nu_c, e_c);
oe_c = [a_c; e_c; i_c; RAAN_c; omega_c; nu_c];

oe_c_osc = mean2osc(oe_c, 1);

[r_c_osc_0_ECI, v_c_osc_0_ECI] = OE2ECI(oe_c_osc);

x_c_osc_0 = [r_c_osc_0_ECI; v_c_osc_0_ECI];

% simulation parameters
n_orbits = 15;
n_steps_per_orbit = 100;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

% Absolute non-linear sim of chief. 
[~, x_c_ECI] = ode113(@func_J2, tspan, x_c_osc_0, options);
x_c_ECI = x_c_ECI';

r_c_ECI = x_c_ECI(1:3,:);
v_c_ECI = x_c_ECI(4:6,:);

% plot orbits in ECI
figure(1);
PlotEarth();
plot3(r_c_ECI(1,:), r_c_ECI(2,:), r_c_ECI(3,:), 'r');
title("Absolute ECI Orbits, with J2");
legend("Earth", "Chief","Location", "best");
xlabel('I (km)'); ylabel('J (km)'); zlabel('K (km)');

%% 2) STM linear model for Quasi Nonsingular ROE with J2

% deputy relative and absolute orbit elements
roe_d = [0; 0.100; 0.050; 0.100; 0.0; 0.200] / a_c ; % [m]
oe_d = ROE2OE(oe_c, roe_d); % deputy mean oe
a_d = oe_d(1);
e_d = oe_d(2);
i_d = oe_d(3);
RAAN_d = oe_d(4);
omega_d = oe_d(5);
nu_d = oe_d(6);
M_d = TrueToMeanAnomaly(nu_d, e_d);

oe_d_mean = oe_d;
oe_d_osc = mean2osc(oe_d_mean, 1);

QNS_roe_d_series_STM = zeros(6, n_iter);
oe_c_series = zeros(6, n_iter);
oe_c_mean_series = zeros(6, n_iter);

% time step.
dt = tspan(2) - tspan(1);

% Set intial roe with the intial deputy roe. 
QNS_roe_d_series_STM(:, 1) = roe_d;

for iter = 1:n_iter
    oe_c_series(:, iter) = ECI2OE(r_c_ECI(:, iter), v_c_ECI(:, iter))';

    oe_c_mean_series(:, iter) = osc2mean(oe_c_series(:, iter), 1);

    if iter < n_iter
        QNS_roe_d_series_STM(:, iter+1) = STM_QNS_ROE_J2(oe_c_series(:,iter), QNS_roe_d_series_STM(:,iter), dt);
    end
end

figure(2);
PlotQNSROE_meters(QNS_roe_d_series_STM, a_c*1000);
subplot(3,1,1);
legend("Relative orbital elements of deputy" ,"Location", "best");
sgtitle("Relative Motion, with J2, STM");


%% 3) Non linear relative orbit propagator of deputy.
QNS_roe_series = zeros(6, n_iter);
QNS_roe_mean_series = zeros(6, n_iter);
oe_d_series = zeros(6, n_iter);
oe_d_mean_series = zeros(6, n_iter);

oe_d_osc = mean2osc([a_d; e_d; i_d; RAAN_d; omega_d; nu_d], 1);

[r_d_osc_0_ECI, v_d_osc_0_ECI] = OE2ECI(oe_d_osc);
x_d_osc_0 = [r_d_osc_0_ECI; v_d_osc_0_ECI];

[~, x_d_ECI] = ode113(@func_J2, tspan, x_d_osc_0, options);
x_d_ECI = x_d_ECI';

r_d_ECI = x_d_ECI(1:3,:);
v_d_ECI = x_d_ECI(4:6,:);

[r_RTN, v_RTN] = ECI2RTN_Vectorized(r_c_ECI, v_c_ECI, r_d_ECI, v_d_ECI);

PlotRTNSpace_meters(1000 * [r_RTN; v_RTN]');

for iter = 1:size(r_c_ECI,2)
    oe_d_series(:, iter) = ECI2OE(r_d_ECI(:, iter), v_d_ECI(:, iter))';
    oe_d_mean_series(:, iter) = osc2mean(oe_d_series(:, iter), 1);
    QNS_roe_series(:, iter) = OE2ROE(oe_c_series(:, iter), oe_d_series(:, iter));  
    QNS_roe_mean_series(:, iter) = OE2ROE(oe_c_mean_series(:, iter), oe_d_mean_series(:, iter));
end

figure(5);
PlotQNSROE_meters(QNS_roe_series, a_c*1000);
PlotQNSROE_meters(QNS_roe_mean_series, a_c*1000);
subplot(3,1,1);
legend("Relative orbital elements of deputy" ,"Location", "best");
sgtitle("Relative Motion, with J2, Non-linear");


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

% Dif. Eqs.
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