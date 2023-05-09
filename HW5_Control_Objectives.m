% Paxton Scott and Chris Osgood
% 5/07/2023
% AA279D HW5

clc; clear; close all;
addpath(genpath("Functions/"));

% Models:
%   1. Absolute non-linear, non-circular absolute model for chief. 
%   2. Relative linear (first order), circular ROE STM with J2. 


%% 1) Initilize orbital elements, constants, and simulation parameters.
mu = 3.986e5;

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

% deputy orbit
roe_d = [0; 0.100; 0; 0.030; 0; 0.030] / a_c; % [m]
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

% simulation parameters.
n_orbits = 200;
n_steps_per_orbit = 60;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);


%% 2) Absolute orbit propagator.
% Absolute non-linear sim of chief. 
[~, x_c_ECI] = ode113(@AbsoluteOrbitWithJ2DiffEq, tspan, x_c_osc_0, options);
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


%% 3) Least squares control solution. 
u_burns = [0, 1, 2, 3];
roe_i = [0, 100, 0, 0, 0, 0];
roe_f = [0, 100, 0, 30, 0, 30];

delta_vs = LS_control_solve(oe_c, roe_f - roe_i, u_burns);

%% 4) STM linear model for Quasi Nonsingular ROE with J2
QNS_roe_d_series_STM = zeros(6, n_iter);
oe_c_series = zeros(6, n_iter);
oe_c_mean_series = zeros(6, n_iter);

% Set intial roe with the intial deputy roe. 
QNS_roe_d_series_STM(:, 1) = roe_d;

% propagate STM
for iter = 1:n_iter
    oe_c_series(:, iter) = ECI2OE(r_c_ECI(:, iter), v_c_ECI(:, iter))';
    oe_c_mean_series(:, iter) = osc2mean(oe_c_series(:, iter), 1);

    if iter < n_iter
        QNS_roe_d_series_STM(:, iter+1) = STM_QNS_ROE_J2(oe_c, QNS_roe_d_series_STM(:,iter), tspan(iter));
    end
    
    
    % Apply Reconfiguration Manuever
    % if current time is a pre-defined manuever time
    %   ApplyDeputyManuever(oe_c_mean_series(:,iter), QNS_roe_d_series(:,iter+1), dv
    
    % Apply Formation Keeping Control
    % if dex has gone out of acceptable band (for our case it will go too
    % positive)
    %   then calculate manuever to push dex to other side of
    %   acceptable band (some pre-defined negative dex)
    %   apply such maneuver following same logic as reconfig manuever above
end

figure(2);
PlotQNSROE_meters(QNS_roe_d_series_STM, a_c*1000);
subplot(3,1,1);
sgtitle("Mean relative orbital elements of deputy, with J2, STM");

