% Paxton Scott and Chris Osgood
% 4/21/23
% AA279D HW3 P2(h)(i)
clc; clear; close all;
addpath(genpath("Functions/"));

% constants
mu = 3.986e5;

%% a) initial orbits
% chief
a_0 = 9000; % km
e_0 = 0.2;
i_0 = deg2rad(98);
RAAN_0 = 0;
omega_0 = deg2rad(90);
nu_0 = 0;
oe_c = [a_0; e_0; i_0; RAAN_0; omega_0; nu_0];

% deputy
a_1 = a_0 + 100;
e_1 = 0.1999;
i_1 = deg2rad(98.01);
RAAN_1 = deg2rad(0.05);
omega_1 = deg2rad(90.05);
nu_1 = deg2rad(0.01);
oe_d = [a_1; e_1; i_1; RAAN_1; omega_1; nu_1];

%% b) use initial conditions to find YA solution to TH equations
[r_ECI_0, v_ECI_0] = OE2ECI(a_0, e_0, i_0, RAAN_0, omega_0, nu_0);
[r_ECI_1, v_ECI_1] = OE2ECI(a_1, e_1, i_1, RAAN_1, omega_1, nu_1);

[r_RTN_d, v_RTN_d] = ECI2RTN(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

rho = norm(r_RTN_d);
r0 = norm(r_ECI_0);

x_RTN_d = [r_RTN_d; v_RTN_d];

K = RTN2YA_IC(x_RTN_d, a_0, e_0, nu_0, 0);

%% c) simulate YA solution
% sim parameters
n_orbits = 15;
n_steps_per_orbit = 30;
n_iter = n_orbits * n_steps_per_orbit;
nu_step = 2 * pi / n_steps_per_orbit;
T = 2 * pi * sqrt(a_0^3 / mu);

t = zeros(1, n_iter);
x_RTN_YA = zeros(6, n_iter);

orbits_passed = -1;
for iter = 1:n_iter % for each orbit
    if mod(iter, n_steps_per_orbit) == 1
        orbits_passed = orbits_passed + 1;
    end
    % compute true anomaly for this timestep
    nu = iter * nu_step;
    % back out time from true anomaly
    t(iter) = TrueAnomalyToTime(nu, a_0, e_0) + orbits_passed * T;
    % propagate using YA solution
    x_RTN_YA(:,iter) = YA2RTN(K, a_0, e_0, nu, t(iter));
end

figure(5);
PlotRTN(t, x_RTN_YA(1:3,:), x_RTN_YA(4:6,:));
sgtitle("YA solution in RTN")

PlotRTNSpaceMultiple(x_RTN_YA', false);

%% e) compute QNS relative OEs
% roe = [da, dlambda, dex, dey, dix, diy]^T
roe = OE2ROE(oe_c, oe_d);

%% f) simulate with linearized geometric mapping (LGM)
M_c_0 = TrueToMeanAnomaly(nu_0, e_0);
M_d_0 = TrueToMeanAnomaly(nu_1, e_1);
dM_0 = M_d_0 - M_c_0; % delta mean anomaly at time 0

t = zeros(1, n_iter);
x_RTN_LGM = zeros(6, n_iter);

orbits_passed = -1;
for iter = 1:n_iter % for each orbit
    if mod(iter, n_steps_per_orbit) == 1
        orbits_passed = orbits_passed + 1;
    end
    % compute true anomaly for this timestep
    nu = iter * nu_step;
    % back out time from true anomaly
    t(iter) = TrueAnomalyToTime(nu, a_0, e_0) + orbits_passed * T;
    % update chief oe with new true anomaly
    oe_c_cur = oe_c;
    oe_c_cur(6) = nu;
    % propagate using linearized geometric mapping solution
    x_RTN_LGM(:,iter) = GeometricMapping_Linear(roe, oe_c_cur, t(iter));
end

figure(6);
PlotRTN(t, x_RTN_LGM(1:3,:), x_RTN_LGM(4:6,:));
sgtitle("Linear geometric mapping solution in RTN");

PlotRTNSpaceMultiple(x_RTN_LGM', false);

%% g) compare ROEs and integration constants
roe
roe_from_YAIC = YAIC2ROE(K, oe_c)

%% h) Produce the true relative position and velocity from analytical or numerical propogation.
tspan = t;

v_ECI_c_expressedInRTN = rECI2RTN([r_ECI_0; v_ECI_0]) * v_ECI_0;

r0 = norm(r_ECI_0);
r0_dot = v_ECI_c_expressedInRTN(1);

theta = nu_0 + omega_0;
theta_dot = norm(cross(r_ECI_0, v_ECI_0))/r0^2;
state0 = [r_RTN_d; theta; r0; v_RTN_d; theta_dot; r0_dot];

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[t, state_out] = ode113(@RelativeMotionDifEqRTN, tspan, state0, options);

x_RTN = [state_out(:,1:3), state_out(:,6:8)];

figure(7);
PlotRTN(t, state_out(:,1:3)', state_out(:,6:8)');
sgtitle("Full nonlinear solution in RTN");

PlotRTNSpaceMultiple(x_RTN, true);

% Propertly label the lines. 
for iter = 1:4
    figure(iter);
    legend({'YA solution', 'Linear Geometric Mapping solution', 'Non-linear numerical solution'})
end

hold off;

% plot errors
error_YA = abs(x_RTN_YA - x_RTN');
error_LGM = abs(x_RTN_LGM - x_RTN');
figure(8);
PlotRTN(t, error_YA(1:3,:), error_YA(4:6,:));
sgtitle("Error in RTN between YA solution and nonlinear solution")
figure(9);
PlotRTN(t, error_LGM(1:3,:), error_LGM(4:6,:));
sgtitle("Error in RTN between linear geometric mapping solution and nonlinear solution")

