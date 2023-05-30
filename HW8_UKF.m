% Chris Osgood
% 5/30/23
% AA279D HW8
% UKF on Relative Motion
clc; clear; close all;
addpath("../Functions/");

rng(273);

%% relative motion setup
% constants
mu = 3.986e5;

%%% initialize chief
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

%%% initialize deputy
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

%% UKF setup

n_dims_state = 6; % number of dimensions of state
m_dims_meas = 12; % number of dimensions of measurement

n_orbits = 30;
n_steps_per_orbit = 100;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
dt = tspan(2) - tspan(1);
orbit_span = (1:n_iter)/n_steps_per_orbit;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

Q = 0.1 * dt * eye(n_dims_state);
R = 0.1 * eye(m_dims_meas);

x_absolute = zeros(m_dims_meas, n_iter); % true absolute state, through simulating dynamics
x_roe = zeros(n_dims_state, n_iter); % true relative state, through simulating dynamics
y = zeros(m_dims_meas, n_iter-1); % measurements

mu = zeros(n_dims_state, n_iter);
Sigma = zeros(n_dims_state, n_dims_state, n_iter);

x_absolute(:,1) = [x_c_osc_0(:); x_d_osc_0(:)];
x_roe(:,1) = roe_d_osc(:);

mu_0 = zeros(n_dims_state, 1); % initial state estimate
Sigma_0 = 0.1 * eye(n_dims_state); % initial covariance estimate



for i = 2:n_iter
    w = mvnrnd(zeros(n_dims_state, 1), Q)'; % process noise
    v = mvnrnd(zeros(m_dims_state, 1), R)'; % measurement noise
    x_absolute(:,i) = f(x(:,i-1), u(:,i-1), dt); %+ w; % dynamics propagation
    y(:,i-1) = g(x(:,i)) + v; % measure
end






%% functions
% nonlinear dynamics
function x_new = f(x_old, u, dt)
    x_new = zeros(size(x_old));
    
end

% nonlinear measurement
function y = g(x)
    y = norm(x(1:2));
end

% generate Jacobian for dynamics
function A = DynamicsJacobian(x, u, dt)
    A = eye(length(x));
    A(1,3) = -dt * u(1) * sin(x(3));
    A(2,3) = dt * u(1) * cos(x(3));
end

% generate Jacobian for measurements
function C = MeasurementJacobian(x)
    p_norm = norm(x(1:2));
    C = [x(1) / p_norm, x(2) / p_norm, 0];
end