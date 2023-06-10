% Chris Osgood
% 6/10/23
% AA279D HW8
% Integrating Navigation and Control
clc; clear; close all;
addpath(genpath("Functions/"));

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
oe_d = ROE2OE(oe_c_mean, roe_d_mean);
a_d = oe_d(1);
e_d = oe_d(2);
i_d = oe_d(3);
RAAN_d = oe_d(4);
omega_d = oe_d(5);
nu_d = oe_d(6);
%M_d = TrueToMeanAnomaly(nu_d, e_d);
% deputy osculating OE
oe_d_mean = oe_d;
oe_d_osc = mean2osc(oe_d_mean, 1);
roe_d_osc = OE2ROE(oe_c_osc, oe_d_osc);
% get deputy initial state for ode113
[r_d_osc_0_ECI, v_d_osc_0_ECI] = OE2ECI(oe_d_osc);
x_d_osc_0 = [r_d_osc_0_ECI; v_d_osc_0_ECI];

%% Lyapunov controller setup
roe_desired = [0; 100; 0; 30; 0; 30] / a_c;
N = 10; % exponent in P matrix
k = 1; % (1/k) factor in front of P matrix
u_lowerbound = 1e-6; % lower bound on control actuation
u_upperbound = 1e-3; % upper bound on control actuation
dlambda_thresh = 0.1;
dlambda_dot = 0.001;
lyap_params = [N, k, u_lowerbound, u_upperbound, dlambda_thresh, dlambda_dot];

%% UKF setup
n_dims_state = 6; % number of dimensions of state
m_dims_meas = 12; % number of dimensions of measurement
p_dims_controlinput = 3; % number of dimensions of control input

n_orbits = 10;
n_steps_per_orbit = 30;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
dt = tspan(2) - tspan(1);
orbit_span = (1:n_iter)/n_steps_per_orbit;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

Q = 0.000001 * eye(n_dims_state);
R = diag([10, 10, 10, .02, .02, .02, 10, 10, 10, .02, .02, .02] / 1e3); % eye(m_dims_meas);

x_absolute = zeros(m_dims_meas, n_iter); % true absolute state (osculating), through simulating dynamics
x_roe = zeros(n_dims_state, n_iter); % true relative state (osculating), through simulating dynamics
y = zeros(m_dims_meas, n_iter-1); % measurements (osculating)
measured_roe = zeros(n_dims_state, n_iter); % the roe values based on measurment.

u = zeros(p_dims_controlinput, n_iter-1); % control inputs (delta-vs in RTN, will be populated throughout propagation)

mu = zeros(n_dims_state, n_iter); % state estimate (estimated osculating roe)
Sigma = zeros(n_dims_state, n_dims_state, n_iter); % covariance estimate

x_absolute(:,1) = [x_c_osc_0(:); x_d_osc_0(:)]; % initial truth state
x_roe(:,1) = roe_d_osc(:); % initial truth state

mu_0 = zeros(n_dims_state, 1); % initial state estimate
Sigma_0 = 0.1 * eye(n_dims_state); % initial covariance estimate

Sigma(:,:,1) = Sigma_0;
mu(:,1) = mu_0;

