% Chris Osgood
% 5/30/23
% AA279D HW8
% UKF on Relative Motion
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
R = diag([10, 10, 10, .02, .02, .02, 10, 10, 10, .02, .02, .02] / 1e5); % eye(m_dims_meas);

x_absolute = zeros(m_dims_meas, n_iter); % true absolute state (osculating), through simulating dynamics
x_roe = zeros(n_dims_state, n_iter); % true relative state (osculating), through simulating dynamics
y = zeros(m_dims_meas, n_iter-1); % measurements (osculating)

u = zeros(p_dims_controlinput, n_iter-1); % control inputs (delta-vs in RTN, will be populated throughout propagation)

mu = zeros(n_dims_state, n_iter); % state estimate (estimated osculating roe)
Sigma = zeros(n_dims_state, n_dims_state, n_iter); % covariance estimate

x_absolute(:,1) = [x_c_osc_0(:); x_d_osc_0(:)]; % initial truth state
x_roe(:,1) = roe_d_osc(:); % initial truth state

mu_0 = zeros(n_dims_state, 1); % initial state estimate
Sigma_0 = 0.1 * eye(n_dims_state); % initial covariance estimate



for iter = 2:n_iter
    %%% ground truth (in ECI) propagation
    % apply control to ground truth
    dx_RTN = [0; 0; 0; u(:,iter-1)];
%     x_and_dv_ECI = RTN2ECI(x_absolute(7:12,iter-1), dx_RTN);
%     dv_ECI = x_and_dv_ECI(4:6);
%     x_absolute(10:12,iter) = x_absolute(10:12,iter-1) + dv_ECI;
    x_absolute(7:12,iter) = RTN2ECI(x_absolute(7:12,iter-1), dx_RTN);
    
    %w = mvnrnd(zeros(n_dims_state, 1), Q)'; % process noise
    v = mvnrnd(zeros(m_dims_meas, 1), R)'; % measurement noise
    
    [~, x_c_ECI] = ode113(@AbsoluteOrbitWithJ2DiffEq, tspan(iter-1:iter), x_absolute(1:6,iter-1), options);
    [~, x_d_ECI] = ode113(@AbsoluteOrbitWithJ2DiffEq, tspan(iter-1:iter), x_absolute(7:12,iter), options);
    x_absolute(1:6,iter) = x_c_ECI(end,:);
    x_absolute(7:12,iter) = x_d_ECI(end,:);
    
    oe_c_osc_cur = ECI2OE(x_absolute(1:3,iter), x_absolute(4:6,iter));
    oe_d_osc_cur = ECI2OE(x_absolute(7:9,iter), x_absolute(10:12,iter));
    
    x_roe(:,iter) = OE2ROE(oe_c_osc_cur, oe_d_osc_cur);
    
    oe_c_mean_cur = osc2mean(oe_c_osc_cur, 1);
    oe_d_mean_cur = osc2mean(oe_d_osc_cur, 1);
    roe_d_mean_cur = OE2ROE(oe_c_mean_cur, oe_d_mean_cur);
    
    % measure
    y(:,iter-1) = x_absolute(:,iter) + v;
    
    
    %%% UKF
    % predict (regular KF style, since our dynamics are linear)
    A = DynamicsJacobian(oe_c_osc_cur, dt); % Jacobian for dynamics
    mu(:,iter) = f(mu(:,iter-1), u(:,iter-1), oe_c_osc_cur, dt); % mean predict
    Sigma(:,:,iter) = A * Sigma(:,:,iter-1) * A' + Q; % covariance predict
%     % predict (UKF style)
%     [points, weights] = UT(mu(:,i-1), Sigma(:,:,i-1));
%     for j = 1:size(points, 2) % loop through sigma-points
%         points(:,j) = f(points(:,j), u(:,i-1), dt); % propagate points through nonlinear dynamics
%     end
%     [mu(:,i), Sigma(:,:,i)] = UTInverse(points, weights); % get predicted mean and covariance
%     Sigma(:,:,i) = Sigma(:,:,i) + Q; % add process noise to covariance
    % update
    [points, weights] = UT(mu(:,iter), Sigma(:,:,iter)); % update sigma-points using predicted mean and covariance
    y_preds = zeros(m_dims_meas, size(points, 2)); % create vector to hold predicted measurements at each sigma-point
    for j = 1:size(points, 2) % loop through sigma-points and measure
        y_preds(:,j) = g(points(:,j), oe_c_osc_cur); % predicted measurements at each sigma point
    end
    [y_hat, Sigma_yy] = UTInverse(y_preds, weights); % get expected measurement and expected measurement covariance
    Sigma_yy = Sigma_yy + R; % add measurement noise to covariance
    Sigma_xy = (weights .* (points - mu(:,iter))) * (y_preds - y_hat)'; % cross covariance between sigma-points and predicted measurements
    K = Sigma_xy / Sigma_yy; % "Kalman gain"
    mu(:,iter) = mu(:,iter) + K * (y(:,iter-1) - y_hat); % mean update
    Sigma(:,:,iter) = Sigma(:,:,iter) - K * Sigma_xy'; % covariance update
    
    %%% compute control for next timestep to use
    u(:,iter) = LyapunovController(roe_d_mean_cur, roe_desired, oe_c_mean_cur, lyap_params);
end



figure; grid on; hold on;
plot(tspan, x_roe(1,:));
plot(tspan, mu(1,:));

figure; grid on; hold on;
plot(tspan, x_roe(2,:));
plot(tspan, mu(2,:));


%% functions
% nonlinear dynamics
function x_new = f(x_old, u, oe_c_osc, dt)
    B = GenerateBControlsMatrix(osc2mean(oe_c_osc,1));

    oe_c_mean = osc2mean(oe_c_osc, 1);
    oe_d_osc = ROE2OE(oe_c_osc, x_old);
    oe_d_mean = osc2mean(oe_d_osc, 1);
    roe_d_mean = OE2ROE(oe_c_mean, oe_d_mean);
    
    A = DynamicsJacobian(oe_c_osc, dt);
    x_new = A * roe_d_mean + B*u;
end


% nonlinear measurement function
function y = g(x, oe_c_osc)
    % map roe of deputy and oe of chief to ECI of chief and deputy.
    % x is roe of deputy.
    oe_d_osc = ROE2OE(oe_c_osc, x);
    [r_c_ECI, v_c_ECI] = OE2ECI(oe_c_osc);
    [r_d_ECI, v_d_ECI] = OE2ECI(oe_d_osc);

    y = [r_c_ECI; v_c_ECI; r_d_ECI; v_d_ECI];
end

% generate Jacobian for dynamics
function A = DynamicsJacobian(oe_c_osc, dt)
    A = GenerateSTMJ2(oe_c_osc, dt);
end

