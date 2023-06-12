% Chris Osgood
% 6/10/23
% AA279D HW8
% Integrating Navigation and Control
clc; clear; close all;
addpath(genpath("Functions/"));

rng(279);

%% relative motion setup
% constants
mu_E = 3.986e5;

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
n_c = sqrt(mu_E / a_c^3);
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

%% sim setup
n_orbits = 15;
n_steps_per_orbit = 100;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu_E);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
dt = tspan(2) - tspan(1);
orbit_span = (1:n_iter)/n_steps_per_orbit;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

%% Lyapunov controller setup
roe_desired = [0; 100; 0; 30; 0; 30] / a_c;
lyap_params.N = 14; % exponent in P matrix
lyap_params.k = 6; % (1/k) factor in front of P matrix
lyap_params.u_lowerbound = 1e-6; % lower bound on control actuation
lyap_params.u_upperbound = 1e-2; % upper bound on control actuation
lyap_params.dlambda_thresh = 0.1 / a_c;
lyap_params.dlambda_dot = 0.005 / a_c;

%% UKF setup
n_dims_state = 6; % number of dimensions of state
m_dims_meas = 12; % number of dimensions of measurement
p_dims_controlinput = 3; % number of dimensions of control input

Q = (0.02 / a_c) * eye(n_dims_state);
R = diag([10, 10, 10, .02, .02, .02, 10, 10, 10, .02, .02, .02] / 1e3); % eye(m_dims_meas);

x_absolute = zeros(m_dims_meas, n_iter); % true absolute state (osculating), through simulating dynamics
x_roe = zeros(n_dims_state, n_iter); % true relative state (osculating), through simulating dynamics
y = zeros(m_dims_meas, n_iter-1); % measurements (osculating)
measured_roe = zeros(n_dims_state, n_iter-1); % the roe values based on measurment.

u = zeros(p_dims_controlinput, n_iter-1); % control inputs (delta-vs in RTN, will be populated throughout propagation)

mu = zeros(n_dims_state, n_iter); % state estimate (estimated osculating roe)
Sigma = zeros(n_dims_state, n_dims_state, n_iter); % covariance estimate

x_absolute(:,1) = [x_c_osc_0(:); x_d_osc_0(:)]; % initial truth state
x_roe(:,1) = roe_d_osc(:); % initial truth state

mu_0 = zeros(n_dims_state, 1); % initial state estimate
Sigma_0 = 0.01 * eye(n_dims_state); % initial covariance estimate

Sigma(:,:,1) = Sigma_0;
mu(:,1) = mu_0;

pre_fit_res = zeros(m_dims_meas, n_iter); % pre fit residual
post_fit_res = zeros(m_dims_meas, n_iter); % post fit residual


for iter = 2:n_iter
    %%% ground truth (in ECI) propagation
    % apply control to ground truth
    dx_RTN = [0; 0; 0; u(:,iter-1)];
    x_absolute(7:12,iter) = RTN2ECI(x_absolute(7:12,iter-1), dx_RTN);
    
    w = mvnrnd(zeros(n_dims_state, 1), Q)'; % process noise (currently only used for input to controller)
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
    
    % add noise for controller
    roe_d_mean_cur_noisy = roe_d_mean_cur + w;
    
    % measure
    y(:,iter-1) = x_absolute(:,iter) + v;

    oe_c_osc_mes = ECI2OE(y(1:3,iter-1), y(4:6,iter-1));
    oe_d_osc_mes = ECI2OE(y(7:9,iter-1), y(10:12,iter-1));

    measured_roe(:,iter-1) = OE2ROE(oe_c_osc_mes, oe_d_osc_mes);

    
    
    %%% UKF
    % predict (regular KF style, since our dynamics are linear)
    A = DynamicsJacobian(oe_c_osc_cur, dt); % Jacobian for dynamics
    mu(:,iter) = f(mu(:,iter-1), u(:,iter-1), oe_c_osc_cur, dt); % mean predict
    Sigma(:,:,iter) = A * Sigma(:,:,iter-1) * A' + Q; % covariance predict
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

    pre_fit_res(:,iter) = (y(:,iter-1) - y_hat);

    Sigma(:,:,iter) = Sigma(:,:,iter) - K * Sigma_xy'; % covariance update
    
    post_fit_res(:,iter) = y(:,iter-1) - g(mu(:,iter), oe_c_osc_cur);
    %%% compute control for next timestep to use
    u(:,iter) = LyapunovController(roe_d_mean_cur_noisy, roe_desired, oe_c_mean_cur, lyap_params);
end

%% plotting
figure;
sgtitle("True and Estimated ROEs vs. time")
subplot(3,2,1); grid on; hold on;
plot(orbit_span, a_c * x_roe(1,:));
plot(orbit_span, a_c * mu(1,:));
plot(orbit_span(2:end), a_c * measured_roe(1,:));
xlabel("time [s]");
ylabel("a\delta a [km]")
subplot(3,2,2); grid on; hold on;
plot(orbit_span, a_c * x_roe(2,:));
plot(orbit_span, a_c * mu(2,:));
plot(orbit_span(2:end), a_c * measured_roe(2,:));
legend("ground truth", "state estimate", "measured", "Location","Best");
xlabel("time [s]");
ylabel("a\delta \lambda [km]")
subplot(3,2,3); grid on; hold on;
plot(orbit_span, a_c * x_roe(3,:));
plot(orbit_span, a_c * mu(3,:));
plot(orbit_span(2:end), a_c * measured_roe(3,:));
xlabel("time [s]");
ylabel("a\delta e_x [km]")
subplot(3,2,4); grid on; hold on;
plot(orbit_span, a_c * x_roe(4,:));
plot(orbit_span, a_c * mu(4,:));
plot(orbit_span(2:end), a_c * measured_roe(4,:));
xlabel("time [s]");
ylabel("a\delta e_y [km]")
subplot(3,2,5); grid on; hold on;
plot(orbit_span, a_c * x_roe(5,:));
plot(orbit_span, a_c * mu(5,:));
plot(orbit_span(2:end), a_c * measured_roe(5,:));
xlabel("time [s]");
ylabel("a\delta i_x [km]")
subplot(3,2,6); grid on; hold on;
plot(orbit_span, a_c * x_roe(6,:));
plot(orbit_span, a_c * mu(6,:));
plot(orbit_span(2:end), a_c * measured_roe(6,:));
xlabel("time [s]");
ylabel("a\delta i_y [km]")

PlotConfidenceInterval(orbit_span, x_roe, mu, Sigma, a_c);
subplot(3,2,1); ylim([-70, 70]);
subplot(3,2,2); ylim([-500, 800]);
subplot(3,2,3); ylim([-60, 60]);
subplot(3,2,4); ylim([-20, 100]);
subplot(3,2,5); ylim([-40, 60]);
subplot(3,2,6); ylim([-20, 80]);
subplot(3,2,1); ylim([-40, 50]);
subplot(3,2,2); ylim([0, 200]);
subplot(3,2,3); ylim([-40, 40]);
subplot(3,2,4); ylim([-10, 70]);
subplot(3,2,5); ylim([-20, 20]);
subplot(3,2,6); ylim([0, 50]);

% Plot pre and post fit residuals.
figure; 
sgtitle("Pre and Post fit residuals verus time")
subplot(2, 1, 1); hold on; grid on;
plot(orbit_span, vecnorm(pre_fit_res(7:9,:)), "o");
plot(orbit_span, vecnorm(post_fit_res(7:9,:)), ".red");
legend("Pre fit residual", "Post fit residual");
xlabel("orbits")
ylabel("norm residual [km]")
ylim([-5, 20]);

subplot(2, 1, 2); hold on; grid on;
plot(orbit_span, vecnorm(pre_fit_res(10:12,:)), "o");
plot(orbit_span, vecnorm(post_fit_res(10:12,:)), ".red");

legend("Pre fit residual", "Post fit residual");
xlabel("orbits")
ylabel("norm residual [km/s]")
ylim([-0.01, 0.05]);

%% functions
% nonlinear dynamics
function x_new = f(x_old, u, oe_c_osc, dt)
    B = GenerateBControlsMatrix(osc2mean(oe_c_osc,1));

    oe_c_mean = osc2mean(oe_c_osc, 1);
    oe_d_osc = ROE2OE(oe_c_osc, x_old);
    oe_d_mean = osc2mean(oe_d_osc, 1);
    roe_d_mean = OE2ROE(oe_c_mean, oe_d_mean);
    
    A = DynamicsJacobian(oe_c_osc, dt);
    x_new = A * roe_d_mean + B * u;
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
