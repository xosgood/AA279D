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

n_orbits = 5;
n_steps_per_orbit = 100;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
dt = tspan(2) - tspan(1);
orbit_span = (1:n_iter)/n_steps_per_orbit;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

Q = (0.02 / a_c) * eye(n_dims_state);
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
Sigma_0 = 0.01 * eye(n_dims_state); % initial covariance estimate
pre_fit_res = zeros(m_dims_meas, n_iter); % pre fit residual
post_fit_res = zeros(m_dims_meas, n_iter); % post fit residual

Sigma(:,:,1) = Sigma_0;
mu(:,1) = mu_0;


for iter = 2:n_iter
    %%% ground truth (in ECI) propagation
    % apply control to ground truth
    dx_RTN = [0; 0; 0; u(:,iter-1)];
%     x_and_dv_ECI = RTN2ECI(x_absolute(7:12,iter-1), dx_RTN);
%     dv_ECI = x_and_dv_ECI(4:6);
%     x_absolute(10:12,iter) = x_absolute(10:12,iter-1) + dv_ECI;
    x_absolute(7:12,iter) = RTN2ECI(x_absolute(7:12,iter-1), dx_RTN);
    
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

    oe_c_osc_mes = ECI2OE(y(1:3,iter-1), y(4:6,iter-1));
    oe_d_osc_mes = ECI2OE(y(7:9,iter-1), y(10:12,iter-1));

    measured_roe(:,iter) = OE2ROE(oe_c_osc_mes, oe_d_osc_mes);

    
    
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

    pre_fit_res(:,iter) = (y(:,iter-1) - y_hat);

    mu(:,iter) = mu(:,iter) + K * (y(:,iter-1) - y_hat); % mean update
    Sigma(:,:,iter) = Sigma(:,:,iter) - K * Sigma_xy'; % covariance update

    post_fit_res(:,iter) = y(:,iter-1) - g(mu(:,iter), oe_c_osc_cur);
    
    %%% compute control for next timestep to use
    u(:,iter) = LyapunovController(roe_d_mean_cur, roe_desired, oe_c_mean_cur, lyap_params);
end


%% plotting
figure;
sgtitle("True and Estimated ROEs vs. time")
subplot(3,2,1); grid on; hold on;
plot(orbit_span, a_c * x_roe(1,:));
plot(orbit_span, a_c * mu(1,:));
plot(orbit_span(1:end), a_c * measured_roe(1,:));

xlabel("time [s]");
ylabel("a\delta a [km]")
subplot(3,2,2); grid on; hold on;
plot(orbit_span, a_c * x_roe(2,:));
plot(orbit_span, a_c * mu(2,:));
plot(orbit_span(1:end), a_c * measured_roe(2,:));
legend("ground truth", "state estimate", "measured", "Location","Best");
xlabel("time [s]");
ylabel("a\delta \lambda [km]")
subplot(3,2,3); grid on; hold on;
plot(orbit_span, a_c * x_roe(3,:));
plot(orbit_span, a_c * mu(3,:));
plot(orbit_span(1:end), a_c * measured_roe(3,:));
xlabel("time [s]");
ylabel("a\delta e_x [km]")
subplot(3,2,4); grid on; hold on;
plot(orbit_span, a_c * x_roe(4,:));
plot(orbit_span, a_c * mu(4,:));
plot(orbit_span(1:end), a_c * measured_roe(4,:));
xlabel("time [s]");
ylabel("a\delta e_y [km]")
subplot(3,2,5); grid on; hold on;
plot(orbit_span, a_c * x_roe(5,:));
plot(orbit_span, a_c * mu(5,:));
plot(orbit_span(1:end), a_c * measured_roe(5,:));
xlabel("time [s]");
ylabel("a\delta i_x [km]")
subplot(3,2,6); grid on; hold on;
plot(orbit_span, a_c * x_roe(6,:));
plot(orbit_span, a_c * mu(6,:));
plot(orbit_span(1:end), a_c * measured_roe(6,:));
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

position_confidence_interval = 1.96 * reshape(sqrt(Sigma(1,1,:) + Sigma(2,2,:) + Sigma(3,3,:)),[1,n_iter]);
velocity_confidence_interval = 1.96 * reshape(sqrt(Sigma(4,4,:) + Sigma(5,5,:) + Sigma(6,6,:)), [1,n_iter]);

figure; 
sgtitle("Pre and Post fit residuals verus time")
subplot(2, 1, 1); hold on; grid on;
plot(orbit_span, vecnorm(pre_fit_res(7:9,:)), ".");
plot(orbit_span, vecnorm(post_fit_res(7:9,:)), ".red");
legend("Pre fit residual", "Post fit residual");
xlabel("orbits")
ylabel("norm residual [km]")
ylim([-5, 15]);

subplot(2, 1, 2); hold on; grid on;
plot(orbit_span, vecnorm(pre_fit_res(10:12,:)), ".");
plot(orbit_span, vecnorm(post_fit_res(10:12,:)), ".red");

legend("Pre fit residual", "Post fit residual");
ylim([-5, 15]);
xlabel("orbits")
ylabel("norm residual [km/s]")
ylim([-0.01, 0.05]);

% plot state estimate error with covariance
alpha = 1.96;
figure;
sgtitle("UKF absolute value error with 95% upper confidence bound");
subplot(3,2,1); grid on; hold on;
plot(orbit_span, a_c * err(1,:));
plot(orbit_span, alpha * a_c * squeeze(sqrt(Sigma(1,1,:))));
xlabel("orbits");
ylabel("a \delta a error [km]");
ylim([0, 10]);
subplot(3,2,2); grid on; hold on;
plot(orbit_span, a_c * err(2,:));
plot(orbit_span, alpha * a_c * squeeze(sqrt(Sigma(2,2,:))));
xlabel("orbits");
ylabel("a \delta \lambda error [km]");
ylim([0, 10]);
legend("error", "upper 95% confidence bound");
subplot(3,2,3); grid on; hold on;
plot(orbit_span, a_c * err(3,:));
plot(orbit_span, alpha * a_c * squeeze(sqrt(Sigma(3,3,:))));
xlabel("orbits");
ylabel("a \delta e_x error [km]");
ylim([0, 10]);
subplot(3,2,4); grid on; hold on;
plot(orbit_span, a_c * err(4,:));
plot(orbit_span, alpha * a_c * squeeze(sqrt(Sigma(4,4,:))));
xlabel("orbits");
ylabel("a \delta e_y error [km]");
ylim([0, 10]);
subplot(3,2,5); grid on; hold on;
plot(orbit_span, a_c * err(5,:));
plot(orbit_span, alpha * a_c * squeeze(sqrt(Sigma(5,5,:))));
xlabel("orbits");
ylabel("a \delta i_x error [km]");
ylim([0, 5]);
subplot(3,2,6); grid on; hold on;
plot(orbit_span, a_c * err(6,:));
plot(orbit_span, alpha * a_c * squeeze(sqrt(Sigma(6,6,:))));
xlabel("orbits");
ylabel("a \delta i_y error [km]");
ylim([0, 5]);

% plot control input history
figure;
sgtitle("Delta-v components vs number of orbits passed");
subplot(4,1,1);
plot(orbit_span, 1000 * u(2,:));
grid on;
xlabel("orbits");
ylabel("tangential delta-v (m/s)");
subplot(4,1,2);
plot(orbit_span, 1000 * u(3,:));
grid on;
xlabel("orbits");
ylabel("normal delta-v (m/s)");
subplot(4,1,3);
plot(orbit_span, 1000 * vecnorm(u));
grid on;
xlabel("orbits");
ylabel("magnitude of delta-v (m/s)");
subplot(4,1,4);
plot(orbit_span, 1000 * cumsum(vecnorm(u), 2));
grid on;
xlabel("orbits");
ylabel("cumulative delta-v (m/s)");

% plot control tracking error
control_tracking_error = x_roe - roe_desired;
figure;
sgtitle("Control tracking error in ROE space");
subplot(3,2,1); grid on; hold on;
plot(orbit_span, a_c * control_tracking_error(1,:));
xlabel("orbits");
ylabel("a \delta a [km]");
subplot(3,2,2); grid on; hold on;
plot(orbit_span, a_c * control_tracking_error(2,:));
xlabel("orbits");
ylabel("a \delta \lambda control tracking error [km]");
subplot(3,2,3); grid on; hold on;
plot(orbit_span, a_c * control_tracking_error(3,:));
xlabel("orbits");
ylabel("a \delta e_x control tracking error [km]");
subplot(3,2,4); grid on; hold on;
plot(orbit_span, a_c * control_tracking_error(4,:));
xlabel("orbits");
ylabel("a \delta e_y control tracking error [km]");
subplot(3,2,5); grid on; hold on;
plot(orbit_span, a_c * control_tracking_error(5,:));
xlabel("orbits");
ylabel("a \delta i_x control tracking error [km]");
subplot(3,2,6); grid on; hold on;
plot(orbit_span, a_c * control_tracking_error(6,:));
xlabel("orbits");
ylabel("a \delta i_y control tracking error [km]");

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

function GraphErrorTerms(tspan, x_true, x_pred, Sigma, n)
    figure; grid on; hold on;
    title("Estimation error vs. time")

    err = abs(x_true(1,:)- x_pred(1,:));
    plot(tspan, err,'red');
    plot(tspan, err - reshape(Sigma(1,1,:),[1,n]), 'blue');
    plot(tspan, err + reshape(Sigma(1,1,:),[1,n]), 'blue');
    %xconf = [tspan tspan(end:-1:1)];
    %yconf = [err - reshape(Sigma(1,1,:),[1,n]), err + reshape(Sigma(1,1,:),[1,n])];
    %p = fill(xconf,yconf, 'red');
    %p.FaceColor = [1 0.8 0.8];      
    %p.EdgeColor = 'none';  
   
    legend("\delta a error", "upper and low confidence bound of 1 variance")
    %legend("\delta a error", "\delta \Omega error", "\delta i_x error", "\delta i_y error", "\delta e_x error", "\delta e_y error");
    xlabel("time [s]")
    ylabel("\delta OE")
end

