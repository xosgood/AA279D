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
nu_c = deg2rad(-15);
M_c = TrueToMeanAnomaly(nu_c, e_c);
oe_c = [a_c; e_c; i_c; RAAN_c; omega_c; nu_c];

oe_c_osc = mean2osc(oe_c, 1);
[r_c_osc_0_ECI, v_c_osc_0_ECI] = OE2ECI(oe_c_osc);
x_c_osc_0 = [r_c_osc_0_ECI; v_c_osc_0_ECI];

% initial deputy orbit
roe_d = [0; 100; 0; 0; 0; 0] / a_c; % [m]
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
n_steps_per_orbit = 500;
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


%% 3) Least squares control solution for reconfiguration from in-train to passive safety ellipse 
u_burns = [0, pi/2, pi];
roe_i = roe_d';
roe_f = [0, 100, 0, 30, 0, 30] / a_c;

delta_vs = LS_control_solve(oe_c, roe_f - roe_i, u_burns);

num_reconfig_burns = length(u_burns);
reconfig_counter = 1;

%% 4) Formation keeping control parameters
adex_max = 10; % km
adlambda_max = 50; % km, plus/minus from desired dlambda
form_keep_man = false; % flag to keep track of if we have a formation keeping maneuver that we need to perform
form_keep_counter = 1;

%% 5) STM linear model for Quasi Nonsingular ROE with J2
QNS_roe_d_series_STM = zeros(6, n_iter);
oe_c_series = zeros(6, n_iter);
oe_c_mean_series = zeros(6, n_iter);
oe_d_mean_series = zeros(6, n_iter);
oe_d_series = zeros(6, n_iter);

% Set intial roe with the intial deputy roe. 
QNS_roe_d_series_STM(:, 1) = roe_d;
oe_d_mean_series(:,1) = oe_d_mean;

cumulative_delta_v = 0;
delta_v_series = zeros(n_iter, 1);
roes_desired = zeros(n_iter, 6);
roes_desired(1:n_steps_per_orbit, :) = ones(n_steps_per_orbit, 1) * roe_i;
roes_desired(n_steps_per_orbit:n_iter, :) = ones(n_iter-n_steps_per_orbit+1, 1) * roe_f;
r_d_RTN = zeros(3, n_iter);
v_d_RTN = zeros(3, n_iter);
r_d_ECI = zeros(3, n_iter);
v_d_ECI = zeros(3, n_iter);

dt = tspan(2) - tspan(1);
% propagate STM
for iter = 1:n_iter
    oe_c_series(:, iter) = ECI2OE(r_c_ECI(:, iter), v_c_ECI(:, iter))';
    oe_c_mean_series(:, iter) = osc2mean(oe_c_series(:, iter), 1);
    delta_v_series(iter) = cumulative_delta_v;
    
    % Create RTN vectors to visulize. 
    %QNS_roe_d_series_STM(2, iter) = 100/7000;
    oe_d(:, iter) = ROE2OE(oe_c_series(:,iter), QNS_roe_d_series_STM(:, iter));
    [r_d_ECI(:, iter), v_d_ECI(:, iter)] = OE2ECI(oe_d(:,iter));
    [r_d_RTN(:, iter), v_d_RTN(:, iter)] =  ECI2RTN(r_c_ECI(:, iter), v_c_ECI(:, iter), r_d_ECI(:, iter), v_d_ECI(:, iter));

    if iter < n_iter
        % this STM propagates mean roe's
        QNS_roe_d_series_STM(:, iter+1) = STM_QNS_ROE_J2(oe_c_mean_series(:,iter), QNS_roe_d_series_STM(:,iter), dt);
        oe_d_mean_series(:,iter+1) = ROE2OE(oe_c_mean_series(:,iter), QNS_roe_d_series_STM(:, iter+1));
    
        % Apply Reconfiguration Manuever
        % compute mean change in mean argument of latitidue for deputy 
        u_c_cur = oe_c_mean_series(5,iter) + oe_c_mean_series(6,iter);
        if (reconfig_counter <= num_reconfig_burns) && (abs(wrapToPi(u_c_cur) - wrapToPi(u_burns(reconfig_counter))) < 1e-2)
            QNS_roe_d_series_STM(:,iter+1) = ApplyDeputyManuever_NearCircular(...
                oe_c_mean_series(:,iter), QNS_roe_d_series_STM(:,iter+1), delta_vs(:,reconfig_counter));
            cumulative_delta_v = cumulative_delta_v + norm(delta_vs(:,reconfig_counter));
            reconfig_counter = reconfig_counter + 1;
        end
        
        % Compute Formation Keeping Control
        adex = a_c * QNS_roe_d_series_STM(3,iter+1);
        adlambda = a_c * QNS_roe_d_series_STM(2,iter+1);
        adlambda_desired = a_c * roe_f(2);
        adlambda_error = adlambda - adlambda_desired;
        if iter > 3 * n_steps_per_orbit ... % only perform formation keeping after we've gone through 1 orbit already (to account for reconfig maneuver)
            && (adex > adex_max || abs(adlambda_error) > adlambda_max) ... % dex or dlambda has gone beyond our dead-band threshold
            && ~form_keep_man % we aren't already in the middle of a formation keeping manuever
            % compute manuever for formation keeping
            u_fk_burns = [wrapToPi(u_c_cur), wrapToPi(u_c_cur + pi), wrapToPi(u_c_cur + 2*pi)];
            roe_i_fk = QNS_roe_d_series_STM(:,iter+1)';
            roe_f_fk = roe_f;
            roe_f_fk(3) = -adex_max / a_c; % push dex to other side of dead band
            if adex > adex_max
                roe_f_fk(4) = QNS_roe_d_series_STM(4,iter+1);
            end
            delta_vs_fk = LS_control_solve(oe_c_mean_series(:, iter), roe_f_fk - roe_i_fk, u_fk_burns);
            form_keep_man = true;
        end
        % Apply Formation Keeping Manuever
        if form_keep_man ...
            && (abs(wrapToPi(u_c_cur) - wrapToPi(u_fk_burns(form_keep_counter))) < 1e-2)
            QNS_roe_d_series_STM(:,iter+1) = ApplyDeputyManuever_NearCircular(...
                oe_c_mean_series(:,iter), QNS_roe_d_series_STM(:,iter+1), delta_vs_fk(:,form_keep_counter));
            cumulative_delta_v = cumulative_delta_v + norm(delta_vs(:,form_keep_counter));
            form_keep_counter = form_keep_counter + 1;
            if form_keep_counter > length(u_fk_burns)
                % then we have completed our formation keeping manuever
                % so we should reset our counter and flag
                form_keep_counter = 1;
                form_keep_man = false;
            end
        end
    end
end

figure(2);
PlotQNSROE_meters(QNS_roe_d_series_STM, a_c*1000);
subplot(3,1,1);
sgtitle("Mean relative orbital elements of deputy, with J2, STM");

figure(3);
plot((1:n_iter)/n_steps_per_orbit, delta_v_series);
grid on;
sgtitle("Delta V vs orbits");
xlabel("orbits");
ylabel("cumulative delta v")

% plot ROE vs time. 
figure(4);
hold on;
plot((1:n_iter)/n_steps_per_orbit, a_c * roes_desired(:,3));
plot((1:n_iter)/n_steps_per_orbit, a_c * QNS_roe_d_series_STM(3,:));
plot((1:n_iter)/n_steps_per_orbit, a_c * roes_desired(:,4));
plot((1:n_iter)/n_steps_per_orbit, a_c * QNS_roe_d_series_STM(4,:));
grid on;
hold off;
xlabel("orbits");
ylabel("ROEs");
legend("a \delta ex desired", "a \delta ex actual", "a \delta ey desired", "a \delta ey actual")
sgtitle("Desired vs. actual ROEs");

figure(5);
hold on;
plot((1:n_iter)/n_steps_per_orbit, a_c * roes_desired(:,5));
plot((1:n_iter)/n_steps_per_orbit, a_c * QNS_roe_d_series_STM(5,:));
plot((1:n_iter)/n_steps_per_orbit, a_c * roes_desired(:,6));
plot((1:n_iter)/n_steps_per_orbit, a_c * QNS_roe_d_series_STM(6,:));
grid on;
hold off;
xlabel("orbits");
ylabel("ROEs");
legend("a \delta ix desired", "a \delta ix actual", "a \delta iy desired", "a \delta iy actual")
sgtitle("Desired vs. actual ROEs");

figure(6);
hold on;
plot((1:n_iter)/n_steps_per_orbit, a_c * roes_desired(:,1));
plot((1:n_iter)/n_steps_per_orbit, a_c * QNS_roe_d_series_STM(1,:));
plot((1:n_iter)/n_steps_per_orbit, a_c * roes_desired(:,2));
plot((1:n_iter)/n_steps_per_orbit, a_c * QNS_roe_d_series_STM(2,:));
grid on;
hold off;
xlabel("orbits");
ylabel("ROEs");
legend("a \delta a desired", "a \delta a actual", "a \delta \lambda desired", "a \delta \lambda actual")
sgtitle("Desired vs. actual ROEs");

% Plot RTN
x_RTN = [r_d_RTN', v_d_RTN'];
PlotRTNSpace(x_RTN);

