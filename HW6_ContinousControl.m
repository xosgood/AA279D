% Chris Osgood
% 5/15/23
% AA279D HW5
% Continous Control
clc; clear; close all;
addpath(genpath("Functions/"));

%% constants
mu = 3.986e5;

%% initialize chief
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

%% initialize deputy
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

%% simulation parameters
n_orbits = 300;
n_steps_per_orbit = 100;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
dt = tspan(2) - tspan(1);
orbit_span = (1:n_iter)/n_steps_per_orbit;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

%% chief absolute orbit full nonlinear propagator 
[~, x_c_ECI] = ode113(@AbsoluteOrbitWithJ2DiffEq, tspan, x_c_osc_0, options);
x_c_ECI = x_c_ECI';
r_c_ECI = x_c_ECI(1:3,:);
v_c_ECI = x_c_ECI(4:6,:);

%% deputy absolute orbit full nonlinear propagator (for ground truth)
[~, x_d_ECI] = ode113(@AbsoluteOrbitWithJ2DiffEq, tspan, x_d_osc_0, options);
x_d_ECI = x_d_ECI';
r_d_ECI = x_d_ECI(1:3,:);
v_d_ECI = x_d_ECI(4:6,:);

%% reconfiguration manuever setup
roe_desired = [0; 100; 0; 30; 0; 30] / a_c;


%% Lyapunov controller setup
N = 10; % exponent in P matrix
k = 1; % (1/k) factor in front of P matrix
u_lowerbound = 1e-6; % lower bound on control actuation
u_upperbound = 1e-3; % upper bound on control actuation
dlambda_thresh = 0.1;
dlambda_dot = 0.001;

%% relative motion propagator using STM (with J2)
roe_d_mean_series = zeros(6, n_iter);
roe_d_osc_series = zeros(6, n_iter);
oe_c_osc_series = zeros(6, n_iter);
oe_c_mean_series = zeros(6, n_iter);
oe_d_osc_series = zeros(6, n_iter);
oe_d_mean_series = zeros(6, n_iter);

% populate initial values
roe_d_mean_series(:,1) = roe_d_mean;
roe_d_osc_series(:,1) = roe_d_osc;
oe_c_osc_series(:,1) = oe_c_osc;
oe_c_mean_series(:,1) = oe_c_mean;
oe_d_osc_series(:,1) = oe_d_osc;
oe_d_mean_series(:,1) = oe_d_mean;

% delta v for plotting
u_series = zeros(2, n_iter);
cumulative_delta_v = 0;
delta_v_cum_series = zeros(n_iter, 1);

% delta roe series for plotting 
delta_roe_series = zeros(6, n_iter);

% rtn for plotting. 
r_d_RTN = zeros(3, n_iter);
v_d_RTN = zeros(3, n_iter);
r_d_ECI_with_control = zeros(3, n_iter);
v_d_ECI_with_control = zeros(3, n_iter);

for iter = 2:n_iter
    % chief ECI from ground truth to OE
    oe_c_osc_series(:,iter) = ECI2OE(r_c_ECI(:, iter), v_c_ECI(:, iter))';
    oe_c_mean_series(:,iter) = osc2mean(oe_c_osc_series(:, iter), 1);
    
    % propagate deputy
    roe_d_mean_series(:,iter) = STM_QNS_ROE_J2(oe_c_mean_series(:,iter), roe_d_mean_series(:,iter-1), dt);
    oe_d_mean_series(:,iter) = ROE2OE(oe_c_mean_series(:,iter), roe_d_mean_series(:, iter));
    oe_d_osc_series(:,iter) = mean2osc(oe_d_mean_series(:,iter), 1);
    roe_d_osc_series(:,iter) = OE2ROE(oe_c_osc_series(:,iter), oe_d_osc_series(:,iter));
    
    % compute Lyapunov control
    Delta_roe = roe_d_mean_series(:,iter) - roe_desired;
    if abs(Delta_roe(2)) >= dlambda_thresh
        % impose a desired da to cause a desired drift in dlambda
        da_des = -2/3 * sign(-Delta_roe(2)) * dlambda_dot / n_c;
        Delta_roe(1) = roe_d_mean_series(1,iter) - da_des;
    end
    Delta_roe_reducedmodel = [Delta_roe(1); Delta_roe(3:6)];
    roe_d_mean_reducedmodel = [roe_d_mean_series(1,iter); roe_d_mean_series(3:6,iter)];
    
    A = A_Control_ReducedModel(oe_c_mean_series(:,iter));
    B = B_Control_ReducedModel(oe_c_mean_series(:,iter));
    P = P_Control_ReducedModel(oe_c_mean_series(:,iter), Delta_roe_reducedmodel, N, k);
    u = - pinv(B) * (A * roe_d_mean_reducedmodel + P * Delta_roe_reducedmodel);
    
    % check lower bound for control effort
    if abs(u(1)) < u_lowerbound
        u(1) = 0;
    end
    if abs(u(2)) < u_lowerbound
        u(2) = 0;
    end
    % check upper bound for control effort
    if abs(u(1)) > u_upperbound
        u(1) = sign(u(1)) * u_upperbound;
    end
    if abs(u(2)) > u_upperbound
        u(2) = sign(u(2)) * u_upperbound;
    end
    
    u_series(:,iter) = u;
    
    delta_v_RTN = [0; u];
    
    % apply manuever
    roe_d_osc_series(:,iter) = ApplyDeputyManuever_NearCircular(...
        oe_c_osc_series(:,iter), roe_d_osc_series(:,iter), delta_v_RTN);
    
    % udpate roe's after doing manuever
    oe_d_osc_series(:,iter) = ROE2OE(oe_c_osc_series(:,iter), roe_d_osc_series(:,iter));
    oe_d_mean_series(:,iter) = osc2mean(oe_d_osc_series(:,iter), 1);
    roe_d_mean_series(:,iter) = OE2ROE(oe_c_mean_series(:,iter), oe_d_mean_series(:,iter));
        
    % total delta-v
    dv_cur = norm(delta_v_RTN);
    cumulative_delta_v = cumulative_delta_v + dv_cur;
    delta_v_cum_series(iter) = cumulative_delta_v;

    % RTN for plotting
    [r_d_ECI_with_control(:, iter), v_d_ECI_with_control(:, iter)] = OE2ECI(oe_d_mean_series(:,iter));
    [r_d_RTN(:, iter), v_d_RTN(:, iter)] =  ECI2RTN(r_c_ECI(:, iter), v_c_ECI(:, iter), r_d_ECI_with_control(:, iter), v_d_ECI_with_control(:, iter));

    % Delta ROE for plotting 
    delta_roe_series(:,iter) = Delta_roe;
    
end


figure(2);
PlotQNSROE_meters(roe_d_mean_series, a_c*1000);
subplot(3,1,1);
sgtitle("Mean relative orbital elements of deputy, with J2, STM");

figure(3);
plot(orbit_span, delta_v_cum_series);
title("Cumulative delta-v vs number of orbits passed");
xlabel("number of orbits");
ylabel("delta-v (km/s)");

figure(4);
plot(orbit_span, u_series(1,:), orbit_span, u_series(2,:));
legend("tangential delta v", "normal delta v")
xlabel("Orbits")
ylabel("\Delta v")

figure(8)
hold on;
plot(orbit_span, delta_roe_series(1,:));
plot(orbit_span, delta_roe_series(2,:));
plot(orbit_span, delta_roe_series(3,:));
plot(orbit_span, delta_roe_series(4,:));
plot(orbit_span, delta_roe_series(5,:));
plot(orbit_span, delta_roe_series(6,:));
hold off;
xlabel("orbits");
ylabel("ROE Error")
legend("\Delta \delta a", "\Delta \delta \lambda", "\Delta \delta e_x", "\Delta \delta e_y", "\Delta \delta i_x", "\Delta \delta i_y")


% Plot RTN
x_RTN = [r_d_RTN', v_d_RTN'];
PlotRTNSpace(x_RTN);

