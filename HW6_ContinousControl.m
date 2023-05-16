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
n_orbits = 30;
n_steps_per_orbit = 300;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_c^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);
dt = tspan(2) - tspan(1);
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
N = 4; % exponent in P matrix
k = 1; % (1/k) factor in front of P matrix
u_lowerbound = 1e-6; % lower bound on control actuation

%% relative motion propagator using STM (with J2)
roe_d_mean_series = zeros(6, n_iter);
roe_d_osc_series = zeros(6, n_iter);
oe_c_osc_series = zeros(6, n_iter);
oe_c_mean_series = zeros(6, n_iter);
oe_d_osc_series = zeros(6, n_iter);
oe_d_mean_series = zeros(6, n_iter);

% populate initial values
roe_d_mean_series(:,1) = roe_d_mean;
roe_d_osc_series(:,1) = roe_d_osc;%%%%%%%%%%%%
oe_c_osc_series(:,1) = oe_c_osc;
oe_c_mean_series(:,1) = oe_c_mean;
oe_d_osc_series(:,1) = oe_d_osc;
oe_d_mean_series(:,1) = oe_d_mean;

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
    Delta_roe_reducedmodel = [Delta_roe(1); Delta_roe(3:6)];
    roe_d_mean_reducedmodel = [roe_d_mean_series(1,iter); roe_d_mean_series(3:6,iter)];
    A = A_Control_ReducedModel(oe_c_mean_series(:,iter));
    B = B_Control_ReducedModel(oe_c_mean_series(:,iter));
    P = P_Control_ReducedModel(oe_c_mean_series(:,iter), Delta_roe_reducedmodel, N, k);
    u = - pinv(B) * (A * roe_d_mean_reducedmodel + P * Delta_roe_reducedmodel);
    
    % check lower bound for control effort
    if u(1) < u_lowerbound
        u(1) = 0;
    end
    if u(2) < u_lowerbound
        u(2) = 0;
    end
    
    delta_v_RTN = [0; u];
    
    % apply manuever
    roe_d_osc_series(:,iter) = ApplyDeputyManuever_NearCircular(...
        oe_c_osc_series(:,iter), roe_d_osc_series(:,iter), delta_v_RTN);
    
    % udpate roe's after doing manuever
    oe_d_osc_series(:,iter) = ROE2OE(oe_c_osc_series(:,iter), roe_d_osc_series(:,iter));
    oe_d_mean_series(:,iter) = osc2mean(oe_d_osc_series(:,iter), 1);
    roe_d_mean_series(:,iter) = OE2ROE(oe_c_mean_series(:,iter), oe_d_mean_series(:,iter));
    
end

