% Paxton Scott and Chris Osgood
% 4/21/23
% AA279D HW3 P1
clc; clear; close all;
addpath(genpath("Functions/"));

%% a) initial orbits
% constants
mu = 3.986e5;

% chief
a_0 = 7000; % km
e_0 = 0.001;
i_0 = deg2rad(98);
RAAN_0 = 0;
omega_0 = deg2rad(90);
nu_0 = 0;
oe_c = [a_0; e_0; i_0; RAAN_0; omega_0; nu_0];

% deputy
a_1 = a_0;
e_1 = 0.001;
i_1 = deg2rad(98.01);
RAAN_1 = deg2rad(0.05);
omega_1 = deg2rad(90.05);
nu_1 = deg2rad(0.01);
oe_d = [a_1; e_1; i_1; RAAN_1; omega_1; nu_1];

% sim parameters
n_orbits = 15;
n_steps_per_orbit = 30;
n_iter = n_steps_per_orbit * n_orbits;
T = 2 * pi * sqrt(a_0^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter);

%% b) convert initial conditions to useful coord frames
[r_ECI_0, v_ECI_0] = OE2ECI(a_0, e_0, i_0, RAAN_0, omega_0, nu_0);
[r_ECI_1, v_ECI_1] = OE2ECI(a_1, e_1, i_1, RAAN_1, omega_1, nu_1);

[r_RTN_d, v_RTN_d] = ECI2RTN(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1);

roe_singular = oe_d - oe_c;

% roe = [da, dlambda, dex, dey, dix, diy]^T
roe = OE2ROE(oe_c, oe_d);

v_ECI_c_expressedInRTN = rECI2RTN([r_ECI_0; v_ECI_0]) * v_ECI_0;

rho = norm(r_RTN_d);
r0 = norm(r_ECI_0);
r0_dot = v_ECI_c_expressedInRTN(1);

assert(rho / r0 <= 0.001);

theta = nu_0 + omega_0;
theta_dot = norm(cross(r_ECI_0, v_ECI_0))/r0^2; % v0 / r0

state0 = [r_RTN_d; theta; r0; v_RTN_d; theta_dot; r0_dot];


%% c) Calculate intergration constants. 

x_RTN_d = [r_RTN_d; v_RTN_d];
K = RTN2HCW_IC(x_RTN_d, a_0, 0);


%% d) Propagate. 

x_RTN_circular = zeros(n_iter, 6);

for iter = 1:size(x_RTN_circular, 1)
    x_RTN_circular(iter,:) = HCW_T2RTN(K, a_0, tspan(iter));
end

% Plot relative position in 3D. 
figure;
plot3(x_RTN_circular(:,1), x_RTN_circular(:,2), x_RTN_circular(:,3));
xlabel("R [km]");
ylabel("T [km]");
zlabel("N [km]");
title("Relative position in RTN over time.")

% Plot relative velocity in 3D. 
figure;
plot3(x_RTN_circular(:,4), x_RTN_circular(:,5), x_RTN_circular(:,6));
xlabel("R [km/s]");
ylabel("T [km/s]");
zlabel("N [km/s]");
title("Relative velocity in RTN over time. ")

% Plot relative position in TR, NR, TN plane.
figure;
subplot(3,1,1);
plot(x_RTN_circular(:,1), x_RTN_circular(:,2));
title("Relative position in RTN over time");
xlabel("R [km]");
ylabel("T [km]");
subplot(3,1,2);
plot(x_RTN_circular(:,1), x_RTN_circular(:,3));
xlabel("R [km]");
ylabel("N [km]");
subplot(3,1,3);
plot(x_RTN_circular(:,2), x_RTN_circular(:,3));
xlabel("T [km]");
ylabel("N [km]");

% Plot relative velocity in TR, NR, TN plane.
figure;
subplot(3,1,1);
plot(x_RTN_circular(:,4), x_RTN_circular(:,6));
title("Relative velocity in RTN over time");
xlabel("R [km/s]");
ylabel("T [km/s]");
subplot(3,1,2);
plot(x_RTN_circular(:,4), x_RTN_circular(:,6));
xlabel("R [km/s]");
ylabel("N [km/s]");
subplot(3,1,3);
plot(x_RTN_circular(:,5), x_RTN_circular(:,6));
xlabel("T [km/s]");
ylabel("N [km/s]");

%% e) Propagate with enforcing cyclical constraint on y_dot.
% Impose artificial constraint to ensure cyclic solution by linear solver.
v_RTN_d(2) = - 2*sqrt(mu/a_0^3)*r_RTN_d(1);

x_RTN_d_cyclic = [r_RTN_d; v_RTN_d];
K_cyclic = RTN2HCW_IC(x_RTN_d_cyclic, a_0, 0);

x_RTN_circular_cyclic = zeros(n_iter, 6);

for iter = 1:size(x_RTN_circular, 1)
    x_RTN_circular_cyclic(iter,:) = HCW_T2RTN(K_cyclic, a_0, tspan(iter));
end

PlotRTNSpace(x_RTN_circular_cyclic);
