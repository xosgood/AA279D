% Paxton Scott and Chris Osgood
% 4/7/23
% AA279D HW1
clc; clear; close all;
addpath("Functions/Absolute_Orbits")

%% a) initial orbital elements
a = 7000; % km
e = 0.001;
i = deg2rad(98);
RAAN = 0;
omega = deg2rad(90);
nu = 0;
% oe_0 = [a; e; i; RAAN; omega; nu];

%% b) convert to ECI
[r_ECI_0, v_ECI_0] = OE2ECI(a, e, i, RAAN, omega, nu);

%% c) simulate orbit with state in ECI
mu = 3.986e5; % (km^3 / s^2)
T = 2 * pi * sqrt(a^3 / mu);
t_f = 10 * T;
tspan = [0, t_f]
state0 = [r_ECI_0; v_ECI_0];

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[t_out, y_out] = ode113(@func, tspan, state0, options);

% First, Plot the Earth
rE = 6378.1; % km
[xE, yE, zE] = ellipsoid(0, 0, 0, rE, rE, rE, 20);
figure
surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black'); 
axis equal;
view(3);
grid on;
hold on;
% now plot orbit
plot3(y_out(:,1), y_out(:,2), y_out(:,3), 'ro');





%% functions
function statedot = func(t, state)
    % State vector is [rx ry rz vx vy vz]’.
    % Although required by form, input value t will go unused here.
    mu = 3.986e5; % (km^3 / s^2)
    r = state(1:3);
    rdot = state(4:6);
    vdot = CentralBodyAccel(mu, r);
    statedot = [rdot;
                vdot];
end
