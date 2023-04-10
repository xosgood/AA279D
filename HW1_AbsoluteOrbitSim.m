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
n_iter = 200;
n_orbits = 10;
mu = 3.986e5; % (km^3 / s^2)
T = 2 * pi * sqrt(a^3 / mu);
t_f = n_orbits * T;
tspan = linspace(0, t_f, n_iter) ; % [0, t_f];
state0 = [r_ECI_0; v_ECI_0];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, y_out] = ode113(@func, tspan, state0, options);
[t_out_J2, y_out_J2] = ode113(@func_J2, tspan, state0, options);

% First, Plot the Earth
rE = 6378.1; % km
[xE, yE, zE] = ellipsoid(0, 0, 0, rE, rE, rE, 20);
figure
surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black'); 
axis equal;
view(3);
grid on;
hold on;
% now plot orbits
plot3(y_out(:,1), y_out(:,2), y_out(:,3), 'r.');
plot3(y_out_J2(:,1), y_out_J2(:,2), y_out_J2(:,3), 'g.');
title("Orbit using ECI position and velocity as state");
legend("Earth", "Excluding J2", "Including J2");
xlabel('I (km)');
ylabel('J (km)');
zlabel('K (km)');

%% d) Keplerian propogation
% set times
t_0_MJD = mjuliandate(2023,07,04,12,0,0); % 07/04/2023 converted to MJD
sec_per_day = 86400;
t_f_MJD = t_0_MJD + n_orbits * T / sec_per_day;

% simulate
M_0 = 0; % initial mean anomaly
n_iter = length(t_out);
geod_station = [0; 0; 0]; % fake value since we don't care about this
kep_data = SimulateOrbitFromOE(a, e, i, RAAN, omega, M_0, geod_station, t_0_MJD, t_f_MJD, n_iter);
r_ECI_kep = kep_data.r_ECI_vec;
v_ECI_kep = kep_data.v_ECI_vec;

% plotting
% plot the Earth
rE = 6378.1; % km
[xE, yE, zE] = ellipsoid(0, 0, 0, rE, rE, rE, 20);
figure
surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black'); 
axis equal; view(3); grid on; hold on;
% plot orbits
plot3(r_ECI_kep(1,:), r_ECI_kep(2,:), r_ECI_kep(3,:), 'r.');
title("Orbit using Keplerian Propagation");
%legend("Earth", "Excluding J2", "Including J2");
xlabel('I (km)');
ylabel('J (km)');
zlabel('K (km)');

% error between keplerian sim and ECI sim
r_err_kep_ECI = abs(r_ECI_kep - y_out(:,1:3)');
v_err_kep_ECI = abs(v_ECI_kep - y_out(:,4:6)');
% plot errors
figure;
subplot(2, 1, 1);
plot(t_out', r_err_kep_ECI, t_out', vecnorm(r_err_kep_ECI));
title("Position Error in ECI Between Keplerian Propagation and ECI Propogation")
ylabel("km");
legend("I", "J", "K", "magnitude", "Location", "best");
subplot(2, 1, 2);
plot(t_out', v_err_kep_ECI, t_out', vecnorm(v_err_kep_ECI));
title("Velocity Error in ECI Between Keplerian Propagation and ECI Propogation")
ylabel("km/s");
xlabel("time (s)");


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

function statedot = func_J2(t, state)
    % State vector is [rx ry rz vx vy vz]’.
    % Although required by form, input value t will go unused here.
    mu = 3.986e5; % (km^3 / s^2) for earth
    R_E = 6378.1; % equatorial radius of earth in km
    J2 = 0.108263e-2; % for earth
    r = state(1:3);
    rdot = state(4:6);
    vdot = CentralBodyAccel(mu, r) + J2AccelECI(J2, mu, R_E, r(1), r(2), r(3));
    statedot = [rdot;
                vdot];
end
