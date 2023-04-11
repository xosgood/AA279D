% Paxton Scott and Chris Osgood
% 4/7/23
% AA279D HW1
clc; clear; close all;
addpath("Functions/")
addpath("Functions/Absolute_Orbits/")

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
n_orbits = 100;
n_iter = 30 * n_orbits;
n_steps = n_iter * n_orbits;
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
figure;
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
r_RTN_kep = kep_data.r_RTN_vec;
v_RTN_kep = kep_data.v_RTN_vec;

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

% change y_out from ECI to RTN
y_out_RTN = zeros(size(y_out));
for i = 1:length(y_out)
    r_ECI_cur = y_out(i,1:3)';
    v_ECI_cur = y_out(i,4:6)';
    R_ECI2RTN_cur = rECI2RTN(y_out(i,:)');
    y_out_RTN(i,1:3) = R_ECI2RTN_cur * r_ECI_cur;
    y_out_RTN(i,4:6) = R_ECI2RTN_cur * v_ECI_cur;
end

% error between keplerian sim and ECI sim
r_err_RTN_kep = abs(r_RTN_kep - y_out_RTN(:,1:3)');
v_err_RTN_kep = abs(v_RTN_kep - y_out_RTN(:,4:6)');
% plot errors
figure;
subplot(2, 1, 1);
plot(t_out', r_err_RTN_kep, t_out', vecnorm(r_err_RTN_kep));
title("Position Error in ECI Between Keplerian Propagation and ECI Propogation")
ylabel("position error [km]");
legend("R", "T", "N", "magnitude", "Location", "best");
subplot(2, 1, 2);
plot(t_out', v_err_RTN_kep, t_out', vecnorm(v_err_RTN_kep));
title("Velocity Error in ECI Between Keplerian Propagation and ECI Propogation")
ylabel("velocity error [km/s]");
xlabel("time [s]");

%% e) Compute and plot Keplerian element. 
PlotOrbitalElements(y_out, mu, "Orbital elements")
PlotOrbitalElements(y_out_J2, mu, "Orbital elements with J2 effect")



%% functions
function PlotOrbitalElements(y, mu, title_string)
    r_eci = y(:, 1:3);
    v_eci = y(:,4:6);
    
    % osculating orbital elements. 
    oe = zeros(length(y), 6);
    
    for i= 1:length(r_eci)
       [a, e, inc, RAAN, omega, nu] = ECI2OE(r_eci(i,:), v_eci(i,:));
       oe(i,:) = [a, e, inc, RAAN, omega, nu];
    end
    
    figure
    sgtitle(title_string)
    subplot(4, 1, 1);
    plot(1:length(y), oe)
    title("Orbital elemens vs. time step")
    legend("a", "e", "i", "RAAN", "Omega", "Nu")
    
    % Calculate angular momentum vector.
    h = zeros(size(r_eci));
    
    for i= 1:length(r_eci)
       h(i,:) = cross(r_eci(i,:), v_eci(i,:));
    end
    
    subplot(4, 1, 2);
    plot(1:length(y), h)
    title("Angular momentum components vs. time steps")
    legend("X", "Y", "Z")
    
    % Calculate eccentricity vector: https://en.wikipedia.org/wiki/Eccentricity_vector. 
    e = EccentricityVector(mu, r_eci, v_eci);
    
    subplot(4, 1, 3);
    plot(1:length(y), e)
    title("Eccentricity vector components vs. time steps")  
    
    % Calculate specific mechanical energy.
    mechanical_energy = zeros(length(y));
    
    for i = 1:length(y)
        mechanical_energy(i) = 0.5 *dot(v_eci(i,:), v_eci(i,:)) + mu/norm(r_eci(i,:));
    end
    
    subplot(4, 1, 4);
    plot(1:length(y), mechanical_energy)
    title("Mechanical energy vs. time steps")
    
end

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