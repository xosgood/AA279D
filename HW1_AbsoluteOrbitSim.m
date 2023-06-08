% Paxton Scott and Chris Osgood
% 4/7/23
% AA279D HW1
clc; clear; close all;
addpath(genpath("Functions/"))

%% a) initial orbital elements
a_c = 7000; % km
e_c = 0.001;
i_c = deg2rad(98);
RAAN_c = 0;
omega_c = deg2rad(90);
nu_c = 0;

oe_c = [a_c; e_c; i_c; RAAN_c; omega_c; nu_c];

oe_c_osc = mean2osc(oe_c, true);

%% b) convert to ECI
[r_ECI_0, v_ECI_0] = OE2ECI(oe_c_osc);

%% c) simulate orbit with state in ECI
n_orbits = 30;
n_iter = 60 * n_orbits;
mu = 3.986e5; % (km^3 / s^2)
T = 2 * pi * sqrt(a_c^3 / mu);
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
n_iter = length(t_out);
geod_station = [0; 0; 0]; % fake value since we don't care about this
kep_data = SimulateOrbitFromOE(oe_c, geod_station, t_0_MJD, t_f_MJD, n_iter);
r_ECI_kep = kep_data.r_ECI_vec;
v_ECI_kep = kep_data.v_ECI_vec;

r_RTN_kep = zeros(size(r_ECI_kep));
v_RTN_kep = zeros(size(r_ECI_kep));

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
for iter = 1:length(y_out)
    r_ECI_cur = y_out(iter,1:3)';
    v_ECI_cur = y_out(iter,4:6)';
    R_ECI2RTN_cur = rECI2RTN(y_out(iter,:)');
    y_out_RTN(iter,1:3) = R_ECI2RTN_cur * r_ECI_cur;
    y_out_RTN(iter,4:6) = R_ECI2RTN_cur * v_ECI_cur;
    
    R_ECI2RTN_kep = rECI2RTN([r_ECI_kep(:,iter); v_ECI_kep(:,iter)]);
    r_RTN_kep(:,iter) = R_ECI2RTN_kep * r_ECI_kep(:,iter);
    v_RTN_kep(:,iter) = R_ECI2RTN_kep * v_ECI_kep(:,iter);
end

% error between keplerian sim and ECI sim, excluding J2 effects.
r_err_RTN_kep = abs(r_RTN_kep - y_out_RTN(:,1:3)');
v_err_RTN_kep = abs(v_RTN_kep - y_out_RTN(:,4:6)');
% plot errors
figure;
subplot(2, 1, 1);
plot(t_out', r_err_RTN_kep, t_out', vecnorm(r_err_RTN_kep));
title("Position Error in RTN Between Keplerian Propagation and ECI Propogation")
ylabel("position error [km]");
legend("R", "T", "N", "magnitude", "Location", "best");
subplot(2, 1, 2);
plot(t_out', v_err_RTN_kep, t_out', vecnorm(v_err_RTN_kep));
title("Velocity Error in RTN Between Keplerian Propagation and ECI Propogation")
ylabel("velocity error [km/s]");
legend("R", "T", "N", "magnitude", "Location", "best");
xlabel("time [s]");

%% e) Compute and plot Keplerian element. 
figure;
plot(t_out', vecnorm(y_out(:,1:3), 2, 2));

PlotOrbitalElements(y_out, mu, "Without J2 effect", t_out);
PlotOrbitalElements(y_out_J2, mu, "With J2 effect", t_out_J2);

%% f) Apply linear mean J2 effects to analytical Keplerian propagation and compare with ECI sim with J2 effects
kep_data_j2 = SimulateOrbitFromOE_WithJ2(oe_c(1), oe_c(2), oe_c(3), oe_c(4), oe_c(5), oe_c(6), geod_station, t_0_MJD, t_f_MJD, n_iter);
r_ECI_kep_J2 = kep_data_j2.r_ECI_vec;
v_ECI_kep_J2 = kep_data_j2.v_ECI_vec;

% plotting
% plot the Earth
rE = 6378.1; % km
[xE, yE, zE] = ellipsoid(0, 0, 0, rE, rE, rE, 20);
figure;
surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black'); 
axis equal; view(3); grid on; hold on;
% plot orbits
plot3(r_ECI_kep_J2(1,:), r_ECI_kep_J2(2,:), r_ECI_kep_J2(3,:), 'r-');
title("Orbit using Keplerian Propagation, Including J2");
xlabel('I (km)');
ylabel('J (km)');
zlabel('K (km)');

% get oe at every timestep for keplerian sim with J2 effect and ECI sim with J2
oe_eciJ2_vec = zeros(6, length(r_ECI_kep_J2));
oe_kepJ2_vec = zeros(6, length(r_ECI_kep_J2));
for iter = 1:length(r_ECI_kep_J2)
    oe_J2 = ECI2OE(y_out_J2(iter,1:3), y_out_J2(iter,4:6));
    oe_kep_J2 = ECI2OE(r_ECI_kep_J2(:,iter), v_ECI_kep_J2(:,iter));

    oe_eciJ2_vec(:,iter) = [oe_J2(1), oe_J2(2), oe_J2(3), wrapToPi(oe_J2(4)), oe_J2(5), oe_J2(6)]';
    oe_kepJ2_vec(:,iter) = [oe_kep_J2(1), oe_kep_J2(2), oe_kep_J2(3), wrapToPi(oe_kep_J2(4)), oe_kep_J2(5), oe_kep_J2(6)]';
end

figure;
sgtitle("Comparing Kep sim and ECI sim, both with J2");
subplot(3,1,1);
plot(t_out, oe_eciJ2_vec(1,:));
hold on;
plot(t_out, oe_kepJ2_vec(1,:));
title("a vs. time");
legend("ECI", " Kep");
ylabel("a [km]");
subplot(3,1,2);
plot(t_out, oe_eciJ2_vec(2:5,:));
hold on;
plot(t_out, oe_kepJ2_vec(2:5,:));
title("e, i, RAAN, omega vs. time");
legend("e, ECI", "i, ECI", "RAAN, ECI", "omega, ECI", ...
       "e, Kep", "i, Kep", "RAAN, Kep", "omega, Kep");
ylabel(["e is unitless;","all angles in radians"]);
subplot(3,1,3);

plot(t_out, oe_eciJ2_vec(6,:));
hold on;
plot(t_out, oe_kepJ2_vec(6,:));
title("nu vs. time");
ylabel("nu [rad]");
xlabel("time [s]");
legend("ECI", " Kep");

% plot errors between oe with J2 from kep and from ECI sims
figure;
sgtitle("Errors: Comparing Kep sim and ECI sim, both with J2");
subplot(4,1,1);
plot(t_out, abs(oe_eciJ2_vec(1,:)-oe_kepJ2_vec(1,:)));
title("a error vs. time");
ylabel("a [km]");
subplot(4,1,2);
plot(t_out, abs(oe_eciJ2_vec(2:4,:)-oe_kepJ2_vec(2:4,:)));
title("e, i, RAAN errors vs. time");
legend("e error", "i error", "RAAN error");
ylabel(["e is unitless;","all angles in radians"]);
subplot(4,1,3);
plot(t_out, abs(oe_eciJ2_vec(5,:)-oe_kepJ2_vec(5,:)));
title("omega error vs. time");
ylabel("omega [rad]");
xlabel("time [s]");
subplot(4,1,4);
plot(t_out, abs(oe_eciJ2_vec(6,:)-oe_kepJ2_vec(6,:)));
title("nu error vs. time");
ylabel("nu [rad]");
xlabel("time [s]");


%% functions
function PlotOrbitalElements(y, mu, title_string, t)
    r_eci = y(:,1:3);
    v_eci = y(:,4:6);
    
    % osculating orbital elements. 
    oe = zeros(length(y), 6);
    
    for i = 1:length(r_eci)
       oe(i,:) = ECI2OE(r_eci(i,:), v_eci(i,:));
    end
    
    figure
    sgtitle(title_string)
    subplot(3, 2, 1)
    plot(t, oe(:,1))
    title("a vs time.")
    ylabel("a [km]")
    xlabel("time [s]")

    subplot(3, 2, 2);
    plot(t, oe(:,2))
    title("Eccentricty vs. time.")
    ylabel("[radians]")
    xlabel("time [s]")

    subplot(3, 2, 3);
    plot(t, oe(:,3))
    title("Inlination vs. time.")
    ylabel("[radians]")
    xlabel("time [s]")

    subplot(3, 2, 4);
    plot(t, oe(:,4))
    title("RAAN vs. time.")
    ylabel("[radians]")
    xlabel("time [s]")

    subplot(3, 2, 5);
    plot(t, oe(:,5))
    title("\omega vs. time.")
    ylabel("[radians]")
    xlabel("time [s]")

    subplot(3, 2, 6);
    plot(t, oe(:,6))
    title("\nu vs. time.")
    ylabel("[radians]")
    xlabel("time [s]")
    
    % Calculate angular momentum vector.
    h = zeros(size(r_eci));
    
    for i = 1:length(r_eci)
       h(i,:) = cross(r_eci(i,:), v_eci(i,:));
    end
    
    figure
    sgtitle(title_string)
    subplot(3, 1, 1);
    plot(t, h)
    title("Angular momentum components vs. time.")
    legend("X", "Y", "Z")
    ylabel("[km*kg/s]")
    
    
    % Calculate eccentricity vector: https://en.wikipedia.org/wiki/Eccentricity_vector. 
    e = zeros(size(r_eci));
    for i= 1:length(r_eci)
        e(i,:) = EccentricityVector(mu, r_eci(i,:), v_eci(i,:));
    end

    subplot(3, 1, 2);
    plot(t, e)
    title("Eccentricity vector components vs. time.")
    legend("X", "Y", "Z")
    ylabel("[km]")
    
    % Calculate specific mechanical energy.
    mechanical_energy = zeros(1, length(y));
    
    for i = 1:length(y)
        mechanical_energy(i) = 0.5 *dot(v_eci(i,:), v_eci(i,:)) + mu/norm(r_eci(i,:));
    end
    
    subplot(3, 1, 3);
    plot(t, mechanical_energy)
    title("Mechanical energy vs. time.")
    ylabel("[km^2*kg/s^2]")
    xlabel("time [s]")
    
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
