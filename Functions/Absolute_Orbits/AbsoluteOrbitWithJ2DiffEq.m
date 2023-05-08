function statedot = AbsoluteOrbitWithJ2DiffEq(t, state)
    % ECI absolute orbit ODE, with J2 perturbations
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