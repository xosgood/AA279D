function [a, e, i, RAAN, omega, nu] = ECI2OE(r_ECI, v_ECI)
    % ECI2OE converts position and velocity vectors in the ECI frame to the
    % 6 keplerian orbital elements. Note that this assumes the central body
    % is the Earth.
    % As of now, this function cannot handle circular/equatorial orbits!
    % Inputs must be in km and km/s respectively. 
    % Output 'a' will be in km, and all angles will be in radians. 
    mu = 3.986e5; % of Earth (km^3 / s^2)
    h = cross(r_ECI, v_ECI);
    W = h / norm(h);
    i = atan2(sqrt(W(1)^2 + W(2)^2), W(3));
    RAAN = atan2(W(1), -W(2));
    p = norm(h)^2 / mu;
    a = 1 / ( (2/norm(r_ECI)) - (norm(v_ECI)^2/mu) );
    n = sqrt(mu / a^3);
    e = sqrt(1 - p/a);
    E = atan2( dot(r_ECI,v_ECI) / (a^2 * n), 1 - norm(r_ECI)/a);
    nu = EccentricToTrueAnomaly(E, e);
    u = atan2( r_ECI(3) / sin(i), r_ECI(1) * cos(RAAN) + r_ECI(2) * sin(RAAN) );
    omega = wrapTo2Pi(u - nu);
    RAAN = wrapTo2Pi(RAAN);
end
