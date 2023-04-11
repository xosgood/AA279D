function [x_vec, y_vec, t_vec, nu_vec, r_vec] = SimulatePerifocalOrbit_WithJ2(a, e, i, t_0, t_f, M_0, n_iter, eps)
    % SIMULATEPERIFOCALORBIT_WITHJ2 simulates a closed orbit in perifocal coordinates
    % from initial time t_0 to final time t_f. It computes n_iter points
    % starting from M_0. 
    % Returns the x, y coordinates in the PQR perifocal coordinate frame. 
    % Assumes central body is the Earth, and it includes J2's mean affect on the mean anomaly.
    % 'a' must be provided in km. The times must be provided in seconds,
    % and all angles must be in radians.
    mu = 3.986e5; % of Earth (km^3 / s^2)
    R_E = 6378.1; % equatorial radius of Earth in km
    J2 = 0.108263e-2; % of Earth
    n = sqrt(mu / a^3);
    p = a * (1 - e^2);
    % compute mean J2 affect on mean anomaly (this formula includes the usual affect from mean motion)
    dMdt = n + (3/4) * n * J2 * (R_E / (a * (1-e^2)))^2 * sqrt(1-e^2) * (3*cos(i)^2 - 1);
    t_vec = linspace(t_0, t_f, n_iter);
    nu_vec = zeros(1,n_iter);
    r_vec = zeros(1,n_iter);
    for iter = 1:n_iter
        M = M_0 + dMdt * (t_vec(iter) - t_0);
        E = MeanToEccentricAnomaly(M, e, eps);
        nu_vec(iter) = EccentricToTrueAnomaly(E, e);
        r_vec(iter) = p / (1 + e * cos(nu_vec(iter)));
    end
    [x_vec, y_vec] = pol2cart(nu_vec, r_vec);
end

