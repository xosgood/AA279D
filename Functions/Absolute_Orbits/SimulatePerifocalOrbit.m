function [x_vec, y_vec, t_vec, nu_vec, r_vec] = SimulatePerifocalOrbit(a, e, t_0, t_f, M_0, n_iter, eps)
    % SIMULATEPERIFOCALORBIT simulates a closed orbit in perifocal coordinates
    % from initial time t_0 to final time t_f. It computes n_iter points
    % starting from M_0. 
    % Returns the x, y coordinates in the PQR perifocal coordinate frame. 
    % Assumes central body is the Earth.
    % 'a' must be provided in km. The times must be provided in seconds,
    % and all angles must be in radians.
    mu = 3.986e5; % of Earth (km^3 / s^2)
    n = sqrt(mu / a^3);
    p = a * (1 - e^2);
    t_vec = linspace(t_0, t_f, n_iter);
    nu_vec = zeros(1,n_iter);
    r_vec = zeros(1,n_iter);
    for i = 1:n_iter
        M = M_0 + n * (t_vec(i) - t_0);
        E = MeanToEccentricAnomaly(M, e, eps);
        nu_vec(i) = EccentricToTrueAnomaly(E, e);
        r_vec(i) = p / (1 + e * cos(nu_vec(i)));
    end
    [x_vec, y_vec] = pol2cart(nu_vec, r_vec);
end

