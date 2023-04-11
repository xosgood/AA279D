function data_sim = SimulateOrbitFromOE_WithJ2(a, e, i, RAAN, omega, M_0, geod_station, t_0_MJD, t_f_MJD, n_iter)
    % SIMULATEORBITFROMOE_WITHJ2 simulates an orbit from t_0 to t_f using the
    % provided orbital elements. 
    % It assumes Earth is the central body, and includes Earth's mean J2 effect.
    % It takes in orbital elements: a in [km], all angles in radians.
    % Times must be provided in MJD.
    % It returns position and velocity vectors in ECI and
    % ECEF, geodetic coordinates, position and velocity vectors in ENU,
    % azimuth and elevation and range, 
    % geod_vec contains [lat_vec; lon_vec; h_vec].
    % We can do this because J2 does not affect in plane motion, just RAAN,
    % omega, and slightly M.
    seconds_per_day = 86400;
    eps = 1e-12;
    mu = 3.986e5; % (km^3 / s^2) for earth
    R_E = 6378.1; % equatorial radius of earth in km
    J2 = 0.108263e-2; % for earth
    T = 2 * pi * sqrt(a^3 / mu);
    n = 2 * pi / T;
    [x_vec, y_vec, t_vec_sec, nu_vec, r_vec] = ...
        SimulatePerifocalOrbit_WithJ2(a, e, i, t_0_MJD*seconds_per_day, t_f_MJD*seconds_per_day, M_0, n_iter, eps);
    t_vec_MJD = t_vec_sec / seconds_per_day;
    r_ECI_vec = zeros(3,n_iter);
    v_ECI_vec = zeros(3,n_iter);
    r_ECEF_vec = zeros(3,n_iter);
    v_ECEF_vec = zeros(3,n_iter);
    geod_vec = zeros(3,n_iter);
    r_ENU_vec = zeros(3,n_iter);
    az_el_rho_vec = zeros(3,n_iter);
    r_RTN_vec = zeros(3,n_iter);
    v_RTN_vec = zeros(3,n_iter);
    for iter = 1:n_iter
        [r_ECI_vec(:,iter), v_ECI_vec(:,iter)] = OE2ECI(a, e, i, RAAN, omega, nu_vec(iter));
        if iter > 1 % skip the very first iteration
            % compute mean J2 effects
            dRAANdt = -(3/2) * n * J2 * (R_E / (a * (1-e^2)))^2 * cos(i);
            dwdt = (3/4) * n * J2 * (R_E / (a * (1-e^2)))^2 * (5 * cos(i)^2 - 1);
            % apply J2 updates
            dt = t_vec_sec(iter) - t_vec_sec(iter - 1);
            RAAN = RAAN + dRAANdt * dt;
            omega = omega + dwdt * dt;
        end
        GMST = MJDToGMST(t_vec_MJD(iter));
        [r_ECEF_vec(:,iter), v_ECEF_vec(:,iter)] = ECI2ECEF(r_ECI_vec(:,iter), v_ECI_vec(:,iter), GMST);
        [geod_vec(1,iter), geod_vec(2,iter), geod_vec(3,iter)] = ECEF2Geod(r_ECEF_vec(:,iter), eps);
        r_ENU_vec(:,iter) = ECEF2ENU(r_ECEF_vec(:,iter), geod_station);
        [az_el_rho_vec(1,iter), az_el_rho_vec(2,iter), az_el_rho_vec(3,iter)] = ENU2AzEl(r_ENU_vec(:,iter));
        R_ECI2RTN = rECI2RTN([r_ECI_vec(:,iter); v_ECI_vec(:,iter)]);
        r_RTN_vec(:,iter) = R_ECI2RTN * r_ECI_vec(:,iter);
        v_RTN_vec(:,iter) = R_ECI2RTN * v_ECI_vec(:,iter);
    end
    data_sim.r_ECI_vec = r_ECI_vec;
    data_sim.v_ECI_vec = v_ECI_vec;
    data_sim.r_ECEF_vec = r_ECEF_vec;
    data_sim.v_ECEF_vec = v_ECEF_vec;
    data_sim.geod_vec = geod_vec;
    data_sim.r_ENU_vec = r_ENU_vec;
    data_sim.az_el_rho_vec = az_el_rho_vec;
    data_sim.r_RTN_vec = r_RTN_vec;
    data_sim.v_RTN_vec = v_RTN_vec;
end

