function data_sim = SimulateOrbitFromOE(a, e, i, RAAN, omega, M_0, geod_station, t_0_MJD, t_f_MJD, n_iter)
    % SIMULATEORBITFROMOE simulates an orbit from t_0 to t_f using the
    % provided orbital elements. 
    % It assumes Earth is the central body. 
    % It takes in orbital elements: a in [km], all angles in radians.
    % Times must be provided in MJD.
    % It returns position and velocity vectors in ECI and
    % ECEF, geodetic coordinates, position and velocity vectors in ENU,
    % azimuth and elevation and range, 
    % geod_vec contains [lat_vec; lon_vec; h_vec].
    % It 
    seconds_per_day = 86400;
    eps = 1e-12;
    [x_vec, y_vec, t_vec_sec, nu_vec, r_vec] = ...
        SimulatePerifocalOrbit(a, e, t_0_MJD*seconds_per_day, t_f_MJD*seconds_per_day, M_0, n_iter, eps);
    t_vec_MJD = t_vec_sec / seconds_per_day;
    r_ECI_vec = zeros(3,n_iter);
    v_ECI_vec = zeros(3,n_iter);
    r_ECEF_vec = zeros(3,n_iter);
    v_ECEF_vec = zeros(3,n_iter);
    geod_vec = zeros(3,n_iter);
    r_ENU_vec = zeros(3,n_iter);
    az_el_rho_vec = zeros(3,n_iter);
    for iter = 1:n_iter
        [r_ECI_vec(:,iter), v_ECI_vec(:,iter)] = OE2ECI(a, e, i, RAAN, omega, nu_vec(iter));
        GMST = MJDToGMST(t_vec_MJD(iter));
        [r_ECEF_vec(:,iter), v_ECEF_vec(:,iter)] = ECI2ECEF(r_ECI_vec(:,iter), v_ECI_vec(:,iter), GMST);
        [geod_vec(1,iter), geod_vec(2,iter), geod_vec(3,iter)] = ECEF2Geod(r_ECEF_vec(:,iter), eps);
        r_ENU_vec(:,iter) = ECEF2ENU(r_ECEF_vec(:,iter), geod_station);
        [az_el_rho_vec(1,iter), az_el_rho_vec(2,iter), az_el_rho_vec(3,iter)] = ENU2AzEl(r_ENU_vec(:,iter));
    end
    data_sim.r_ECI_vec = r_ECI_vec;
    data_sim.v_ECI_vec = v_ECI_vec;
    data_sim.r_ECEF_vec = r_ECEF_vec;
    data_sim.v_ECEF_vec = v_ECEF_vec;
    data_sim.geod_vec = geod_vec;
    data_sim.r_ENU_vec = r_ENU_vec;
    data_sim.az_el_rho_vec = az_el_rho_vec;
end

