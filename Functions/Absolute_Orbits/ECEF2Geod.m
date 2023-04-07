function [lat_d, lon_d, h_d] = ECEF2Geod(r_ECEF, eps)
    % ECEF2GEOD converts a position vector in ECEF coordinates to 
    % latitude, longitude, and altitude in geodetic coordinates. 
    % Specifically, this function uses an oblate ellipsoid model of the Earth
    % with an eccentricity of 0.0818.
    % All lat and lon variables will be returned in radians.
    % eps is tolerance that the iterative method should converge at. 
    R_E = 6378.137; % equatorial radius of Earth (km)
    e_E = 0.0818; % eccentricity of ellipsoid of Earth
    r_x = r_ECEF(1);
    r_y = r_ECEF(2);
    r_z = r_ECEF(3);
    
    [lat_c, lon_c] = ECEF2Geoc(r_ECEF);
    
    lon_d = lon_c; % longitude stays the same
    
    lat_d = lat_c ; % initial guess for lattitude
    N = 0; % initialize
    err = eps + 1;
    while err > eps
        lat_d_prev = lat_d;
        N = R_E / sqrt(1 - e_E^2 * sin(lat_d)^2);
        lat_d = atan2(r_z + N * e_E^2 * sin(lat_d), sqrt(r_x^2 + r_y^2));
        err = abs(lat_d - lat_d_prev);
    end
    h_d = sqrt(r_x^2 + r_y^2) / cos(lat_d) - N;
end
