function r_ECEF = Geod2ECEF(lat, lon, h)
    % GEOC2ECEF converts from geodetic latitude, longitude, and altitude
    % to a position vector in ECEF frame. 
    % Assumes Earth is oblate ellipsoid with eccentricity e_E.
    % r_ECEF will be returned in km. 
    R_E = 6378.137; % [km]
    e_E = 0.0818;
    N = R_E / sqrt(1 - e_E^2 * sin(lat)^2);
    r_ECEF = [(N + h) * cos(lat) * cos(lon);
              (N + h) * cos(lat) * sin(lon);
              (N * (1 - e_E^2) + h) * sin(lat)];
end
