function r_ECEF = Geoc2ECEF(lat, lon, h)
    % GEOC2ECEF converts from geocentric latitude, longitude, and altitude
    % to a position vector in ECEF frame. 
    % Assumes Earth is sphere.
    % r_ECEF will be returned in km. 
    R_E = 6378.137; % [km]
    r_ECEF = (R_E + h) * [cos(lat) * cos(lon);
                          cos(lat) * sin(lon);
                               sin(lat)      ];
end
