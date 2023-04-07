function [lat_d, lon_d, h_d] = Geoc2Geod(lat_c, lon_c, h_c)
    % GEOC2GEOD converts geocentric latitude, longitude, and height to
    % geodetic latitude, longitude, and height. 
    % Assumes central body is Earth.
    % All lat and lon in radians. 
    % All heights in km.
    r_ECEF = Geoc2ECEF(lat_c, lon_c, h_c);
    [lat_d, lon_d, h_d] = ECEF2Geod(r_ECEF, 1e-12);
end
