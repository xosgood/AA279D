function [lat, lon, h] = ECEF2Geoc(r_ECEF)
    % ECEF2GEOC converts a position vector in ECEF frame to geocentric coordinates. 
    % r_ECEF must be given in km. 
    % lat and lon will be returned in radians.
    % h will be returned in km.
    % h is an altitude about the geocentric radius of the Earth. 
    R_Earth = 6378.137; % (km)
    r_x = r_ECEF(1);
    r_y = r_ECEF(2);
    r_z = r_ECEF(3);
    lat = atan2(r_z, sqrt(r_x^2 + r_y^2));
    lon = atan2(r_y, r_x);
    h = norm(r_ECEF) - R_Earth;
end
