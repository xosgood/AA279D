function r_ENU = ECEF2ENU(r_sat_ECEF, geod_station)
    % ECEF2ENU converts a position vector from ECEF coordinate frame to
    % ENU coordinate frame. 
    % Position of sat must be provided in km.
    % geoc_station contains the geocentric latitude, longitude, and altitude, 
    % and must be in the form [lat, lon, h] in radians and km respectively.
    % Returns a position vector in ENU frame in km.
    lat = geod_station(1);
    lon = geod_station(2);
    h = geod_station(3);
    r_station_ECEF = Geod2ECEF(lat, lon, h);
    
    % remember rotz and rotx use degrees and do active instead of passive rotations
    Rot_ECEF_to_ENU = rotx(-rad2deg(pi/2 - lat)) * rotz(-rad2deg(lon + pi/2)); 
%     Rot = [     -sin(lon)      ,       cos(lon)       ,     0;
%            -sin(lat) * cos(lon), -sin(lat) * sin(lon) , cos(lat);
%             cos(lat) * cos(lon),   cos(lat) * sin(lon), sin(lat)];
    r_ENU = Rot_ECEF_to_ENU * (r_sat_ECEF - r_station_ECEF);
end

