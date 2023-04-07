function [az, el, rho] = ENU2AzEl(r_ENU)
    % ENU2AZEL converts a position vector in ENU frame to azimuth and
    % elevation angles, as well as distance rho.
    % Angles will be in radians.
    r_E = r_ENU(1);
    r_N = r_ENU(2);
    r_U = r_ENU(3);
    rho = norm(r_ENU);
    az = atan2(r_E, r_N);
    el = atan2(r_U, sqrt(r_E^2 + r_N^2));
end

