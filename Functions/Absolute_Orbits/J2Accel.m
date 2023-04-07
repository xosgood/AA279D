function f = J2Accel(mu, J2, R, r, i, u)
    % J2ACCEL
    % R is the equatorial radius of the central body.
    % r is the radius of the orbit, the norm of the position vector, in [km].
    % i is inclination, and u is argument of latitude, both in radians.
    % Returns acceleration due to J2 in the RTN frame.
    term = 3 * mu * J2 * R^2 / (2 * r^4);
    f_R = -term * (1 - 3 * sin(i)^2 * sin(u)^2);
    f_T = -term * sin(i)^2 * sin(u) * cos(u);
    f_N = -term * sin(i) * cos(i) * sin(u);
    f = [f_R; f_T; f_N];
end

