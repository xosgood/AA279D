function f = J2AccelECI(J2, mu, R_E, x, y, z)
    r = norm([x, y, z]);
    
    g_x = - mu * x / r^3 * (-(3/2)*J2*(R_E/r)^2 * (5 * z^2 / r^2 - 1));
    g_y = - mu * y / r^3 * (-(3/2)*J2*(R_E/r)^2 * (5 * z^2 / r^2 - 1));
    g_z = - mu * z / r^3 * (-(3/2)*J2*(R_E/r)^2 * (5 * z^2 / r^2 - 3));
    
    f = [g_x; g_y; g_z];
end