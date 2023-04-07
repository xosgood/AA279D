function f = J2AccelECI(mu, J2, a, x, y, z)
    r = sqrt(x^2 + y^2 + z^2);
    g_x = - mu * x / r^3 * (1 + 3/2 * J2 * a^2 / r^2 - 15/2 * J2 * a^2 * z^2 / r^4);
    g_y = - mu * y / r^3 * (1 + 3/2 * J2 * a^2 / r^2 - 15/2 * J2 * a^2 * z^2 / r^4);
    g_z = - mu * z / r^3 * (1 + 9/2 * J2 * a^2 / r^2 - 15/2 * J2 * a^2 * z^2 / r^4);
    
    f = [g_x; g_y; g_z];
end