function B = GenerateBControlsMatrix(oe)
    % Returns a matrix that maps between delta v and ROE elements.
    % Arguments:
    %   oe: (6x1) orbital elements of chief [a, e, i, RAAN, omega, nu]^T.
    mu = 3.986e5;

    a = oe(1);
    n = sqrt(mu / a^3);
    omega = oe(5);
    nu = oe(6);
    u = omega + nu;

    B = [0, 2, 0;
         -2, 0, 0;
         sin(u), 2*cos(u), 0;
         -cos(u), 2*sin(u), 0;
         0, 0, cos(u);
         0, 0, sin(u)] / (n * a);
end
