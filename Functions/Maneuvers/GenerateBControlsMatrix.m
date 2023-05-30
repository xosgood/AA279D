function B = GenerateBControlsMatrix(u, oe)
    % Returns a matrix that maps between delta v and ROE elements.
    % Arguments:
    %   oe: (6x1) orbital elements of chief [a, e, i, RAAN, omega, nu]^T.
    
    a = oe(1);
    mu = 3.986e5;
    n = sqrt(mu / a^3);

    B = [0, 2, 0;
         -2, 0, 0;
         sin(u), 2*cos(u), 0;
         -cos(u), 2*sin(u), 0;
         0, 0, cos(u);
         0, 0, sin(u)] / (n * a);

end
