function x_RTN = GeometricMapping_Linear(roe, oe, t)
    % GEOMETRICMAPPING_LINEAR Linear geometric mapping from relative
    % orbital elements to position and velocity in RTN. For eccentric
    % orbits. Linearized (approximation from Taylor expansion). Valid for
    % small ROEs. 
    % Formula comes from equation 2.47 in Willis's PhD Thesis.
    % Arguments:
    %   roe: (6x1) relative orbital elements (quasi-non singular) of deputy
    %       - [da, dlambda, dex, dey, dix, diy]^T
    %   oe: (6x1) orbital elements of chief
    %   t: time since initial time
    % Returns:
    %   x_RTN: (6x1) position and velocity in RTN
    a = oe(1);
    e = oe(2);
    i = oe(3);
    RAAN = oe(4);
    w = oe(5);
    f = oe(6);
    
    [ex, ey] = DecomposeEccentricity(e, w);
    
    mu = 3.986e5;
    n = sqrt(mu / a^3);
    u = f + w;
    k = 1 + ex * cos(u) + ey * sin(u);
    kp = -ex * sin(u) + ey * cos(u);
    eta = sqrt(1 - e^2);
    
    A = [a * eta^2 * eye(3), zeros(3);
         zeros(3), a * n / eta * eye(3)];
    
    bx1 = 1/k + 3/2 * kp * n / eta^3 * t;
    bx2 = -kp / eta^3;
    bx3 = (ex * ((k - 1) / (1 + eta)) - cos(u)) / eta^3;
    bx4 = (ey * ((k - 1) / (1 + eta)) - sin(u)) / eta^3;
    bx6 = kp / eta^3 * cot(i);
    by1 = -3/2 * k * n / eta^3 * t;
    by2 = k / eta^3;
    by3 = ((1 + 1/k) * sin(u) + ey/k + k/eta * (ey / (1+eta))) / eta^2;
    by4 = -((1 + 1/k) * cos(u) + ex/k + k/eta * (ex / (1+eta))) / eta^2;
    by6 = (1 / k - k / eta^3) * cot(i);
    bz5 = sin(u) / k;
    bz6 = - cos(u) / k;
    bxd1 = kp/2 + 3/2 * k^2 * (1 - k) * n / eta^3 * t;
    bxd2 = k^2 / eta^3 * (k - 1);
    bxd3 = k^2 / eta^3 * (eta * sin(u) + ey * ((k - 1) / (1 + eta)));
    bxd4 = -k^2 / eta^3 * (eta * cos(u) + ex * ((k - 1) / (1 + eta)));
    bxd6 = -k^2 / eta^3 * (k - 1) * cot(i);
    byd1 = -3/2 * k * (1 + k * kp * n / eta^3 * t);
    byd2 = k^2 / eta^3 * kp;
    byd3 = (1 + k^2/eta^3) * cos(u) + ex*k/eta^2 * (1 + k/eta * ((k-1) / (1+eta)));
    byd4 = (1 + k^2/eta^3) * sin(u) + ey*k/eta^2 * (1 + k/eta * ((k-1) / (1+eta)));
    byd6 = -(1 + k^2/eta^3) * kp * cot(i);
    bzd5 = cos(u) + ex;
    bzd6 = sin(u) + ey;
    
    B = [bx1, bx2, bx3, bx4, 0, bx6;
         by1, by2, by3, by4, 0, by6;
         0, 0, 0, 0, bz5, bz6;
         bxd1, bxd2, bxd3, bxd4, 0, bxd6;
         byd1, byd2, byd3, byd4, 0, byd6;
         0, 0, 0, 0, bzd5, bzd6];
    
    x_RTN = A * B * roe;
end

