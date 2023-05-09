function roe_new = ApplyDeputyManuever_NearCircular(oe, roe, dv)
    % Apply an impulsive delta-v manuever to deputy.
    % Arguments:
    %   oe: (6x1) orbital elements of chief [a, e, i, RAAN, omega, nu]^T.
    %   roe: (6x1) quasi non-singular relative orbital elements of deputy
    %       [da, dlambda, dex, dey, dix, diy]^T.
    %   dv: (3x1) delta-v, expressed in RTN frame of deputy
    roe = roe(:); % make sure roe is a column vector
    dv = dv(:); % make sure dv is a column vector
    
    mu = 3.986e5;
    
    a = oe(1);
    e = oe(2);
    i = oe(3);
    RAAN = oe(4);
    omega = oe(5);
    nu = oe(6);
    [ex, ey] = DecomposeEccentricity(e, omega);
    eta = sqrt(1 - e^2);
    u = omega + nu;
    n = sqrt(mu / a^3);
    
%     Gamma = [2 / eta * e * sin(nu), 2 / eta * (1 + e*cos(nu)), 0;
%              -2 * eta^2 / (1 + e * cos(nu)), 0, 0;
%              eta*sin(theta), eta*((2 + e*cos(nu))*cos(theta) + ex) / (1+e*cos(nu)), eta*ey*sin(theta) / (tan(i)*(1+e*cos(nu)));
%              -eta*cos(theta), eta*((2 + e*cos(nu))*sin(theta) + ey) / (1+e*cos(nu)), -eta*ex*sin(theta) / (tan(i)*(1+e*cos(nu)));
%              0, 0, eta * cos(theta) / (1 + e * cos(nu));
%              0, 0, eta * sin(theta) / (1 + e * cos(nu))] / (n * a * 1000);
    Gamma = [0, 2, 0;
             -2, 0, 0;
             sin(u), 2*cos(u), 0;
             -cos(u), 2*sin(u), 0;
             0, 0, cos(u);
             0, 0, sin(u)] / (n * a * 1000);
    
    delta_roe = Gamma * dv;
    
    roe_new = roe + delta_roe;
end
