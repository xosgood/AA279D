function oe_d = ROE2OE(oe_c, roe_d)
    % ROE2OE converts relative orbital elements of deputy to keplerian orbital
    % elements.
    % oe in form [a, e, i, RAAN, omega, nu]^T
    % roe in form [da, dlambda, dex, dey, dix, diy]^T
    
    % Pull out cheif elements
    a_c = oe_c(1);
    e_c = oe_c(2);
    i_c = oe_c(3);
    RAAN_c = oe_c(4);
    omega_c = oe_c(5);
    nu_c = oe_c(6);
    M_c = TrueToMeanAnomaly(nu_c, e_c);

    [ex_c, ey_c] = DecomposeEccentricity(e_c, omega_c);

    % Argument of lattitude.
    u_c = nu_c + omega_c;
    
    da = roe_d(1);
    dlambda = roe_d(2);
    dex = roe_d(3);
    dey = roe_d(4);
    dix = roe_d(5);
    diy = roe_d(6);
    
    a_d = a_c + a_c * da;
    
    i_d = i_c + dix;
    
    ex_d = ex_c + dex;
    ey_d = ey_c + dey;
    e_d = sqrt(ex_d^2 + ey_d^2);
    
    omega_d = atan2(ey_d, ex_d);
    
    dRAAN = diy / sin(i_c);
    RAAN_d = RAAN_c + dRAAN;
    
    domega = omega_d - omega_c;
    dM = dlambda - domega - dRAAN * cos(i_c);
    M_d = M_c + dM;
    nu_d = MeanToTrueAnomaly(M_d, e_d);
    
    oe_d = [a_d; e_d; i_d; RAAN_d; omega_d; nu_d];
end

