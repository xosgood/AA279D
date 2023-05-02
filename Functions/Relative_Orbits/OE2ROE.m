function roe = OE2ROE(oe_c, oe_d)
    % OE2ROE converts orbital elements of a deputy to quasi non-singular
    % relative orbital elements, based off of the chief's orbital elements.
    % All angles must be provided in radians.
    % oe in form [a, e, i, RAAN, omega, nu]^T.
    % roe in form [da, dlambda, dex, dey, dix, diy]^T.
    
    % Pull out chief elements
    a_c = oe_c(1);
    e_c = oe_c(2);
    i_c = wrapToPi(oe_c(3));
    RAAN_c = wrapToPi(oe_c(4));
    omega_c = wrapToPi(oe_c(5));
    nu_c = oe_c(6);
    M_c = TrueToMeanAnomaly(nu_c, e_c);
    
    % Pull out deputy elements
    a_d = oe_d(1);
    e_d = oe_d(2);
    i_d = wrapToPi(oe_d(3));
    RAAN_d = wrapToPi(oe_d(4));
    omega_d = wrapToPi(oe_d(5));
    nu_d = oe_d(6);
    M_d = TrueToMeanAnomaly(nu_d, e_d);
    
    dM = wrapToPi(M_d - M_c);
    domega = wrapToPi(omega_d - omega_c);
    dRAAN = wrapToPi(RAAN_d - RAAN_c);
    
    % Form relative orbital elements
    da = (a_d - a_c) / a_c;
    dlambda = dM + domega + dRAAN * cos(i_c);
    dex = e_d * cos(omega_d) - e_c * cos(omega_c);
    dey = e_d * sin(omega_d) - e_c * sin(omega_c);
    dix = wrapToPi(i_d - i_c);
    diy = dRAAN * sin(i_c);
    
    % Pack up relative orbital element
    roe = [da; dlambda; dex; dey; dix; diy];
end

