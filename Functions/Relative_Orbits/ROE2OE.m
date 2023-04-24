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
    
    % initialize deputy elements
    oe_d = zeros(size(oe_c));
    
    dRAAN = roe_d(5) / a_c * sin(i_c);
    
%     /* Compute difference in RAAN */
%     dO = roe(5) / (ac * sin(ic));
%     /* Compute difference in mean argument of latitude */
%     du = (roe(1) / ac) - (dO * cos(ic));
% 
%     /* Compute deputy keplerian orbit elements */
%     ad = ac + roe(0);
%     exd = exc + (roe(2) / ac);
%     eyd = eyc + (roe(3) / ac);
%     id = ic + (roe(4) / ac);
%     Od = Oc + dO;
%     ud = uc + du;
% 
%     orbit_d = { ad, exd, eyd, id, Od, ud };
end

