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

    ex_c = e_c * cos(omega_c);
    ey_c = e_c * sin(omega_c);

    % Argument of lattitude.
    u_c = nu_c + omega_c;
    
    % initialize deputy elements
    oe_d = zeros(size(oe_c));
    
    dRAAN = roe_d(5) / a_c * sin(i_c);
    du = (roe_d(1) / a_c) - (dRAAN * cos(i_c));

    a_d = a_c + roe_d(1);
    ex_d = ex_c + roe_d(3) / a_c;
    ey_d = ey_c + roe_d(4) / a_c;
    i_d = i_c + roe_d(5) / a_c;
    RAAN_d = RAAN_c + dRAAN;
    u_d = u_c + du;
    
    % e as a function of e_x and e_y.
    e_d = sqrt(ex_d^2 + ey_d^2);

    % omega as a function of ey_d, ex_d.
    omega_d = atan2(ey_d, ex_d);

    % true anomoly as a function of argument of lattitude and preapsis. 
    nu_d = u_d - omega_d;

    oe_d = [a_d; e_d; i_d; RAAN_d; omega_d; nu_d];

    
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

