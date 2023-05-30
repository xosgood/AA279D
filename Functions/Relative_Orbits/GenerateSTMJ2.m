function A = GenerateSTMJ2(oe_c, tau)

    % For arguments takes:
    %   oe_c - kepler orbital elements of chief.
    %   roe_qns - relative orbital elements of deputy relative to chief.
    
    mu = 3.986e5; % (km^3 / s^2) for earth
    R_E = 6378.1; % equatorial radius of earth in km
    J2 = 0.108263e-2; % for earth

    a = oe_c(1);
    e = oe_c(2);
    i = oe_c(3);
    RAAN = oe_c(4);
    omega = oe_c(5);
    nu = oe_c(6);

    n = sqrt(1 - e^2);
    kappa = 3 * J2 * R_E^2 * sqrt(mu) / (4 * a^3.5 * n^4);
    E = 1 + n;
    P = 3*cos(i)^2;
    Q = 5*cos(i)^2 -1; 
    S = sin(2*i);
    T = sin(i)^2;
    F = 4 + 3*n;
    G = 1/n^2;
    
    omega_dot = kappa*Q;
    omega_f = omega + omega_dot*tau;

    e_xi = e*cos(omega + nu*RAAN);
    e_yi = e*sin(omega + nu*RAAN);
    e_xf = e*cos(omega_f + nu*RAAN);
    e_yf = e*sin(omega_f + nu*RAAN);
   
    A_21 = -(1.5*n + 3.5*kappa*E*P*tau);
    A_23 = kappa*e_xi*e_yf*G*Q*tau;
    A_24 = kappa*e_xi*e_yf*G*Q*tau;
    A_25 = -kappa*F*S*tau; 

    A_31 = 3.5*kappa*e_yf *Q*tau;
    A_33 = cos(omega_dot *tau) - 4*kappa*e_xi*e_yf*G*Q*tau;
    A_34 = -sin(omega_dot *tau) - 4*kappa*e_yi*e_yf*G*Q*tau;
    A_35 = 5*kappa*e_yf*S*tau;

    A_41 = -3.5 * kappa * e_xf * Q * tau;
    A_43 = sin(omega_dot * tau) + 4*kappa*e_xi*e_xf*G*Q*tau;
    A_44 = cos(omega_dot * tau) + 4*kappa*e_yi*e_xf*G*Q*tau;
    A_45 = -5*kappa*e_xf*S*tau;

    A_61 = 3.5 *kappa * S * tau;
    A_63 = -4*kappa*e_xi*G*S*tau;
    A_64 = -4*kappa*e_yi*G*S*tau;
    A_65 = 2*kappa*T*tau;

    A = [1, 0, 0, 0, 0, 0; 
         A_21, 1, A_23, A_24, A_25, 0;
         A_31, 0, A_33, A_34, A_35, 0;
         A_41, 0, A_43, A_44, A_45, 0;
         0, 0, 0, 0, 1, 0;
         A_61, 0, A_63, A_64, A_65, 1];


end