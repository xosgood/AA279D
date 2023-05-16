function A = A_Control_ReducedModel(oe_c)
    % Returns the reduced plant matrix given the orbital elements of the
    % chief. 

    % We will use the definition of A as defined in Steindorf, et al.
    %   Lyapunov Artificial Potential paper ("CONSTRAINED LOW-THRUST SATELLITE FORMATION-FLYING USING RELATIVE ORBIT ELEMENTS".

    % Arguments:
    %   oe_c: (6x1) keplerian orbital elements of the chief [a, e, i, RAAN, omega, nu]^T
    
    mu = 3.986e5; % (km^3 / s^2) for earth
    R_E = 6378.1; % equatorial radius of earth in km
    J2 = 0.108263e-2; % for earth
    
    a = oe_c(1);
    e = oe_c(2);
    i = oe_c(3);
    omega = oe_c(5);
    
    eta = sqrt(1 - e^2);
    kappa = 3 * J2 * R_E^2 * sqrt(mu) / (4 * a^3.5 * eta^4);
    
    e_x = e*cos(omega);
    e_y = e*sin(omega);
    
    C = sin(omega);
    D = cos(omega);
    G = 1 / eta^2;
    Q = 5*cos(i)^2 -1; 
    S = sin(2*i);
    T = sin(i)^2;
    
    A = zeros(5, 5);
    
    A(2, 1) = 7/2 * e_y * Q;
    A(2, 2) = - (4* e_x * e_y * G + C)*Q;
    A(2, 3) = - (1 + 4*e_y^2 *G - D)*Q;
    A(2, 4) = 5*e_y*S; 
    
    A(3, 1) = -7/2 *e_x *Q;
    A(3, 2) = (1 + 4*e_x^2*G - D)*Q;
    A(3, 3) = (4*e_x*e_y*G -C)*Q;
    A(3, 4) = -5*e_x * S;
    
    A(5, 1) = 7/2 *S; 
    A(5, 2) = -4*e_x*G*S;
    A(5, 3) = -4*e_y*G*S;
    A(5, 4) = 2*T;

    A = kappa*A;

end

