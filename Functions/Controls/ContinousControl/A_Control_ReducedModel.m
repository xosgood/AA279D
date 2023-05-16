function A = A_Control_ReducedModel(oe_c)
    %A_CONTROL_REDUCEDMODEL Summary of this function goes here
    %   Detailed explanation goes here
    
    mu = 3.986e5; % (km^3 / s^2) for earth
    R_E = 6378.1; % equatorial radius of earth in km
    J2 = 0.108263e-2; % for earth
    
    a = oe_c(1);
    e = oe_c(2);
    i = oe_c(3);
    RAAN = oe_c(4);
    omega = oe_c(5);
    nu = oe_c(6);
    
    n = sqrt(mu / a^3);
    eta = sqrt(1 - e^2);
    kappa = 3 * J2 * R_E^2 * sqrt(mu) / (4 * a^3.5 * eta^4);
    
    e_x = e*cos(omega);
    e_y = e*sin(omega);
    
    C = sin(omega);
    D = cos(omega);
    E = 1 + eta;
    F = 4 + 3*eta;
    G = 1 / eta^2;  
    P = 3*cos(i)^2 - 1;
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

