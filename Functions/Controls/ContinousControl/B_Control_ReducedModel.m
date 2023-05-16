function B = B_Control_ReducedModel(oe_c)
    %B_CONTROL_REDUCEDMODEL Summary of this function goes here
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
    
    B = zeros(5, 2);

    B(1, 1) = 2/eta * (1 + e*cos(nu));
    
    B(2, 1) = eta * ((2 + e*cos(nu)) * cos(omega + nu) + e_x) / (1 + e*cos(nu));
    B(2, 2) = (eta * e_y)/ tan(i) * sin(w + nu) / (1 + e*cos(nu));

    B(3, 1) = eta * ((2 + e*cos(nu)) * sin(omega+nu) + e) /(1 + e*cos(nu));
    B(3, 2) = (-eta * e_x)/ tan(i) * sin(w + nu) / (1 + e*cos(nu));

    B(4, 2) = eta * cos(omega + nu) / (1 + e*cos(nu));
    B(5, 2) = eta * sin(omega + nu) / (1 + e*cos(nu));

    B = B/(a*n);
end 

