function B = B_Control_ReducedModel(oe_c)
    % Returns the reduced control input matrix given 
    % the orbital elements of the chief. 

    % We will use the definition of B as defined in Steindorf, et al.
    %   Lyapunov Artificial Potential paper ("CONSTRAINED LOW-THRUST SATELLITE FORMATION-FLYING USING RELATIVE ORBIT ELEMENTS".

    % Arguments:
    %   oe_c: (6x1) keplerian orbital elements of the chief [a, e, i, RAAN, omega, nu]^T

    mu = 3.986e5; % (km^3 / s^2) for earth
    
    a = oe_c(1);
    e = oe_c(2);
    i = oe_c(3);
    w = oe_c(5);
    nu = oe_c(6);
    
    n = sqrt(mu / a^3);
    eta = sqrt(1 - e^2);
    
    e_x = e*cos(w);
    e_y = e*sin(w);
    
    B = zeros(5, 2);

    B(1, 1) = 2/eta * (1 + e*cos(nu));
    
    B(2, 1) = eta * ((2 + e*cos(nu)) * cos(w + nu) + e_x) / (1 + e*cos(nu));
    B(2, 2) = (eta * e_y)/ tan(i) * sin(w + nu) / (1 + e*cos(nu));

    B(3, 1) = eta * ((2 + e*cos(nu)) * sin(w + nu) + e_y) /(1 + e*cos(nu));
    B(3, 2) = (-eta * e_x)/ tan(i) * sin(w + nu) / (1 + e*cos(nu));

    B(4, 2) = eta * cos(w + nu) / (1 + e*cos(nu));
    B(5, 2) = eta * sin(w + nu) / (1 + e*cos(nu));

    B = B/(a*n);
end 

