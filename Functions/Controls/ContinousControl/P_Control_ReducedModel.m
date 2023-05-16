function P = P_Control_ReducedModel(oe_c, Delta_roe_reducedmodel, N, k)
    % Compute the P proportional gain matrix for Lyapunov continous control.
    % As it pertains to the Reduced Model (removing dlambda and dv_R).
    % We will use the definition of P as defined in Steindorf, et al.
    %   Lyapunov Artificial Potential paper ("CONSTRAINED LOW-THRUST SATELLITE FORMATION-FLYING USING RELATIVE ORBIT ELEMENTS".
    % Arguments:
    %   oe_c: (6x1) keplerian orbital elements of the chief [a, e, i, RAAN, omega, nu]^T
    %   Delta_roe_reducedmodel: (5x1) change in roe of deputy from current
    %       state to reference state.
    %   N: scalar exponent of cos() functions, must be a positive even integer.
    %   k: inverse of scalar to be placed in front P matrix, must be positive.
    % Returns:
    %   P: (5x5) proportional gain matrix as defined by Steindorf, et al.
    
    % unpack chief oe
    e = oe_c(2);
    w = oe_c(5);
    f = oe_c(6);
    M = TrueToMeanAnomaly(f, e);
    
    % unpack deputy roe
    Delta_dex = Delta_roe_reducedmodel(2);
    Delta_dey = Delta_roe_reducedmodel(3);
    Delta_dix = Delta_roe_reducedmodel(4);
    Delta_diy = Delta_roe_reducedmodel(5);
    
    phi = w + M; % mean argument of latitude
    
    phi_IP = atan2(Delta_dey, Delta_dex); % ideal in-plane maneuver location
    phi_OOP = atan2(Delta_diy, Delta_dix); % ideal out-of-plane maneuver location
    
    J = phi - phi_IP;
    H = phi - phi_OOP;
    
    IP_entry = cos(J)^N;
    OOP_entry = cos(H)^N;
    
    P = (1/k) * diag([IP_entry, IP_entry, IP_entry, OOP_entry, OOP_entry]);
end

