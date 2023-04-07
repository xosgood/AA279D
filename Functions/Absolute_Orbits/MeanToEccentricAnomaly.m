function E = MeanToEccentricAnomaly(M, e, eps)
    % MEANTOECCENTRICANOMALY converts mean anomaly to eccentric anomaly
    % using Newton-Raphson method. 
    % M is given in radians. e is eccentricity and must be in the range [0, 1). 
    % eps is tolerance that the Newton-Raphson method should converge at. 
    % E is returned in radians from [0, 2pi]
    M = wrapTo2Pi(M);
    iter_lim = 100; % maximum number of iterations to perform
    i = 1;
    
    % make initial guess
    if (e < 1.0) && (e >= 0.8)
        E_old = pi;
    elseif (e >= 0.0) && (e < 0.8) % 0 <= e < 0.9
        E_old = M;
    else
        error("Error: MeanToEccentricAnomaly -- Invalid Eccentricity Provided");
    end
    
    err = eps + 1; % initialize error
    E = 0; % initialize E
    % start loop
    while(err > eps && i <= iter_lim)
        E = E_old - ( (E_old - e * sin(E_old) - M) / (1 - e * cos(E_old)) );
        err = abs(E - E_old);
        E_old = E;
        i = i + 1;
    end
    
    % case that it did not converge in iter_lim
    if i > iter_lim
        error("Error: MeanToEccentricAnomaly -- Newton's method did not converge");
    end
    
    E = wrapTo2Pi(E); % just a precaution
end

