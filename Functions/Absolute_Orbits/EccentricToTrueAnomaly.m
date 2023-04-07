function nu = EccentricToTrueAnomaly(E, e)
    % ECCENTRICTOTRUEANOMALY converts from eccentric anomaly to mean
    % anomaly. All angles must be in radians. 
    E = wrapTo2Pi(E);
    nu = 2 * atan( sqrt((1+e)/(1-e)) * tan(E/2) );
    if (E >= 0) && (E <= pi) % is in the top half plane
        nu = wrapToPi(nu);
    elseif (E > pi) && (E <= 2*pi) % is in the bottom half plane
        nu = wrapTo2Pi(nu);
    else
        error("Error: EccentricToTrueAnomaly -- Something wrong with wrapTo2Pi Function.");
    end
end
