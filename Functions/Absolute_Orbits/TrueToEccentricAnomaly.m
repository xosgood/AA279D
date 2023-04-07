function E = TrueToEccentricAnomaly(nu, e)
    % TRUETOECCENTRICANOMALY converts true anomaly to mean anomaly.
    % All angles must be in radians. 
    nu = wrapTo2Pi(nu);
    E = 2 * atan( sqrt((1-e)/(1+e)) * tan(nu/2) );
    if (nu >= 0) && (nu <= pi) % is in the top half plane
        E = wrapToPi(E);
    elseif (nu > pi) && (nu <= 2*pi) % is in the bottom half plane
        E = wrapTo2Pi(E);
    else
        error("Error: TrueToEccentricAnomaly -- Something wrong with wrapTo2Pi Function.");
    end
end
