function M = EccentricToMeanAnomaly(E, e)
    % ECCENTRICTOMEANANOMALY converts from eccentric anomaly ('E') to mean
    % anomaly ('M'), using the eccentricity ('e') of the orbit. 
    % All angles must be in radians.
    M = E - e * sin(E);
end

