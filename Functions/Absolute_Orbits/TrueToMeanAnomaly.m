function M = TrueToMeanAnomaly(nu, e)
    % TRUETOMEANANOMALY converts from mean anomaly to true anomaly. 
    % All angles must be in radians. 
    E = TrueToEccentricAnomaly(nu, e);
    M = EccentricToMeanAnomaly(E, e);
    M = wrapTo2Pi(M);
end
