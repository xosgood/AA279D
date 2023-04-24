function t = TrueAnomalyToTime(nu, a, e)
    % TRUEANOMALYTOTIME Converts a true anomaly to time since periapsis.
    mu = 3.986e5;
    n = sqrt(mu / a^3);
    
    M = TrueToMeanAnomaly(nu, e);
    t = M / n;
end

