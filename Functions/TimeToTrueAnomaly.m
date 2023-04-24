function nu = TimeToTrueAnomaly(t, e)
    % TIMETOTRUEANOMALY Converts time from periapsis to true anomaly.
    mu = 3.986e5;
    n = sqrt(mu/a^3);
    
    M = n * t;
    nu = MeanToTrueAnomaly(M, e);
end

