function nu = MeanToTrueAnomaly(M, e)
    % MEANTOTRUEANOMALY Converts mean anomaly to true anomaly
    nu = EccentricToTrueAnomaly(MeanToEccentricAnomaly(M, e, 1e-12), e);
end

