function M = true2mean(f, e)
% true2mean solves Kepler's equation for mean anomaly
%
% Inputs:
%     f - true anomaly [rad]
%     e - eccentricity of orbit
%
% Outputs:
%     M - mean anomaly [rad]

E = true2ecc(f,e);
M = ecc2mean(E,e);

end