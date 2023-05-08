function [r_RTN, v_RTN] = ECI2RTN(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1)
    % ECI2RTN converts a state (position and velocity) from ECI to RTN.
    % "_0" denotes chief, "_1" denotes deputy.
    % All vectors must be provided as column vectors.
    % Reference: https://sisl.github.io/SatelliteDynamics.jl/latest/modules/reference_systems/#SatelliteDynamics.sECItoRTN
    
    % rotation matrix from ECI to RTN
    R_ECI2RTN = rECI2RTN([r_ECI_0; v_ECI_0]);
    
    % turn position into RTN
    rho = r_ECI_1 - r_ECI_0;
    r_RTN = R_ECI2RTN * rho;
    
    % turn velocity into RTN
    f_dot = norm(cross(r_ECI_0, v_ECI_0)) / norm(r_ECI_0)^2;
    omega = [0; 0; f_dot];
    rho_dot = v_ECI_1 - v_ECI_0;
    v_RTN = R_ECI2RTN * rho_dot - cross(omega, r_RTN);
end

