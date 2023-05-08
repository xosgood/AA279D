function [x1_ECI] = ECI2RTN(x0_ECI, x1_RTN)
    % Arugments:
    %   - x0_ECI: intertial state in the ECI frame of cheif.
    %   - x1_RTN: relative RTN position of deputy to chief.
    % Returns
    %   - x1_ECI: inertial state of deputy in ECI frame. 

    R_rtn2eci = rRTN2ECI(x0_ECI);

    % Transform position.
    r_ECI = x0_ECI(1:3);
    r_RTN = x1_RTN(1:3);
    r1_ECI = R_rtn2eci*r_RTN + r_ECI;

    % Transform velocity.
    v_ECI = x0_ECI(4:6);
    v_RTN = x1_RTN(4:6);
    f_dot = norm(cross(r_ECI,v_ECI)) / norm(r_ECI)^2;
    omega = [0; 0; f_dot];
    v1_ECI = R_rtn2eci*(v_RTN + cross(omega, r_RTN)) + v_ECI;

    x1_ECI = [r1_ECI; v1_ECI];
    
end