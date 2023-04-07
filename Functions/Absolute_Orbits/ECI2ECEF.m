function [r_ECEF, v_ECEF] = ECI2ECEF(r_ECI, v_ECI, GMST)
    % ECI2ECEF converts position and velocity vectors in the ECI frame to
    % the ECEF frame. 
    R_ECI_to_ECEF = rotz(-rad2deg(GMST));
    r_ECEF = R_ECI_to_ECEF * r_ECI;
    v_ECEF = R_ECI_to_ECEF * v_ECI;
end

