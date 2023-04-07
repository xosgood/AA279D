function R_CRF_to_TRF = RotCRFToTRF(GMST)
    % ROTCRFTOTRF takes in GMST in radians, and outputs a rotation matrix
    % that will rotate a vector from Celestial Reference Frame (CRF) to 
    % Terrestrial Reference Frame (TRF).
    % Note that this function neglects nutation, precession, polar motion,
    % and equation of equinoxes for the Earth. 
    R_CRF_to_TRF = rotz(-rad2deg(GMST));
end
