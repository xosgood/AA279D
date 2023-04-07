function GMST = MJDToGMST(MJD_UT1)
    % MJDTOGMST takes in UT1 time in MJD format and converts it to
    % Greenwich Mean Sidereal Time, expressed in radians. 
    % Note that this formula uses an approximation that doesn't account for
    % small secular variations in the Earth's rotation. 
    d = MJD_UT1 - 51544.5;
    GMST = 280.4606 + 360.9856473 * d; % (degrees)
    GMST = deg2rad(GMST);
    GMST = wrapTo2Pi(GMST);
end
