function [r_ECI, v_ECI] = OE2ECI(a, e, i, RAAN, omega, nu)
    % OE2ECI takes in 6 Keplerian orbital elements and returns the
    % postition vector and velocity vector in ECI coordinate frame.
    % It assumes the central body is the Earth.
    % 'a' must be provided in km. Angles must be provided in radians. 
    mu = 3.986e5; % (km^3 / s^2)
    E = TrueToEccentricAnomaly(nu, e);
    n = sqrt(mu / a^3);
    r_PQW = [a * (cos(E) - e);
             a * sqrt(1-e^2) * sin(E);
             0];
    v_PQW = ((a * n) /  (1 - e * cos(E))) * [-sin(E);
                                             sqrt(1-e^2) * cos(E);
                                             0];
    %%% Note: for some stupid reason rotz and rotx work in degrees ugh
            % and it does the active rotation instead of
            % passive, so i need to flip all the negative
            % signs, double ugh, I'll write my own functions in
            % the future
    i_deg = rad2deg(i);
    RAAN_deg = rad2deg(RAAN);
    omega_deg = rad2deg(omega);
    if i == 0 && e == 0         % equatorial and circular -> don't use RAAN or omega
        R_PQW_to_ECI = 1;
    elseif i == 0 && e ~= 0     % equatorial and non-circular -> don't use RAAN or i
        R_PQW_to_ECI = rotz(omega_deg);
    elseif e == 0 && i ~= 0     % circular and inclined 
        R_PQW_to_ECI = rotz(RAAN_deg) * rotx(i_deg); 
    else                        % non-circular and inclined
        R_PQW_to_ECI = rotz(RAAN_deg) * rotx(i_deg) * rotz(omega_deg); 
    end

    r_ECI = R_PQW_to_ECI * r_PQW;
    v_ECI = R_PQW_to_ECI * v_PQW;
end
