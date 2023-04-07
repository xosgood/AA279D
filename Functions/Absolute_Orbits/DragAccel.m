function f = DragAccel(r_ECI, v_ECI, B)
    % DRAGACCEL takes in ECI position and velocity vectors as well as a 
    % ballistic coefficient and converts them to an acceleration due to 
    % drag force f_Drag in ECI frame.
    % This function assumes the Earth is spherical and that the Earth and 
    % its atmosphere are one rigid body, i.e. the atmosphere rotates 
    % perfectly with the Earth.
    % If position and velocity are given in [km] and [km/s], respectively, 
    % this function will then return acceleration in [km/s^2]. 
    % Note B must be provided in [m^2/kg]
    R_E = 6378.137; % [km]
    omega_E = [0; 0; 7.2921e-5]; % angular rotation rate of Earth [rad/s]
    h = norm(r_ECI) - R_E;
    rho = AtmosphericDensity(h); % in [kg/m^3]
    v_rel = v_ECI - cross(omega_E, r_ECI);
    v_rel = v_rel * 1000; % convert [km/s] to [m/s]
    f = -0.5 * B * rho * norm(v_rel) * v_rel;
    f = f / 1000; % convert [m/s^2] to [km/s^2]
end

