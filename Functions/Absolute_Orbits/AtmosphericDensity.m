function rho = AtmosphericDensity(h, h_0, H, rho_0)
    % ATMOSPHERICDENSITY uses a simple exponential atmopheric density model
    % to convert a geocentric altitude to density. The model is the one
    % described in the lecture slides. 
    % h must be provided in [km].
    % rho is returned in [kg/m^3]. 
    if nargin == 1 % User should provide either only 1 input h, or all 4 inputs
        h_0 = 0; % sea level
        H = 10; % [km]
        rho_0 = 1.225; % [kg/m^3]
    else
        assert(nargin == 4, ...
            'AtmosphericDensity -- User must input 1 or 4 parameters, not any other number of input parameters.');
    end
    rho = rho_0 * exp(-((h - h_0) / H));
end
