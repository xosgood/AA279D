function f = SRPAccel(r, B_SRP)
    % SRPACCEL computes the accleration vector due to solar radiation
    % pressure. 
    % Assumes satellite is around Earth.
    % r should be in [km].
    % B_SRP = (C_SRP * A_SRP) / mass of sat, in [m^2 / kg]
    % Returns acceleration in [km/s]
    p_SRP = 4.57e-6; % [N/m^2]
    f = p_SRP * B_SRP * (r / norm(r));
    f = f / 1000; % convert to [km/s]
end

