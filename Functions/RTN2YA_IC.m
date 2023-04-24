function K = RTN2YA_IC(x_RTN, a, e, f, t)
    % RTN2YA_IC Convert position and velocity in RTN to integration
    % constants K using the Yamanaka-Ankersen (YA) solutions to the
    % Tschauner-Hempel (TH) equations for short range elliptical orbits.
    % Arguments:
    %   x_RTN must be a 6x1 vector [r_RTN; v_RTN].
    %   e is the eccentricity of the chief orbit.
    %   f is the true anomaly of the chief orbit.
    % Returns:
    %   K will be a 6x1 vector.
    mu = 3.986e5;
    n = sqrt(mu/a^3);
    k = 1 + e * cos(f);
    kp = -e * sin(f);
    tau = n * t / eta^3;
    eta = sqrt(1 - e^2);
    A = [a * eta^2 * eye(3), zeros(3);
         zeros(3), a * n / eta];
    B = [1/k + 3/2 * kp * tau, 
end

