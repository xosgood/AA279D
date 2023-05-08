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
    eta = sqrt(1 - e^2);
    tau = n * t / eta^3;
    A = [a * eta^2 * eye(3), zeros(3);
         zeros(3), a * n / eta * eye(3)];
    B = [1/k + 3/2 * kp * tau, sin(f), cos(f), 0, 0, 0;
         -3/2 * k * tau, (1 + 1/k) * cos(f), -(1 + 1/k) * sin(f), 1/k, 0, 0;
         0, 0, 0, 0, 1/k * sin(f), 1/k * cos(f);
         kp/2 - 3/2 * k^2 * (k-1) * tau, k^2 * cos(f), -k^2 * sin(f), 0, 0, 0;
         -3/2 * (k + k^2 * kp * tau), -(k^2 - 1) * sin(f), -e - (k^2 + 1)*cos(f), -kp, 0, 0;
         0, 0, 0, 0, e + cos(f), -sin(f)];
    
    K = (A*B) \ x_RTN;
end

