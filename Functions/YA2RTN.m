function x_RTN = YA2RTN(K, a, e, f, t)
    % YA2RTN converts true anomaly (and time) to RTN using the YA solutions to the TH
    % equations.
    % Arguments:
    %   K are the integration constants
    %   a, e, f are the semi-major axis, eccentricity, and true anomaly of
    %   the chief, respectively
    % Returns:
    %   x_RTN is the 6x1 position and velocity of the deputy in RTN,
    %       in the form [r_RTN; v_RTN].
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
    
    x_RTN = A * B * K;
end

