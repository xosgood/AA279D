function delta_v = LS_control_solve(oe_c, delta_roe, u_controls)
    % oe_c -- mean elements of the chief.
    % roe_i -- intitial relative orbital elements.
    % roe_f -- final relative orbital elements. 
    % u_controls -- mean argument of lattitude to execute maneuver. 
    % u_f -- final mean argument of lattitude. We might take this out. 

    mu = 3.986e5; % (km^3 / s^2) for earth
    R_E = 6378.1; % equatorial radius of earth in km
    J2 = 0.108263e-2; % for earth

    N = length(u_controls);

    M = zeros(6, 3*N);

    a = oe_c(1);
    e = oe_c(2);
    i = oe_c(3);
    
    n = sqrt(a / mu^3);
    eta = sqrt(1 - e^2);
    kappa = 3 * J2 * R_E^2 * sqrt(mu) / (4 * a^3.5 * eta^4);
    P = 3*cos(i)^2;
    S = sin(2*i);
    Q = 5*cos(i)^2 -1; 
    E = 1 + eta;
    F = 4 + 3*eta;

    omega_dot = kappa*Q;


    for iter = 1:N
        A = zeros(6, 6);
        B = zeros(6, 3);

        % assume that mean argument of latitude of maneuever. 
        delta_u = u_controls(iter) - u_controls(1);

        tau = delta_u / (n + kappa*(eta * P + Q));

        A(1,1) = 1;
        A(2, 1) = - (7 * kappa * E * P + 3*n) * tau / 2;
        A(2, 2) = 1; 
        A(2, 5) = - kappa * F * S * tau;
        A(3, 3) = cos(omega_dot * tau);
        A(3, 4) = -sin(omega_dot * tau);
        A(3, 5) = 0;
        A(3, 6) = 0;
        A(4, 3) = sin(omega_dot * tau);
        A(4, 4) = cos(omega_dot * tau);
        A(5, 5) = 1;
        A(6, 1) = 3.5 * kappa * S * tau;

        B = zeros(6, 3);
        B(1, 2) = 2;
        B(2, 1) = -2;
        B(3, 1) = sin(u_controls(iter));
        B(3, 2) = 2*cos(u_controls(iter));
        B(4, 1) = -cos(u_controls(iter));
        B(4, 2) = 2*sin(u_controls(iter));
        B(5, 3) = cos(u_controls(iter));
        B(6, 3) = sin(u_controls(iter));

        B = 1/(n*a) * B;

        M(:, (1+3*(iter-1)):3*(iter)) = A*B;

    end

    disp(size(M));
    disp(size(delta_roe));

    delta_v = lsqr(M, delta_roe');

end