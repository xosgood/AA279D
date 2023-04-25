function x_RTN = GeometricMapping_Linear(roe, oe, M_0, dM_0)
    % GEOMETRICMAPPING_LINEAR Linear geometric mapping from relative
    % orbital elements to position and velocity in RTN. For eccentric
    % orbits. Linearized (approximation from Taylor expansion. Valid for
    % small ROEs. 
    % Arguments:
    %   roe: (6x1) relative orbital elements (quasi-non singular) of deputy
    %       - [da, dlambda, dex, dey, dix, diy]^T
    %   oe: (6x1) orbital elements of chief
    %   M_0: initial mean anomaly of the chief [radians]
    %   dM_0: initial delta mean anomaly between chief and deputy
    % Returns:
    %   x_RTN: (6x1) position and velocity in RTN
    %% unpack arguments
    da = roe(1);
    dl = roe(2);
    dex = roe(3);
    dey = roe(4);
    dix = roe(5);
    diy = roe(6);
    
    a = oe(1);
    e = oe(2);
    i = oe(3);
    RAAN = oe(4);
    w = oe(5);
    f = oe(6);
    
    %% getting position in RTN
    [r_ECI_c, v_ECI_c] = OE2ECI(a, e, i, RAAN, w, f);
    r_c = norm(r_ECI_c);
    
    [ex, ey] = DecomposeEccentricity(e, w);
    
    u = f + w;
    k = 1 + ex * cos(u) + ey * sin(u);
    kp = -ex * sin(u) + ey * cos(u);
    eta = sqrt(1 - e^2);
    
    x_bar = da - k * kp / eta^3 * dl - dex / eta^3 * k * cos(u) - dex / eta^3 * k * sin(u) ...
          + k / eta^3 * ((k + 1) / (1 + eta)) * (ex * dex + ey * dey) + k * kp / eta^3 * diy * cot(i);
    y_bar = k^2 / eta^3 * dl + dex / eta^2 * (1 + k) * sin(u) - dey / eta^2 * (1 + k) * cos(u) ...
          + 1 / eta^2 * (eta + k^2 / (1+eta)) * (ey * dex + ex * dey) + (1 - k^2 / eta^3) * diy * cot(i);
    z_bar = dix * sin(u) - diy * cos(u);
    
    x = x_bar * r_c;
    y = y_bar * r_c;
    z = z_bar * r_c;
    
    %% getting velocity in RTN
    dw = atan2(dey, dex) - w;
    de = sqrt(dex^2 + dey^2);
    di = dix; % sqrt(dix^2 + diy^2); % ????????
    dRAAN = diy / sin(i);
    
    % compute dM (delta mean anomaly)
    % dM = dl - dw - diy * cot(i); % WRONG
    M = TrueToMeanAnomaly(f, e);
    dM = dM_0 + 3/2 * da/a * (M - M_0);

    d_u = sqrt(e^2 * dM^2 / eta^2 + de^2);
    d_w = sqrt(di^2 + sin(i)^2 * dRAAN^2);
    f_u = atan2(e * dM, -eta * de);
    % f_v = f_u - pi/2;
    theta_w = atan2(di, -sin(i) * dRAAN);
    theta = u; % ??????????
    
    u = da / a - e * de / (2 * eta^2) + d_u / eta^2 * (cos(f - f_u) + e/2 * cos(2*f - f_u));
    v = (1 - e^2 / 2) * dM / eta^3 + dw + cos(i) * dRAAN - d_u / eta^2 * (2 * sin(f - f_u) + e/2 * sin(2*f - f_u));
    w = d_w * cos(theta - theta_w);
    
    %% pack return variable
    x_RTN = [x; y; z; u; v; w];
end

