function roe = YAIC2ROE(K, oe)
    % YAIC2ROE Converts YA constants to quasi non-singular relative orbital
    % elements. 
    % Reference: Willis' PhD thesis equation 2.48.
    % Arguments:
    %   K: (6x1) integration constants of YA solution (from something like
    %       a function call to RTN2YA_IC), with true anomaly as independent variable.
    %   oe: (6x1) orbital elements of the chief. 
    % Returns:
    %   roe: (6x1) quasi non-singular orbital elements of deputy
    a = oe(1);
    e = oe(2);
    i = oe(3);
    RAAN = oe(4);
    w = oe(5);
    f = oe(6);
    
    [ex, ey] = DecomposeEccentricity(e, w);
    eta = sqrt(1 - e^2);
    
    R_f_to_u = [1, 0, 0;
                0, cos(w), sin(w);
                0, -sin(w), cos(w)];
    R = [R_f_to_u, zeros(3);
         zeros(3), R_f_to_u];
    K_u = R * K;
    
    A = [1, 0, 0, 0, 0, 0;
         0, -ex * (eta + 1/(1+eta)), ey * (eta + 1/(1+eta)), 1, 0, 0;
         0, ex * ey, ex^2 - 1, -ey, 0, -ey * cot(i);
         0, ey^2 - 1, ex * ey, ex, 0, ex * cot(i);
         0, 0, 0, 0, 1, 0;
         0, 0, 0, 0, 0, -1];
    
    roe = A * K_u;
end

