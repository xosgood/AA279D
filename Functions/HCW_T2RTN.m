function x_RTN = HCW_T2RTN(K, a, t)
    % K - the integration constants.
    % a - the initial semimajor axis.
    % t - (t-t_0).
    
    mu = 3.986e5;
    n = sqrt(mu/a^3);

    A = [a*eye(3,3), zeros(3,3); zeros(3,3), a*n*eye(3,3)];

    B = [1, sin(n*t), cos(n*t), 0, 0, 0;
        -3/2*n*t, 2*cos(n*t), -2*sin(n*t), 1, 0, 0;
        0, 0, 0, 0, sin(n*t), cos(n*t); 
        0, cos(n*t), -sin(n*t), 0, 0, 0;
        -3/2, -2*sin(n*t), -2*cos(n*t), 0, 0, 0;
        0, 0, 0, 0, cos(n*t), -sin(n*t)];

    x_RTN = A*B*K;
end