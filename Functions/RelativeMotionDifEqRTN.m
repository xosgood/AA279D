function statedot = RelativeMotionDifEqRTN(t, state)
    % State vector is [x_rtn y_rtn z_rtn theta0 r0 vx_rtn vy_rtn vz_rtn omega0 dr0] .
    
    mu = 3.986e5; % (km^3 / s^2) for earth

    x = state(1);
    y = state(2);
    z = state(3);
    theta0 = state(4);
    r0 = state(5);

    vx = state(6);
    vy = state(7);
    vz = state(8);
    omega0 = state(9);
    dr0 = state(10);
      
    
    dr0_dot = r0 * omega0^2 - mu/r0^2;
    omega0_dot = - 2*dr0*omega0/r0;

    vx_dot = 2*omega0*vy - omega0_dot*y + omega0^2 *x - mu*(r0+x) / ((r0 + x)^2 + y^2 + z^2)^(3/2) + mu/r0^2; 
    vy_dot = -2*omega0*vx - omega0_dot*x + omega0^2*y - mu*y/((r0 + x)^2 + y^2 + z^2)^(3/2);
    vz_dot = - mu*z / ((r0 + x)^2 + y^2 + z^2)^(3/2);

    statedot = [vx; vy; vz; omega0; dr0; vx_dot; vy_dot; vz_dot; omega0_dot; dr0_dot];
end