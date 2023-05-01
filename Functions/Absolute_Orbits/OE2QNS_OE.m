function oe_qns = OE2QNS_OE(oe)
    a = oe(1);
    e = oe(2);
    i = oe(3);
    RAAN = oe(4);
    omega = oe(5);
    nu = oe(6);

    u = omega + nu;
    ex = e * cos(omega);
    ey = e * sin(omega);

    oe_qns = [a, u, ex, ey, i, RAAN];
    
end
