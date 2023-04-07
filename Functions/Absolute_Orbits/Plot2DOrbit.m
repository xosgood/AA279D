function [x,y] = Plot2DOrbit(p, e)
    nus = deg2rad([0:1:360]);
    r = p ./ (1 + e * cos(nus));
    k = find(r>0);
    r_valid = r(k);
    nus_valid = nus(k);
    [x,y] = pol2cart(nus_valid, r_valid);
    plot(x,y);
    grid on;
    grid minor;
    axis equal;
    xlabel("x");
    ylabel("y");
end