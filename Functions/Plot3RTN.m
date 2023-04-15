function Plot3RTN(r_RTN, a)
    %PLOT3RTN Plots RTN position in 3D.
    % Input vector should be size (3xN) where N is the number of timesteps.
    % If a is zero, or negative, then plot dimensionalized quantity,
    % otherwise plot nondimensionalized quantities, nondimensionalized by a.
    if a > 0
        r_RTN = r_RTN / a;
        xstr = "x/a (R)";
        ystr = "y/a (T)";
        zstr = "z/a (N)";
    else
        xstr = "x (R) [km]";
        ystr = "y (T) [km]";
        zstr = "z (N) [km]";
    end
    plot3(r_RTN(1,:), r_RTN(2,:), r_RTN(3,:));
    xlabel(xstr);
    ylabel(ystr);
    zlabel(zstr);
end

