function PlotQNSROE_meters(roe_qns, a_c)
    % Plot x/y for eccentricity vector.
    % Plot x/y for inclination vector.
    % Plot dlambda/da.
	% oe_qns is (6xN).
    % a_c is the semi-major axis of the chief.
    da = roe_qns(1,:);
    dlambda = roe_qns(2,:);
    dex = roe_qns(3,:);
    dey = roe_qns(4,:);
    dix = roe_qns(5,:);
    diy = roe_qns(6,:);
    
    subplot(3,1,1);
    plot(a_c * dex, a_c * dey);
    axis equal;
    hold on;
    xlabel("a \delta e_x [m]");
    ylabel("a \delta e_y [m]");
    title("Relative eccentricity vector");
    
    subplot(3,1,2);
    plot(a_c * dix, a_c * diy);
    axis equal;
    hold on;
    xlabel("a \delta i_x [m]");
    ylabel("a \delta i_y [m]");
    title("Relative inclination vector");

    subplot(3,1,3);
    plot(a_c * dlambda, a_c * da);
    axis equal;
    hold on;
    xlabel("a \delta \lambda [m]");
    ylabel("a \delta a [m]");
    title("Relative mean longitude and semi-major axis");
end

