function PlotRTNSpace(x_RTN)
    % Plot relative position in 3D. 
    figure;
    plot3(x_RTN(:,1), x_RTN(:,2), x_RTN(:,3));
    grid on;
    xlabel("R [km]");
    ylabel("T [km]");
    zlabel("N [km]");
    title("Relative position in RTN over time.")
    axis equal;
    
    % Plot relative velocity in 3D. 
    figure;
    plot3(x_RTN(:,4), x_RTN(:,5), x_RTN(:,6));
    grid on;
    xlabel("R [km/s]");
    ylabel("T [km/s]");
    zlabel("N [km/s]");
    title("Relative velocity in RTN over time. ")
    axis equal;
    
    % Plot relative position in TR, NR, TN plane.
    figure;
    subplot(3,1,1);
    plot(x_RTN(:,2), x_RTN(:,1));
    axis equal; grid on;
    title("Relative position in RTN over time");
    xlabel("T [km]");
    ylabel("R [km]");
    subplot(3,1,2);
    plot(x_RTN(:,3), x_RTN(:,1));
    axis equal; grid on;
    xlabel("N [km]");
    ylabel("R [km]");
    subplot(3,1,3);
    plot(x_RTN(:,2), x_RTN(:,3));
    axis equal; grid on;
    xlabel("T [km]");
    ylabel("N [km]");
    
    % Plot relative velocity in TR, NR, TN plane.
    figure;
    subplot(3,1,1);
    plot(x_RTN(:,5), x_RTN(:,4));
    axis equal; grid on;
    title("Relative velocity in RTN over time");
    xlabel("T [km/s]");
    ylabel("R [km/s]");
    subplot(3,1,2);
    plot(x_RTN(:,6), x_RTN(:,4));
    axis equal; grid on;
    xlabel("N [km/s]");
    ylabel("R [km/s]");
    subplot(3,1,3);
    plot(x_RTN(:,5), x_RTN(:,6));
    axis equal; grid on;
    xlabel("T [km/s]");
    ylabel("N [km/s]");
end
