function PlotRTNPaxton(x_RTN)
    % Plot relative position in 3D. 
    figure;
    plot3(x_RTN(:,1), x_RTN(:,2), x_RTN(:,3));
    xlabel("R [km]");
    ylabel("T [km]");
    zlabel("N [km]");
    title("Relative position in RTN over time.")
    
    % Plot relative velocity in 3D. 
    figure;
    plot3(x_RTN(:,4), x_RTN(:,5), x_RTN(:,6));
    xlabel("R [km/s]");
    ylabel("T [km/s]");
    zlabel("N [km/s]");
    title("Relative velocity in RTN over time. ")
    
    % Plot relative position in TR, NR, TN plane.
    figure;
    subplot(3,1,1);
    plot(x_RTN(:,1), x_RTN(:,2));
    title("Relative position in RTN over time");
    xlabel("R [km]");
    ylabel("T [km]");
    subplot(3,1,2);
    plot(x_RTN(:,1), x_RTN(:,3));
    xlabel("R [km]");
    ylabel("N [km]");
    subplot(3,1,3);
    plot(x_RTN(:,2), x_RTN(:,3));
    xlabel("T [km]");
    ylabel("N [km]");
    
    
    % Plot relative velocity in TR, NR, TN plane.
    figure;
    subplot(3,1,1);
    plot(x_RTN(:,4), x_RTN(:,6));
    title("Relative velocity in RTN over time");
    xlabel("R [km/s]");
    ylabel("T [km/s]");
    subplot(3,1,2);
    plot(x_RTN(:,4), x_RTN(:,6));
    xlabel("R [km/s]");
    ylabel("N [km/s]");
    subplot(3,1,3);
    plot(x_RTN(:,5), x_RTN(:,6));
    xlabel("T [km/s]");
    ylabel("N [km/s]");
end