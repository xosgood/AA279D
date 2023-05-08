function PlotRTNSpaceMultiple(x_RTN, last_call)
    % Plot relative position in 3D. 
    figure(1);
    plot3(x_RTN(:,1), x_RTN(:,2), x_RTN(:,3));
    xlabel("R [km]");
    ylabel("T [km]");
    zlabel("N [km]");
    title("Relative position in RTN over time.")
    if (~last_call)
        hold on;
    end
    
    % Plot relative velocity in 3D. 
    figure(2);
    plot3(x_RTN(:,4), x_RTN(:,5), x_RTN(:,6));
    xlabel("R [km/s]");
    ylabel("T [km/s]");
    zlabel("N [km/s]");
    title("Relative velocity in RTN over time. ")
    if (~last_call)
        hold on;
    end
    
    % Plot relative position in TR, NR, TN plane.
    figure(3);
    subplot(3,1,1);
    plot(x_RTN(:,2), x_RTN(:,1));
    title("Relative position in RTN over time");
    xlabel("T [km]");
    ylabel("R [km]");
    if (~last_call)
        hold on;
    end
    subplot(3,1,2);
    plot(x_RTN(:,3), x_RTN(:,1));
    xlabel("N [km]");
    ylabel("R [km]");
    if (~last_call)
        hold on;
    end
    subplot(3,1,3);
    plot(x_RTN(:,2), x_RTN(:,3));
    xlabel("T [km]");
    ylabel("N [km]");
    if (~last_call)
        hold on;
    end
    
    % Plot relative velocity in TR, NR, TN plane.
    figure(4);
    subplot(3,1,1);
    plot(x_RTN(:,5), x_RTN(:,4));
    title("Relative velocity in RTN over time");
    xlabel("T [km/s]");
    ylabel("R [km/s]");
    if (~last_call)
        hold on;
    end
    subplot(3,1,2);
    plot(x_RTN(:,6), x_RTN(:,4));
    xlabel("N [km/s]");
    ylabel("R [km/s]");
    if (~last_call)
        hold on;
    end
    subplot(3,1,3);
    plot(x_RTN(:,5), x_RTN(:,6));
    xlabel("T [km/s]");
    ylabel("N [km/s]");
    if (~last_call)
        hold on;
    end
end