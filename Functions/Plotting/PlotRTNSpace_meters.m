function PlotRTNSpace_meters(x_RTN)
%%% Input shall be given in meters
    % Plot relative position in 3D. 
    figure;
    plot3(x_RTN(:,1), x_RTN(:,2), x_RTN(:,3));
    xlabel("R [m]");
    ylabel("T [m]");
    zlabel("N [m]");
    %title("Relative position in RTN over time.")
    axis equal;
    
    % Plot relative velocity in 3D. 
%     figure;
%     plot3(x_RTN(:,4), x_RTN(:,5), x_RTN(:,6));
%     xlabel("R [m/s]");
%     ylabel("T [m/s]");
%     zlabel("N [m/s]");
%     title("Relative velocity in RTN over time. ")
%     axis equal;
    
    % Plot relative position in TR, NR, TN plane.
    figure;
    subplot(3,1,1);
    plot(x_RTN(:,2), x_RTN(:,1));
    axis equal;
    %title("Relative position in RTN over time");
    xlabel("T [m]");
    ylabel("R [m]");
    subplot(3,1,2);
    plot(x_RTN(:,3), x_RTN(:,1));
    axis equal;
    xlabel("N [m]");
    ylabel("R [m]");
    subplot(3,1,3);
    plot(x_RTN(:,2), x_RTN(:,3));
    axis equal;
    xlabel("T [m]");
    ylabel("N [m]");
    
    % Plot relative velocity in TR, NR, TN plane.
%     figure;
%     subplot(3,1,1);
%     plot(x_RTN(:,5), x_RTN(:,4));
%     axis equal;
%     title("Relative velocity in RTN over time");
%     xlabel("T [m/s]");
%     ylabel("R [m/s]");
%     subplot(3,1,2);
%     plot(x_RTN(:,6), x_RTN(:,4));
%     axis equal;
%     xlabel("N [m/s]");
%     ylabel("R [m/s]");
%     subplot(3,1,3);
%     plot(x_RTN(:,5), x_RTN(:,6));
%     axis equal;
%     xlabel("T [m/s]");
%     ylabel("N [m/s]");
end