function PlotRTN(t, r_RTN, v_RTN)
    % PLOTRTN 
    % Each vector should be size (3xN) where N is the number of timesteps
    
    %%% Plot state vs time
    
    % plot positions in RTN
    subplot(3,3,1);
    plot(t, r_RTN(1,:));
    ylabel("position in R [km]");
    xlabel("time [s]");
    subplot(3,3,4);
    plot(t, r_RTN(2,:));
    ylabel("position in T [km]");
    xlabel("time [s]");
    subplot(3,3,7);
    plot(t, r_RTN(3,:));
    ylabel("position in N [km]");
    xlabel("time [s]");
    % plot velocities in RTN
    subplot(3,3,2);
    plot(t, v_RTN(1,:));
    ylabel("velocity in R [km/s]");
    xlabel("time [s]");
    subplot(3,3,5);
    plot(t, v_RTN(2,:));
    ylabel("velocity in T [km/s]");
    xlabel("time [s]");
    subplot(3,3,8);
    plot(t, v_RTN(3,:));
    ylabel("velocity in N [km/s]");
    xlabel("time [s]");
    
    %%% Plot relative motion in projected planes
    % plot TR
    subplot(3,3,3);
    plot(r_RTN(2,:), r_RTN(1,:));
    xlabel("T [km]");
    ylabel("R [km]");
    axis equal
    % plot NR
    subplot(3,3,6);
    plot(r_RTN(3,:), r_RTN(1,:));
    xlabel("N [km]");
    ylabel("R [km]");
    axis equal
    % plot TN
    subplot(3,3,9);
    plot(r_RTN(2,:), r_RTN(3,:));
    xlabel("T [km]");
    ylabel("N [km]");
    axis equal
end

