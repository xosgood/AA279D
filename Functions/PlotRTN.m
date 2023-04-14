function PlotRTN(t, r_RTN, v_RTN)
    % PLOTRTN 
    % Each vector should be size (3xN) where N is the number of timesteps
    
    %%% Plot state vs time
    
    % plot positions in RTN
    subplot(3,2,1);
    plot(t, r_RTN(1,:));
    ylabel("position in R [km]");
    subplot(3,2,3);
    plot(t, r_RTN(2,:));
    ylabel("position in T [km]");
    subplot(3,2,5);
    plot(t, r_RTN(3,:));
    ylabel("position in N [km]");
    xlabel("time [s]");
    % plot velocities in RTN
    subplot(3,2,2);
    plot(t, v_RTN(1,:));
    ylabel("velocity in R [km/s]");
    subplot(3,2,4);
    plot(t, v_RTN(2,:));
    ylabel("velocity in T [km/s]");
    subplot(3,2,6);
    plot(t, v_RTN(3,:));
    ylabel("velocity in N [km/s]");
    xlabel("time [s]");
    
    %%% Plot relative motion in projected planes
    % TODO
    % plot RT
    
    % plot RN
    
    % plot TN
    
end

