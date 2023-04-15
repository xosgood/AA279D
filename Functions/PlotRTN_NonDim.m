function PlotRTN_NonDim(t, r_RTN, v_RTN, a)
    %PLOTRTN_NONDIM
    r_RTN = r_RTN / a;
    
    % plot positions in RTN
    subplot(3,3,1);
    plot(t, r_RTN(1,:));
    ylabel("x/a (R)");
    xlabel("time [s]");
    subplot(3,3,4);
    plot(t, r_RTN(2,:));
    ylabel("y/a (T)");
    xlabel("time [s]");
    subplot(3,3,7);
    plot(t, r_RTN(3,:));
    ylabel("z/a (N)");
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
    % plot RT
    subplot(3,3,3);
    plot(r_RTN(1,:), r_RTN(2,:));
    xlabel("x/a (R)");
    ylabel("y/a (T)");
    % plot RN
    subplot(3,3,6);
    plot(r_RTN(1,:), r_RTN(3,:));
    xlabel("x/a (R)");
    ylabel("z/a (N)");
    % plot TN
    subplot(3,3,9);
    plot(r_RTN(2,:), r_RTN(3,:));
    xlabel("y/a (T)");
    ylabel("z/a (N)");
end

