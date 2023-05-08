function PlotOEvsTime(t, oe, ylabels)
    % PLOTOEVSTIME Plot oe over time. 6 subplots. oe should be a (6xN).
    for i = 1:6
        subplot(6,1,i);
        plot(t, oe(i,:));
        hold on;
        ylabel(ylabels(i));
    end
    xlabel("time (s)");
end

