function PlotConfidenceInterval(orbit_span, x, mu, Sigma, a_c)
    % x is ground truth, mu is state estimate. 
    % a_c is chief's semi-major axis (for scaling purposes).
    a = 1.96; % 1.96 std dev is 95% confidence interval
    
    x = a_c * x;
    mu = a_c * mu;
    Sigma = a_c^2 * Sigma;
    
    orbit_span_conf = [orbit_span orbit_span(end:-1:1)];
    mu_conf = [mu(1,:) + a * sqrt(reshape(Sigma(1,1,:), 1,length(Sigma(1,1,:)))), mu(1,end:-1:1) - a * sqrt(reshape(Sigma(1,1,end:-1:1), 1,length(Sigma(1,1,end:-1:1))));
               mu(2,:) + a * sqrt(reshape(Sigma(2,2,:), 1,length(Sigma(2,2,:)))), mu(2,end:-1:1) - a * sqrt(reshape(Sigma(2,2,end:-1:1), 1,length(Sigma(2,2,end:-1:1))));
               mu(3,:) + a * sqrt(reshape(Sigma(3,3,:), 1,length(Sigma(3,3,:)))), mu(3,end:-1:1) - a * sqrt(reshape(Sigma(3,3,end:-1:1), 1,length(Sigma(3,3,end:-1:1))));
               mu(4,:) + a * sqrt(reshape(Sigma(4,4,:), 1,length(Sigma(4,4,:)))), mu(4,end:-1:1) - a * sqrt(reshape(Sigma(4,4,end:-1:1), 1,length(Sigma(4,4,end:-1:1))));
               mu(5,:) + a * sqrt(reshape(Sigma(5,5,:), 1,length(Sigma(5,5,:)))), mu(5,end:-1:1) - a * sqrt(reshape(Sigma(5,5,end:-1:1), 1,length(Sigma(5,5,end:-1:1))));
               mu(6,:) + a * sqrt(reshape(Sigma(6,6,:), 1,length(Sigma(6,6,:)))), mu(6,end:-1:1) - a * sqrt(reshape(Sigma(6,6,end:-1:1), 1,length(Sigma(6,6,end:-1:1))))];

    figure; grid on; hold on;
    sgtitle("Estimation error vs. time");
    ylabels = ["a\delta a [km]", "a\delta \lambda [km]", "a\delta e_x [km]", "a\delta e_y [km]", "a\delta i_x [km]", "a\delta i_y [km]"];
    
    for i = 1:6
        subplot(3,2,i);
        grid on; hold on;
        plot(orbit_span, x(i,:));
        plot(orbit_span, mu(i,:));
        
        p = fill(orbit_span_conf, mu_conf(i,:), 'm');
        p.FaceAlpha = 0.2;
        p.EdgeColor = 'none';
        
        xlabel("orbits");
        ylabel(ylabels(i));
    end
    
    subplot(3,2,2);
    legend("True", "Estimated", "95% confidence interval");
end


