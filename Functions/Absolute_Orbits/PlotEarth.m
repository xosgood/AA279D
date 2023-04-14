function PlotEarth()
    % plot the Earth
    rE = 6378.1; % km (equatorial Earth radius)
    [xE, yE, zE] = ellipsoid(0, 0, 0, rE, rE, rE, 20);
    surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black');
    axis equal; view(3); grid on; hold on;
end

