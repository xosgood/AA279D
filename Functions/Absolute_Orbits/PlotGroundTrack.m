function PlotGroundTrack(lat_vec, lon_vec)
    % PLOTGROUNDTRACK plots a ground track given lat lon data (in radians)
    
    % Load and plot MATLAB built-in Earth topography data
    figure;
    load('topo.mat', 'topo');
    topoplot = [topo(:, 181:360), topo(:, 1:180)];
    contour(-180:179, -90:89, topoplot, [0, 0], 'black');
    axis equal;
    grid on;
    xlim([-180, 180]);
    ylim([-90, 90]);
    xlabel('Longitude [\circ]');
    ylabel('Latitude [\circ]');
    
    hold on;
    
    plot(rad2deg(lon_vec), rad2deg(lat_vec), 'b.');
    
end

