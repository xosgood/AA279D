function [e_x, e_y] = DecomposeEccentricity(e, w)
    % DECOMPOSEECCENTRICITY Decomposes the eccentricity into e_x and e_y,
    % using the argument of periapsis.
    e_x = e * cos(w);
    e_y = e * sin(w);
end

