function f = ThirdBodyAccel(mu, r_13, r_23)
    % THIRDBODYACCEL
    % mu is that of the third body, in [km^3/s^2].
    % r_13 is the inertial position vector from the central body (1) to 
    % the third body, [km]. 
    % r_23 is the inertial position vector from the satellite (2) to the 
    % third body, [km].
    f = mu * ((r_23 / norm(r_23)^3) - (r_13 / norm(r_13)^3));
end

