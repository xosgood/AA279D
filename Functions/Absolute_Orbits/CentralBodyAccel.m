function f = CentralBodyAccel(mu, r)
    % CENTRALBODYACCEL
    f = -(mu / norm(r)^3) * r;
end

