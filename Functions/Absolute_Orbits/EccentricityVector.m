function e = EccentricityVector(mu, r, v)
    e_r = dot(v, v)/mu - 1/norm(r);        %vector component in the r direction.
    e_v = -dot(r, v)/mu;                    %vector compoenent in the v direction. 

    e= e_r.*r + e_v.*v;
end