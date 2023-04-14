function R_ECI2RTN = rECI2RTN(vec_ECI)
    %ECI2RTN convert a vector in ECI to RTN frame
    % vec_ECI should be [r; v]_ECI
    assert(all(size(vec_ECI) == [6, 1]), "Vector provided to rECI2RTN() must be size 6x1")
    
    r = vec_ECI(1:3);
    v = vec_ECI(4:6);
    n = cross(r, v);

    R = r / norm(r);
    N = r / norm(n);
    T = cross(N, R);

    R_ECI2RTN = [R, T, N]';
end

