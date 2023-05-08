function [r_RTN, v_RTN] = ECI2RTN_Vectorized(r_ECI_0, v_ECI_0, r_ECI_1, v_ECI_1)
    % ECI2RTN_VECTORIZED Turn time series of ECI states to RTN states.
    % Each vector should be size (3xN) where N is the number of timesteps
    r_RTN = zeros(size(r_ECI_1));
    v_RTN = zeros(size(r_ECI_1));
    for iter = 1:size(r_ECI_0,2)
        [r_RTN(:,iter), v_RTN(:,iter)] = ECI2RTN(r_ECI_0(:,iter), v_ECI_0(:,iter), r_ECI_1(:,iter), v_ECI_1(:,iter));
    end
end

