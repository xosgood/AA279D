function [x1_ECI] = ECI2RTN_Vectorized(x0_ECI, x1_RTN)
    % Arugments:
    %   - x0_ECI: intertial state in the ECI frame of cheif.
    %   - x1_RTN: relative RTN position of deputy to chief.
    % Returns
    %   - x1_ECI: inertial state of deputy in ECI frame.  

    % ECI2RTN_VECTORIZED Turn time series of ECI Chief states and RTN deputy states to ECI deputy states.
    % Each vector should be size (3xN) where N is the number of timesteps

    x1_ECI = zeros(size(x0_ECI));

    for iter = 1:size(x0_ECI,2)
        [x1_ECI(:,iter)] = RTN2ECI(x0_ECI(:,iter), x1_RTN(:,iter));
    end
end