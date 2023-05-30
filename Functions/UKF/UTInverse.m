function [mu, Sigma] = UTInverse(points, weights)
    % UTInverse is the Inverse of the Unscented Transform.
    % Arguments:
    %   points: (n x 2n+1) sigma-points resulting from the unscented transform
    %   weights: (1 x 2n+1) weights corresponding to each sigma-point
    % Returns:
    %   mu: (n x 1) mean vector
    %   Sigma: (n x n) covariance matrix
    mu = sum(points .* weights, 2);
    Sigma = (weights .* (points - mu)) * (points - mu)';
end

