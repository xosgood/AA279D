function [points, weights] = UT(mu, Sigma, lambda)
    % UT is the Unscented Transform (aka sigma-point transform).
    % Arguments:
    %   mu: (n x 1) mean vector
    %   Sigma: (n x n) covariance matrix
    %   lambda: scalar hyperparameter that helps define the weights
    % Returns:
    %   points: (n x 2n+1) sigma-points resulting from the unscented transform
    %   weights: (1 x 2n+1) weights corresponding to each sigma-point
    if ~exist('lambda','var') % make lambda an optional argument
        lambda = 2; % default lambda
    end
    
    mu = mu(:); % ensure mu is a column vector
    
    n = length(mu);
    
    Sigma_sqrt = chol(Sigma)'; % Cholesky decomposition to get matrix square root
    
    points = zeros(n, 2 * n + 1);
    weights = zeros(1, 2 * n + 1);
        
    points(:,1) = mu;
    weights(1) = lambda / (lambda + n);
    for i = 2:n+1
        points(:,i) = mu + sqrt(lambda + n) * Sigma_sqrt(:,i-1);
        points(:,i+n) = mu - sqrt(lambda + n) * Sigma_sqrt(:,i-1);
        weights(i) = 1 / (2 * (lambda + n));
        weights(i+n) = 1 / (2 * (lambda + n));
    end
end

