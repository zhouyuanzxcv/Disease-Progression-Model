function log_X_theta = calc_logprob_X_given_theta_t(X, ts, As, sigmas, nus, mus, Sigmas)
% Calculate logarithm of the time series with t-distributed noise residues
% given different cluster parameters
% Input:
%   X - cell array of N sequences of data, each with size D x Ti
%   ts - cell array of N sequences of time, each with size D x Ti
%   As - transition matrix of K clusters, size: D x D x K
%   sigmas - noise standard deviation for K clusters, size: 1 x K
%   nus - degree of freedom for K clusters, size: 1 x K
%   mus - mean of first data point for K clusters, size: K x D
%   Sigmas - covariance of first data point for K clusters, size: D x D x K
% Output:
%   log_X_theta - log of probability of each series for each cluster, size:
%                 N x K
[V, delta_tau] = calculate_velocity(X, ts);

D = size(X{1},1);
K = size(As, 3);
log_X_theta = zeros(length(X), K);

for k = 1:K
    for i = 1:length(X)
        Ti = size(X{i},2);
        
        % for the first time point, it is assumed to be generated from a
        % Gaussian distribution
        tmp = logmvn(X{i}(:,1)', mus(k,:), Sigmas(:,:,k));
        if Ti == 1     
            tmp2 = tmp;
        else
            % for the remaining time points, they are from t-distributions
            % with diagonal scale matrices.
            tmp = tmp + (Ti - 1) * gammaln((nus(k) + D) / 2);
            tmp = tmp - (Ti - 1) * gammaln(nus(k)/2);
            tmp = tmp - D * sum(log(delta_tau{i}*sigmas(k)));
            tmp = tmp - (D/2) * (Ti - 1) * log(nus(k) * pi);
            
            % Delta t^2 from the scale matrix is merged into the numerater
            % to create the V.
            tmp1 = sum((V{i} - As(:,:,k)*X{i}(:,1:end-1)).^2, 1) ...
                / (nus(k)*sigmas(k)^2) + 1;
            tmp2 = tmp - (nus(k) + D)/2 * sum(log(tmp1));
        end
        log_X_theta(i,k) = tmp2;
    end
end

end