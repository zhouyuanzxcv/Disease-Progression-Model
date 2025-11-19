function log_X_theta = calc_logprob_X_given_theta(X, ts, As, sigmas, mus, Sigmas)
% Calculate logarithm of Gaussian given different cluster parameters
% Input:
%   X - cell array of N sequences of data, each with size D x Ti
%   ts - cell array of N sequences of time, each with size D x Ti
%   As - transition matrix of K clusters, size: D x D x K
%   sigmas - noise standard deviation for K clusters, size: 1 x K
%   mus - mean of first data point for K clusters, size: K x D
%   Sigmas - covariance of first data point for K clusters, size: D x D x K
% Output:
%   log_X_theta - log of probability of each sequence for each cluster, size:
%                 N x K
[V, delta_tau] = calculate_velocity(X, ts);

D = size(X{1},1);
K = size(As, 3);
log_X_theta = zeros(length(X), K);

for k = 1:K
    for i = 1:length(X)
        num_scans = size(X{i},2);
        
        % implementation 1
%         mus1 = zeros(num_scans, D);
%         Sigmas1 = zeros(D, D, num_scans);
%         mus1(1,:) = mus(k,:);
%         Sigmas1(:,:,1) = Sigmas(:,:,k);
%         
%         tmp = As(:,:,k) * X{i}(:,1:end-1);
%         tmp = tmp .* repmat(delta_tau{i}, [D 1]);
%         tmp = tmp + X{i}(:,1:end-1);
%         mus1(2:end,:) = tmp';
%         
%         tmp = repmat(sigmas(k)^2 * eye(D), [1 1 num_scans-1]);
%         tmp2 = repmat(reshape(delta_tau{i}.^2, [1 1 num_scans-1]), ...
%             [size(tmp,1), size(tmp,2)]);
%         Sigmas1(:,:,2:end) = tmp .* tmp2;
%         
%         tmp1 = sum(logmvn(X{i}', mus1, Sigmas1));
        
        % implementation 2
        tmp = logmvn(X{i}(:,1)', mus(k,:), Sigmas(:,:,k));
        if num_scans == 1
            tmp2 = tmp;
        else
            tmp = tmp - (num_scans - 1) * D * (0.5*log(2*pi)+log(sigmas(k)));
            tmp = tmp - D * sum(log(delta_tau{i}));
            tmp2 = tmp - 1/ (2*sigmas(k)^2) * sum(sum((V{i} - As(:,:,k)*X{i}(:,1:end-1)).^2));
        end
%         
%         assert(abs(tmp1 - tmp2) < 1e-9*abs(tmp1));
        
        % implementation 2 is faster than 1 with 0.64 secs vs. 1.25 secs
        log_X_theta(i,k) = tmp2;
    end
end

end