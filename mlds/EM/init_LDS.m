function [A_ks, sigmas, pis, mus, Sigmas] = init_LDS(X, ts, K, options)
[V, delta_tau] = calculate_velocity(X, ts);

[A_ks,ind] = init_As(X, V, delta_tau, K, options);

% initialize noise std
errors = [];
for k = 1:K
    for i = 1:length(X)
        errors(i,k) = mean(mean((V{i} - A_ks(:,:,k)*X{i}(:,1:end-1)).^2));
    end
end

for k = 1:K
    sigmas(k) = sqrt(mean(errors(ind == k, k)));
end

% when ind has a cluster with only one sample, sigmas will be 0 for that
% sample
if any(sigmas == 0)
    error('Initialization has one cluster with too few samples');
end

% other parameters
pis = ones(1,K) * (1/K);

if options.use_single_Sigma_for_first_point
    [mu, Sigma] = calc_first_point_statistics(X);
    mus = repmat(mu, [K 1]);
    Sigmas = repmat(Sigma, [1 1 K]);
else
    [mus, Sigmas] = init_first_scan_params(X, K);
end

end

function [mus, Sigmas] = init_first_scan_params(X, K)
first_scans = [];
for i = 1:length(X)
    first_scans = cat(1, first_scans, X{i}(:,1)');
end
mus = mean(first_scans);
Sigmas = cov(first_scans);
mus = repmat(mus, [K 1]);
Sigmas = repmat(Sigmas, [1 1 K]);

end

