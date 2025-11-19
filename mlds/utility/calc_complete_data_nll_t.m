function val = calc_complete_data_nll_t(X, ts, As, sigmas, nus, mus, Sigmas, pis, z)
%CALC_COMPLETE_DATA_NLL Summary of this function goes here
%   Detailed explanation goes here
K = size(As,3);

if length(nus) == 1
    nus = repmat(nus, [1 K]);
end

log_X_theta = calc_logprob_X_given_theta_t(X, ts, As, ...
    sigmas, nus, mus, Sigmas);
for k = 1:K
    inds = (z == k);
    ll(inds) = log_X_theta(inds,k) + log(pis(k));
end
val = -sum(ll); % negative log-likelihood
end

