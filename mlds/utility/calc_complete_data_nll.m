function val = calc_complete_data_nll(X, ts, As, sigmas, mus, Sigmas, pis, z)
%CALC_COMPLETE_DATA_NLL Summary of this function goes here
%   Detailed explanation goes here
log_X_theta = calc_logprob_X_given_theta(X, ts, As, ...
    sigmas, mus, Sigmas);
K = size(As,3);
for k = 1:K
    inds = (z == k);
    ll(inds) = log_X_theta(inds,k) + log(pis(k));
end
val = -sum(ll); % negative log-likelihood
end

