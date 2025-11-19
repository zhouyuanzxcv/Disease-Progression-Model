function value = calc_neg_loglikelihood(X, ts, As, sigmas, mus, Sigmas, pis)
% As - D x D x K
% sigmas - 1 x K 
% pis - 1 x K
% mus - K X D
% Sigmas - D x D x K
log_X_theta = calc_logprob_X_given_theta(X, ts, As, ...
    sigmas, mus, Sigmas);
tmp = calc_log_mixture(log_X_theta, pis);
value = -sum(tmp); % negative log-likelihood

end
