function value = calc_neg_loglikelihood_t(X, ts, As, sigmas, nus, mus, ...
    Sigmas, pis, cluster_size)
% As - D x D x K
% sigmas - 1 x K 
% nus - 1 x K or 1 x 1
% pis - 1 x K
% mus - K X D
% Sigmas - D x D x K

inds_nonempty = (cluster_size ~= 0);
As = As(:, :, inds_nonempty);
sigmas = sigmas(:, inds_nonempty);
if length(nus) ~= 1
    nus = nus(:, inds_nonempty);
end
pis = pis(:, inds_nonempty);
mus = mus(inds_nonempty,:);
Sigmas = Sigmas(:,:,inds_nonempty);


K = size(As, 3);
if length(nus) == 1
    nus = repmat(nus, [1 K]);
end

log_X_theta = calc_logprob_X_given_theta_t(X, ts, As, ...
    sigmas, nus, mus, Sigmas);
tmp = calc_log_mixture(log_X_theta, pis);
value = -sum(tmp); % negative log-likelihood

end
