function test_r_sk = predict_subtype(test_data, test_ts, pis, ...
    As, sigmas, nus, mus_first, Sigmas_first)
%PREDICT_SUBTYPE Summary of this function goes here
%   Detailed explanation goes here
K = size(As,3);
if size(nus, 2) == 1
    nus = repmat(nus, [1 K]);
end
log_X_theta = calc_logprob_X_given_theta_t(test_data, test_ts, As, ...
        sigmas, nus, mus_first, Sigmas_first);
test_r_sk = calc_gamma_E_step(log_X_theta, pis);
end

