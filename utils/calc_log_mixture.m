function value = calc_log_mixture(N_nk, w)
% calculate the logarithm of a mixture model
% Input:
%   N_nk : N by K matrix containing log(p(x_n|theta_k))
%   w : 1 by K vector containing pis.
% Output:
%   value : N by 1 vector of log(\sum_k{\pi_k*p(x_n|theta_k)})
[N,K] = size(N_nk);
% Need to avoid too large negative logarithm
N_nk_max = max(N_nk,[],2);

% eval_N1 = -sum(log(sum((ones(N,1)*w_k) .* exp(N_nk), 2)));
eval_N = log(sum(repmat(w,N,1) .* exp(N_nk - repmat(N_nk_max,1,K)), 2));
value = N_nk_max + eval_N;

end
