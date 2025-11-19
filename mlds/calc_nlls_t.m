function nll = calc_nlls_t(train_data, train_ts, model, extra)
K = max(extra.zs_all(:));
num_iter = size(extra.zs_all,1);

mus_first = model.mus_first;
Sigmas_first = model.Sigmas_first;

As_all = extra.As_all;
sigma2s_all = extra.sigma2s_all;
nus_all = extra.nus_all;
pis_all = extra.pis_all;
cluster_sizes_all = extra.cluster_sizes_all;

nll = []; 

for iter = 1:num_iter
    nll(iter) = calc_neg_loglikelihood_t(train_data, train_ts, As_all(:,:,:,iter), ...
        sqrt(sigma2s_all(iter,:)), nus_all(iter,:), mus_first, ...
        Sigmas_first, pis_all(iter,:), cluster_sizes_all(iter,:));
end

end


