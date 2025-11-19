function [model, r_ik, extra] = fit_LDS_gibbs_t(X, ts, K, options)
%FIT_LDS_GIBBS Summary of this function goes here
%   Detailed explanation goes here
t_start = tic;

alpha = parse_param(options,'alpha',K);
num_iter = parse_param(options,'num_iter',1000);
% predict_method = parse_param(options,'predict_method','ConjugatePrior');
predict_params = parse_param(options,'predict_params',[]);
burn_in_ratio = parse_param(options,'burn_in_ratio',0.5);
use_multiple_nu = parse_param(options,'use_multiple_nu',0);
save_weights_all = parse_param(options,'save_weights_all',0);

[V, delta_tau] = calculate_velocity(X, ts);
[mu_first, Sigma_first] = calc_first_point_statistics(X);
predict_params.mu_first = mu_first;
predict_params.Sigma_first = Sigma_first;
mus_first = repmat(mu_first, [K 1]);
Sigmas_first = repmat(Sigma_first, [1 1 K]);

N = length(X);
D = size(X{1}, 1);

% initialize component assignment
z = ceil(rand(1,N)*K);
weights = {};
for i = 1:length(X)
    weights{i} = ones(1, size(X{i},2) - 1);
end

ksai0 = 1e-3;
tau0 = 1e-3;

nu = 30;

for k = 1:K
    inds = find(z==k);
    
    predict_params.weights = weights(inds);
    predict_params.nu = nu;
    clusters(k) = ClusterConjugatePrior_t(inds, X(inds), ts(inds), ...
        predict_params);
end

zs_all = zeros(num_iter, N);

As_all = zeros(D, D, K, num_iter);
As_mean_all = zeros(D, D, K, num_iter);

sigma2s_all = zeros(num_iter, K);
sigma2s_mean_all = zeros(num_iter, K);

pis_all = zeros(num_iter, K);
cluster_sizes_all = zeros(num_iter, K);
weights_all = cell(num_iter, N);

if use_multiple_nu
    nus_all = zeros(num_iter, K);
else
    nus_all = zeros(num_iter, 1);
end

nll = [];
for iter = 1:num_iter
    %% sample z
    rand_inds = randperm(N);
    for i = 1:N
        ind = rand_inds(i);
        curr_k = z(ind);
        
        w = clusters(curr_k).removeData(ind);
        
        % calculate the posterior of z_i given the rest zs
        p = log(alpha/K + get_cluster_sizes(clusters));
        if use_multiple_nu
            for k = 1:K
                p(k) = p(k) + sum(log_gamma(w{1}, clusters(k).nu/2, clusters(k).nu/2));
            end
        end
        
        for k = 1:K
            p(k) = p(k) + clusters(k).calcLogPredictive(X(ind), ts(ind), w);
        end
        p = exp(p - max(p));
        p = p / sum(p);

        % sample z_i
        u = rand;
        new_k = 1+sum(u>cumsum(p));
        
        % add data back
        z(ind) = new_k;
        clusters(new_k).addData(ind, X(ind), ts(ind), w);
    end
    % store generated sample
    zs_all(iter, :) = z;
    
    %% sample As, sigmas, pis
    for k = 1:length(clusters)
        [A, sigma2, A_mean, sigma2_mean] = clusters(k).sampleAsAndSigmas();
        As_all(:,:,k,iter) = A;
        sigma2s_all(iter,k) = sigma2;
        
        As_mean_all(:,:,k,iter) = A_mean;
        sigma2s_mean_all(iter,k) = sigma2_mean;
    end
    
    pis_all(iter,:) = sample_dirichlet(alpha/K + get_cluster_sizes(clusters));
    cluster_sizes_all(iter,:) = get_cluster_sizes(clusters);
    
    %% sample weights
    for k = 1:K
        weights(clusters(k).Ids) = clusters(k).sampleWeights();
    end
    weights_all(iter, :) = weights;
    
    %% sample nus
    if use_multiple_nu
        for k = 1:K
            nus_all(iter, k) = clusters(k).sampleNu();
        end
    else
        nu = sample_nu(ksai0, tau0, clusters);
        nus_all(iter) = nu;
    end
    
    %% calculate nll
%     nll = [nll, calc_complete_data_nll_t(X, ts, As_all(:,:,:,iter), ...
%         sqrt(sigma2s_all(iter,:)), nus_all(iter,:), mus_first, Sigmas_first, ...
%         pis_all(iter,:), zs_all(iter,:))];
    
    nll(iter) = calc_neg_loglikelihood_t(X, ts, As_all(:,:,:,iter), ...
        sqrt(sigma2s_all(iter,:)), nus_all(iter,:), mus_first, ...
        Sigmas_first, pis_all(iter,:), cluster_sizes_all(iter,:));
    
    if mod(iter, round(num_iter/10)) == 0
        disp(['Processed iteration ',num2str(iter)]);
    end
end

keep_ind = ceil(num_iter*burn_in_ratio)+1;

r_ik = calc_prob_z(zs_all(keep_ind:end,:));

model = [];
model.As = mean(As_mean_all(:,:,:,keep_ind:end), 4);
model.sigmas = mean(sqrt(sigma2s_all(keep_ind:end,:)), 1);
model.nus = mean(nus_all(keep_ind:end,:), 1);
model.pis = mean(pis_all(keep_ind:end,:), 1);
model.mus_first = mus_first;
model.Sigmas_first = Sigmas_first;

for i = 1:size(weights_all,2)
    tmp = cell2mat(weights_all(keep_ind:end,i));
    weights{i} = mean(tmp, 1);
    weights_std{i} = std(tmp);
end

model.weights = weights;
model.weights_std = weights_std;


extra = [];
extra.nll = nll;
extra.As_all = As_all;
extra.sigma2s_all = sigma2s_all;
extra.nus_all = nus_all;
extra.pis_all = pis_all;
extra.zs_all = zs_all;
if save_weights_all
    extra.weights_all = weights_all;
end
extra.cluster_sizes_all = cluster_sizes_all;
extra.keep_ind = keep_ind;

t_elapsed = toc(t_start);
disp(['Final cluster sizes are: ',num2str(cluster_sizes_all(end,:))]);
disp(['Elapsed time is ',num2str(t_elapsed),' seconds']);

end

function nu = sample_nu(ksai0, tau0, clusters)
ksai1 = ksai0;
tau1 = tau0;
K = length(clusters);
for k = 1:K
    ksai1 = ksai1 + clusters(k).SuffStats.TiMinus1;
    for i = 1:length(clusters(k).weights)
        tau1 = tau1 + 0.5 * sum(log(clusters(k).weights{i}) - clusters(k).weights{i});
    end
end     
nu = sample_T_DOF(ksai1, tau1, 1);
for k = 1:K
    clusters(k).nu = nu;
end

end

function nums = get_cluster_sizes(clusters)
for k = 1:length(clusters)
    nums(k) = length(clusters(k).X);
end

end

