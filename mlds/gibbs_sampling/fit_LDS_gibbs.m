function [model, r_ik, extra] = fit_LDS_gibbs(X, ts, K, options)
%FIT_LDS_GIBBS Summary of this function goes here
%   Detailed explanation goes here
t_start = tic;

alpha = parse_param(options,'alpha',1);
num_iter = parse_param(options,'num_iter',3000);
predict_method = parse_param(options,'predict_method','ConjugatePrior');
predict_params = parse_param(options,'predict_params',[]);
burn_in_ratio = parse_param(options,'burn_in_ratio',0.4);

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

for k = 1:K
    inds = find(z==k);
    switch predict_method
        case 'DeltaMAP'
            clusters(k) = ClusterDeltaMAP(inds, X(inds), ts(inds), ...
                predict_params);
        case 'ConjugatePrior'
            clusters(k) = ClusterConjugatePrior(inds, X(inds), ts(inds), ...
                predict_params);
        otherwise
    end
end

zs_all = zeros(num_iter, N);
As_all = zeros(D, D, K, num_iter);
sigma2s_all = zeros(num_iter, K);
pis_all = zeros(num_iter, K);
cluster_sizes_all = zeros(num_iter, K);

nll = [];
for iter = 1:num_iter
    %% sample z
    rand_inds = randperm(N);
    for i = 1:N
        ind = rand_inds(i);
        curr_k = z(ind);
        
        clusters(curr_k).removeData(ind);
        
        % calculate the posterior of z_i given the rest zs
        p = log(alpha/K + get_cluster_sizes(clusters));
        for k = 1:K
            p(k) = p(k) + clusters(k).calcLogPredictive(X(ind), ts(ind));
        end
        p = exp(p - max(p));
        p = p / sum(p);

        % sample z_i
        u = rand;
        new_k = 1+sum(u>cumsum(p));
        
        % add data back
        z(ind) = new_k;
        clusters(new_k).addData(ind, X(ind), ts(ind));
    end
    % store generated sample
    zs_all(iter, :) = z;
    
    %% sample As, sigmas, pis
    for k = 1:length(clusters)
        [A, sigma2] = clusters(k).sampleAsAndSigmas();
        As_all(:,:,k,iter) = A;
        sigma2s_all(iter,k) = sigma2;
    end
    
    pis_all(iter,:) = sample_dirichlet(alpha/K + get_cluster_sizes(clusters));
    cluster_sizes_all(iter,:) = get_cluster_sizes(clusters);
    
    %% calculate nll
%     nll = [nll, calc_complete_data_nll(X, ts, As_all(:,:,:,iter), ...
%         sqrt(sigma2s_all(iter,:)), mus_first, Sigmas_first, ...
%         pis_all(iter,:), zs_all(iter,:))];
    nll(iter) = calc_neg_loglikelihood(X, ts, ...
        As_all(:,:,:,iter), sqrt(sigma2s_all(iter,:)), mus_first, ...
        Sigmas_first, pis_all(iter,:));
    
    if mod(iter, round(num_iter/10)) == 0
        disp(['Processed iteration ',num2str(iter)]);
    end
end

keep_ind = ceil(num_iter*burn_in_ratio)+1;

r_ik = calc_prob_z(zs_all(keep_ind:end,:));

model = [];
model.As = mean(As_all(:,:,:,keep_ind:end), 4);
model.sigmas = mean(sqrt(sigma2s_all(keep_ind:end,:)), 1);
model.pis = mean(pis_all(keep_ind:end,:), 1);
model.mus_first = mus_first;
model.Sigmas_first = Sigmas_first;

extra = [];
extra.nll = nll;
extra.As_all = As_all;
extra.sigma2s_all = sigma2s_all;
extra.pis_all = pis_all;
extra.zs_all = zs_all;
extra.cluster_sizes_all = cluster_sizes_all;
extra.keep_ind = keep_ind;

t_elapsed = toc(t_start);
disp(['Final cluster sizes are: ',num2str(cluster_sizes_all(end,:))]);
disp(['Elapsed time is ',num2str(t_elapsed),' seconds']);

end


function nums = get_cluster_sizes(clusters)
for k = 1:length(clusters)
    nums(k) = length(clusters(k).X);
end
    
end