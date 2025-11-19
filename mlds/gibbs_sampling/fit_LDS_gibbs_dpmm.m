function [model, r_ik, extra] = fit_LDS_gibbs_dpmm(X, ts, K, options)
%FIT_LDS_GIBBS Summary of this function goes here
%   Detailed explanation goes here
t_start = tic;

alpha = parse_param(options,'alpha',1);
num_iter = parse_param(options,'num_iter',1000);
% predict_method = parse_param(options,'predict_method','ConjugatePrior');
predict_params = parse_param(options,'predict_params',[]);

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
    clusters(k) = ClusterConjugatePrior(inds, X(inds), ts(inds), ...
        predict_params);
end
% cluster(K+1) handles calculation of new cluster
predict_params.dim = D;
clusters(K+1) = ClusterConjugatePrior([], {}, {}, predict_params);

zs_all = zeros(num_iter, N);
As_all = cell(1, num_iter);
sigma2s_all = cell(1, num_iter);
pis_all = cell(1, num_iter);

nll = [];
for iter = 1:num_iter
    %% sample z
    rand_inds = randperm(N);
    for i = 1:N
        ind = rand_inds(i);
        curr_k = z(ind);
        
        clusters(curr_k).removeData(ind);
        
        % remove cluster if empty
        if isempty(clusters(curr_k).X)
            delete(clusters(curr_k));
            clusters(curr_k) = [];
            mus_first(curr_k,:) = [];
            Sigmas_first(:,:,curr_k) = [];
            K = K-1;
            inds_reduce = find(z>curr_k);
            z(inds_reduce) = z(inds_reduce) - 1;
        end
        
        % calculate the posterior of z_i given the rest zs
        p = log([get_cluster_sizes(clusters(1:end-1)), alpha]);
        for k = 1:K+1
            p(k) = p(k) + clusters(k).calcLogPredictive(X(ind), ts(ind));
        end
        p = exp(p - max(p));
        p = p / sum(p);

        % sample z_i
        u = rand;
        new_k = 1+sum(u>cumsum(p));
        
        % initialize new cluster
        if new_k == K+1
            K = K + 1;
            clusters(K+1) = ClusterConjugatePrior([], {}, {}, predict_params);
            mus_first = cat(1, mus_first, mu_first);
            Sigmas_first = cat(3, Sigmas_first, Sigma_first);
        end
        
        % add data back
        z(ind) = new_k;
        clusters(new_k).addData(ind, X(ind), ts(ind));
    end
    % store generated sample
    zs_all(iter, :) = z;
    
    %% sample As, sigmas, pis
    for k = 1:K
        [A, sigma2] = clusters(k).sampleAsAndSigmas();
        As_all{iter}(:,:,k) = A;
        sigma2s_all{iter}(k) = sigma2;
    end
    
    pis_all{iter} = sample_dirichlet(alpha/K + get_cluster_sizes(clusters(1:end-1)));
    
    %% calculate nll
    nll = [nll, calc_complete_data_nll(X, ts, As_all{iter}, ...
        sqrt(sigma2s_all{iter}), mus_first, Sigmas_first, ...
        pis_all{iter}, zs_all(iter,:))];
    
    if mod(iter, round(num_iter/10)) == 0
        disp(['Processed iteration ',num2str(iter)]);
    end
end

burn_in_ratio = 0.4;
keep_ind = ceil(num_iter*burn_in_ratio)+1;

r_ik = calc_prob_z(zs_all(keep_ind:end,:));

model = [];
model.As = mean(cat(4, As_all{keep_ind:end}), 4);
model.sigmas = mean(sqrt(cell2mat(sigma2s_all(keep_ind:end)')), 1);
model.pis = mean(cell2mat(pis_all(keep_ind:end)'), 1);
model.mus_first = mus_first;
model.Sigmas_first = Sigmas_first;

extra = [];
extra.nll = nll;

t_elapsed = toc(t_start);
disp(['Elapsed time is ',num2str(t_elapsed),' seconds']);

end


function nums = get_cluster_sizes(clusters)
for k = 1:length(clusters)
    nums(k) = length(clusters(k).X);
end
    
end