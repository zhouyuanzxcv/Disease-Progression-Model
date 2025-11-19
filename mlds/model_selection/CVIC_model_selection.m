function K_optimal = CVIC_model_selection()
%ESTIMATE_NUM_COMP Summary of this function goes here
%   Detailed explanation goes here


Kmax = 10;

method = 'Gibbs_finite';
% method = 'EM';
dataset_name = 'processed_images';

dataset = load_dataset(dataset_name, [], [], []);
train_data = dataset.train_data;     
N = length(train_data);
get_number_samples_with_length(train_data, dataset.train_ts);

num_folds = 10;
if 0
    Ks = (3:4);
    
    inds = crossvalind('Kfold',N,num_folds);
    K_and_folds = cartprod(Ks, (1:num_folds));
    Ks_cart = K_and_folds(:,1);
    folds_cart = K_and_folds(:,2);
    parfor tuple_idx = 1:size(K_and_folds,1)
       try
            K = Ks_cart(tuple_idx);
            j = folds_cart(tuple_idx);
            
            test_inds = find(inds == j);
            train_inds = find(inds ~= j);
            file_name = get_filename(dataset_name, method, K, j);
            
            options = [];
            options.method = method;
            options.dataset_name = dataset_name;
            options.K = K;
            options.train_inds = train_inds;
            options.test_inds = test_inds;
            options.result_filename = file_name;
            
            run_algo(options);
            close all
        catch me
            disp(me.message)
        end
    end
end

nll = [];
num_clusters = [];
for K = 1:Kmax
    for j = 1:num_folds
        file_name = get_filename(dataset_name, method, K, j);
        if ~exist(file_name, 'file')
            nll(K,j) = nan;
            continue;
        end
        load(file_name);
        dataset = load_dataset(dataset_name, train_inds, test_inds, []);

        test_data = dataset.test_data;     
        test_ts = dataset.test_ts;         
        
        disp(['Cluster sizes for ',file_name,' : ', ...
            num2str(extra.cluster_sizes_all(end,:))]);
        num_clusters(K,j) = length(find(extra.cluster_sizes_all(end,:) ~= 0));
        nll(K,j) = calc_likelihood(test_data, test_ts, model, extra);
    end
end

nll1 = nansum(nll,2);

figure,

yyaxis left
plot(-nll1, 'linewidth', 1);
xlabel('Number of clusters (K)');
ylabel('Log-likelihood');

yyaxis right
plot(mean(num_clusters, 2), 'linewidth', 1);
ylabel('Number of nonempty clusters in average');


[~,K_optimal] = min(nll1);
K_optimal


end

function val = calc_likelihood(data, ts, model, extra)
szs = extra.cluster_sizes_all(end,:);
As = model.As(:,:,szs ~= 0);
sigmas = model.sigmas(szs ~= 0);
pis = model.pis(szs ~= 0);
mus_first = model.mus_first(szs ~= 0,:);
Sigmas_first = model.Sigmas_first(:,:,szs ~= 0);

val = calc_neg_loglikelihood(data, ts, As, sigmas, mus_first, ...
    Sigmas_first, pis);

end

function val = calc_marginal_likelihood(data, ts, model, extra)
As_all = extra.As_all(:,:,:,extra.keep_ind:end);
sigma2s_all = extra.sigma2s_all(extra.keep_ind:end,:);
pis_all = extra.pis_all(extra.keep_ind:end,:);

nlls = [];
for s = 1:size(pis_all,1)
    nll = calc_neg_loglikelihood(data, ts, ...
        As_all(:,:,:,s), sqrt(sigma2s_all(s,:)), model.mus_first, ...
        model.Sigmas_first, pis_all(s,:));
    nlls(s) = nll;
end
C = -min(nlls);
val = log(length(nlls)) + C - log(sum(exp(nlls + C)));

end

function filename1 = add_fold_postfix(filename, j)
filename1 = [filename(1:end-4),'_CV_',num2str(j),'.mat'];
end

function filename = get_filename(dataset_name, method, K, j)
filename = ['analysis/result_',dataset_name,'_',method,'_CV_',num2str(K),...
    '_',num2str(j),'.mat'];
end

function get_number_samples_with_length(X, ts)
Ts = zeros(1, length(X));
for i = 1:length(X)
    Ts(i) = size(X{i},2);
end
T_max = max(Ts);
length(find(Ts == 2))
length(find(Ts == 3))
length(find(Ts == 4))
length(find(Ts == 5))
length(find(Ts == 6))

[V, delta_tau] = calculate_velocity(X, ts);
delta_tau_all = cell2mat(delta_tau);
length(find(delta_tau_all > 0.5 & delta_tau_all < 1.5))
length(find(delta_tau_all > 1.5 & delta_tau_all < 2.5))
length(find(delta_tau_all > 2.5 & delta_tau_all < 3.5))
length(find(delta_tau_all > 3.5 & delta_tau_all < 4.5))

end
