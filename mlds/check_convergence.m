function [outputArg1,outputArg2] = check_convergence(extra)
%CHECK_CONVERGENCE Summary of this function goes here
%   Detailed explanation goes here
num_runs = 5;

if nargin < 1
    K = 4;
    result_file = ['analysis/result_processed_images_Gibbs_finite_t_',num2str(K)];
    
    if 0
        parfor run = 1:num_runs
            options = [];
            options.num_iter = 1500;
            options.K = K;

            result_file1 = [result_file,'_',num2str(run)];
            options.result_filename = result_file1;
            run_algo(options);
        end
    end
        
    load([result_file,'_1']);

    dataset = load_dataset(dataset_name, train_inds, test_inds, []);
    assert(all(strcmp(labels, dataset.labels)));
end

train_data = dataset.train_data;     
train_ts = dataset.train_ts;         
train_ids = dataset.train_ids;   

% test_data = dataset.test_data;       
% test_ts = dataset.test_ts;           
% test_ids = dataset.test_ids;         
% train_inds = dataset.train_inds;     
% test_inds = dataset.test_inds;       
% labels = dataset.labels;       

% nll1 = calc_nll(train_data, train_ts, extra);
% nll2 = calc_nll_complete_data(train_data, train_ts, extra);
% 
% figure,semilogx(nll1);
% figure,semilogx(nll2);



parfor run = 1:num_runs
    result_file1 = [result_file, '_', num2str(run)];
    S = load(result_file1);
    extra = S.extra; 
    model = S.model;
    
    nll(:,run) = calc_nlls_t(train_data, train_ts, model, extra);
end

figure;
for run = 1:num_runs
    subplot(2,3,run);
    semilogx(nll(:,run));
end

for run = 1:num_runs
    result_file1 = [result_file, '_', num2str(run)];
    load(result_file1);
    [M,ind] = max(r_sk, [], 2);

    show_seq_quiver_plot(train_data, train_ids, labels, struct('classes',ind));
end

chosen = find_optimal_converged_chains(nll);
chosen

[num_attempts, ratios, means] = assess_mixing(nll(:,chosen));

figure,plot(num_attempts, ratios);


end


function nll = calc_nll_complete_data(train_data, train_ts, extra)
K = max(extra.zs_all(:));
num_iter = size(extra.zs_all,1);

[mu_first, Sigma_first] = calc_first_point_statistics(train_data);
mus_first = repmat(mu_first, [K 1]);
Sigmas_first = repmat(Sigma_first, [1 1 K]);

As_all = extra.As_all;
sigma2s_all = extra.sigma2s_all;
nus_all = extra.nus_all;
pis_all = extra.pis_all;
zs_all = extra.zs_all;

nll = []; % nll

for iter = 1:num_iter
    nll(iter) = calc_complete_data_nll_t(train_data, train_ts, As_all(:,:,:,iter), ...
        sqrt(sigma2s_all(iter,:)), nus_all(iter,:), mus_first, Sigmas_first, ...
        pis_all(iter,:), zs_all(iter,:));
end

end


