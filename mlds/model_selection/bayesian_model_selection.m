function K_optimal = bayesian_model_selection
%BAYESIAN_MODEL_SELECTION Summary of this function goes here
%   Detailed explanation goes here
Kmax = 10;
logIs = [];

method = 'Gibbs_finite';
dataset_name = 'processed_images';

if 0
    for K = 1:4
%         try
            options = [];
            options.method = method;
            options.dataset_name = dataset_name;
            options.K = K;
            
            run_algo(options);
            
            close all
%         catch me
%             disp(me.message)
%         end
    end
end

if 0
    for K = 1:Kmax
        file_name = ['analysis/result_',dataset_name,'_',method,'_',num2str(K),'.mat'];
        if ~exist(file_name, 'file')
            logIs(K) = nan;
            continue;
        end
        load(file_name);
        dataset = load_dataset(dataset_name, train_inds, test_inds, []);
        
        train_data = dataset.train_data;
        train_ts = dataset.train_ts;
        
        % filter out the empty clusters
        disp(['Cluster sizes for ',file_name,' : ', ...
            num2str(extra.cluster_sizes_all(end,:))]);
        
        inds_nonempty = (extra.cluster_sizes_all(end,:) ~= 0);
        As_all = extra.As_all(:, :, inds_nonempty, extra.keep_ind:end);
        sigma2s_all = extra.sigma2s_all(extra.keep_ind:end, inds_nonempty);
        pis_all = extra.pis_all(extra.keep_ind:end, inds_nonempty);
        mus_first = model.mus_first(inds_nonempty,:);
        Sigmas_first = model.Sigmas_first(:,:,inds_nonempty);
        
        nlls = [];
        for s = 1:size(pis_all,1)
            nll = calc_neg_loglikelihood(train_data, train_ts, ...
                As_all(:,:,:,s), sqrt(sigma2s_all(s,:)), mus_first, ...
                Sigmas_first, pis_all(s,:));
            nlls(s) = nll;
        end
        C = -min(nlls);
        logIs(K) = log(length(nlls)) + C - log(sum(exp(nlls + C)));
        final_cluster_sizes(K) = length(find(extra.cluster_sizes_all(end,:) ~= 0));
    end
    save('analysis/bayesian_model_selection_result.mat','logIs','final_cluster_sizes'); 
else
    load('analysis/bayesian_model_selection_result.mat');
end

%% show log-likelihood
figure,
yyaxis left
plot(logIs, 'linewidth', 1);
xlabel('Number of clusters (K)');
ylabel('Log-likelihood');

yyaxis right
plot(final_cluster_sizes, 'linewidth', 1);
ylabel('Number of nonempty clusters');


[~,K_optimal] = max(logIs);
K_optimal

end

