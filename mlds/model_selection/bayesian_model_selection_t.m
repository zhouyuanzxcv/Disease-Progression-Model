function K_optimal = bayesian_model_selection_t(dataset_name)
%BAYESIAN_MODEL_SELECTION Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
%     dataset_name = 'processed_images';
    dataset_name = 'processed_images_flipped_shifted_dilate_1_2021';
    % dataset_name = 'csv_file_363';
    % dataset_name = 'csv_file';
end

Kmax = 6;

logIs = [];
entropies = [];

method = 'multiple_Gibbs_finite_t';

savefile = ['analysis/',dataset_name,'_bayesian_model_selection_result_t.mat'];

if 1
    parfor K = 1:Kmax
%         try
            options = [];
            options.method = method;
            options.dataset_name = dataset_name;
            options.K = K;
            options.result_filename = ['analysis/result_',dataset_name,...
                '_',method,'_bayesian_model_selection_',num2str(K),'.mat'];
            
            run_algo(options);
            
            close all
%         catch me
%             disp(me.message)
%         end
    end
end

if 1
    for K = 1:Kmax
        file_name = ['analysis/result_',dataset_name,'_',method,...
            '_bayesian_model_selection_',num2str(K),'.mat'];
        if ~exist(file_name, 'file')
            logIs(K) = nan;
            continue;
        end
        load(file_name);
        dataset = load_dataset(dataset_name, train_inds, test_inds, []);
        
        train_data = dataset.train_data;
        train_ts = dataset.train_ts;
        
        % filter out the empty clusters
        disp(['Cluster sizes at the last iteration for ',file_name,' : ', ...
            num2str(extra.cluster_sizes_all(end,:))]);
        
        nlls = calc_nlls_t(train_data, train_ts, model, extra);
        nlls = nlls(extra.keep_ind:end);
        
        C = -min(nlls);
        logIs(K) = log(length(nlls)) + C - log(sum(exp(nlls + C)));
        final_cluster_sizes(K) = length(find(extra.cluster_sizes_all(end,:) ~= 0));
        
        entropies(K) = calc_entropy(r_sk);

    end
    save(savefile,'logIs','final_cluster_sizes'); 
else
    load(savefile);
end

%% show log-likelihood
fh = figure;
yyaxis left
plot(logIs, '-s', 'linewidth', 1);
xlabel('Number of subtypes (K)');
ylabel('Log-likelihood');

yyaxis right
plot(final_cluster_sizes, '--*', 'linewidth', 1);
ylabel('Number of nonempty subtypes');


[~,K_optimal] = max(logIs);
K_optimal

set(fh,'position',[680 741 308 237]);

end

function e = calc_entropy(ps)
for i = 1:size(ps,1)
    p = ps(i,:);
    p(p == 0) = [];
    es(i) = sum(-p.*log(p));
end
e = sum(es);

end