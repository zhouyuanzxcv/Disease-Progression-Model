function run_algo(options)
%RUN_ALGO Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    options = [];
    options.dataset_name = './data/demo/synthetic_noise_2.mat';


    % This is the dataset for PPMI datscans. The name is long but it does
    % reflect the procedure we used to process the image. Because the
    % original datscan is in radiological space (left side of the image is 
    % actual right side of the person), we need to flip it. Also, because
    % the image and the atlas are centered at the 45th column, we need to
    % shift the images such that it is flipped along the 45th column. The
    % atlas is dilated by 1 voxel to capture the partial volume effect.
    
    % options.dataset_name = 'processed_images_flipped_shifted_dilate_1_2021';
    
%     options.dataset_name = 'csv_file_363';
%     options.dataset_name = 'csv_file';

%     options.dataset_name = '../high_dimension/processed_images_fine_regions_cluster_2.mat';


    options.save_weights_all = 0;
    options.use_parfor = 1;
    options.K = 3;
    
    % This is the method proposed in the paper. The 'multiple' in the name
    % means we run multiple (5) chains.
    options.method = 'multiple_Gibbs_finite_t';

    % change it to zero to use a full transition matrix
    options.predict_params.use_centrosym = 1;

    % If A is not a centrosymmetric matrix, specific the basis E for A,
    % i.e. A(:) = E_for_A * a. The default basis is an identity matrix,
    % i.e. no constraint on the transition matrix.
    options.predict_params.E_for_A = [];
    
%     options.method = 'EM';
%     options.alpha = 0.1 * options.K;
end

% close all

addpath('gibbs_sampling');
addpath('utility');
addpath('analysis');
addpath('updrs');
addpath('model_selection');
addpath('EM');

addpath('../utils');

dataset_name = parse_param(options, 'dataset_name', 'processed_images'); % processed_images
K = parse_param(options,'K',3);
method = parse_param(options,'method','multiple_Gibbs_finite_t'); % Gibbs_finite
use_parfor = parse_param(options,'use_parfor',0);

train_inds = parse_param(options,'train_inds',[]);
test_inds = parse_param(options,'test_inds',[]);

[filepath,name,ext] = fileparts(dataset_name);
result_filename = parse_param(options, 'result_filename', ...
    ['results/result_',name,'_',method,'_',num2str(K),'.mat']);

dataset = load_dataset(dataset_name, train_inds, test_inds, options);

train_data = dataset.train_data;     
train_ts = dataset.train_ts;         
train_ids = dataset.train_ids;       
test_data = dataset.test_data;       
test_ts = dataset.test_ts;           
test_ids = dataset.test_ids;         
train_inds = dataset.train_inds;     
test_inds = dataset.test_inds;       
labels = dataset.labels;       

%% Fit
D = size(train_data{1},1);

start_t = tic;

switch method
    case 'Kmeans'
        [model, r_sk, extra] = fit_LDS_kmeans(train_data, train_ts, K, options);
    case 'multiple_EMs'
        models = struct('As',{},'sigmas',{},'pis',{},'mus_first',{},'Sigmas_first',{});
        options.init_As_random_noise = 1;
        run_times = 100;
        for run_id = 1:run_times
            try
                [model, r_sk, extra] = fit_LDS(train_data, train_ts, K, options);
                models(end+1) = model;
            catch err
                disp(['Iter ', num2str(run_id), ' err: ', err.message]);
            end
        end
        
        models1 = permute_mixture_models(models);
        
        As_multiple = reshape([models1.As], [D D length(models) K]);
        As_multiple = permute(As_multiple, [1 2 4 3]);
        As_multiple = reshape(As_multiple, [D*D*K length(models)])';
        [coeff,score,latent] = pca(As_multiple);
        figure, scatter3(score(:,1), score(:,2), score(:,3));
    case 'EM'
        % use EM algorithm
        [model, r_sk, extra] = fit_LDS(train_data, train_ts, K, options);
    case 'EM_wo_centrosym'
        options.use_centrosym = 0;
        [model, r_sk, extra] = fit_LDS(train_data, train_ts, K, options);
    case 'Gibbs_finite'
        % use Gibbs sampling for finite mixture
        [model, r_sk, extra] = fit_LDS_gibbs(train_data, train_ts, K, options);
    case 'Gibbs_finite_wo_centrosym'
        options.predict_params.use_centrosym = 0;
        [model, r_sk, extra] = fit_LDS_gibbs(train_data, train_ts, K, options);
    case 'Gibbs_finite_t'
        [model, r_sk, extra] = fit_LDS_gibbs_t(train_data, train_ts, K, options);
    case 'multiple_Gibbs_finite_t'
        num_runs = 5;
        options.num_iter = 1500;
        
        models = {};
        r_sks = {};
        extras = {};
        while 1
            curr_len = length(models);
            if use_parfor
                parfor run_id = 1:num_runs
                    [model, r_sk, extra] = fit_LDS_gibbs_t(train_data, train_ts, K, options);
                    models{curr_len + run_id} = model;
                    r_sks{curr_len + run_id} = r_sk;
                    extras{curr_len + run_id} = extra;
                end
            else
                for run_id = 1:num_runs
                    [model, r_sk, extra] = fit_LDS_gibbs_t(train_data, train_ts, K, options);
                    models{curr_len + run_id} = model;
                    r_sks{curr_len + run_id} = r_sk;
                    extras{curr_len + run_id} = extra;
                end
            end
            [model, r_sk, extra] = merge_models_t(models, r_sks, extras);
            if ~isempty(model)
                break;
            else
                disp('No more than 2 chains converged to the same distribution');
            end
        end
    case 'Gibbs_dpmm'
        % use dirichlet process mixture model
        K = 1;
        [model, r_sk, extra] = fit_LDS_gibbs_dpmm(train_data, train_ts, K, options);
end

fprintf('The method takes %.1f seconds\n', toc(start_t));

[M,ind] = max(r_sk, [], 2);
show_seq_quiver_plot(train_data, train_ids, labels, struct('classes',ind,...
    'all_elapsed_days',{train_ts}));
% show_seq_quiver_plot(train_data, train_ids, labels, struct('classes',ind,...
%     'labels_expected',{{'LP_f','RP_f','LP_b','RP_b'}},'all_elapsed_days',{train_ts}));

figure,semilogx(extra.nll); xlabel('Iteration number'); ylabel('NLL');

%% predict test
test_predict = predict_LDS(test_data, test_ts, 2, model, extra);
show_seq_quiver_plot(test_data, test_ids, labels, struct('predicted_data', {test_predict}));
disp('Prediction error: ');
nanmean(calc_pred_error(test_data, test_predict))

disp(['pis: ',num2str(model.pis)]);
disp(['sigmas: ',num2str(model.sigmas)]);
for k = 1:K
    [V,D] = eigs(model.As(:,:,k));
    V
    d = diag(D)'
end

if isfield(dataset, 'extra') && isfield(dataset.extra, 'class_inds')
    [~, class_inds_est] = max(r_sk, [], 2);
    ari = rand_index(dataset.extra.class_inds, class_inds_est);
    disp('ARI between the ground truth clustering and the estimated clustering')
    ari
end

save(result_filename, 'dataset_name', 'method', 'K', 'model', ...
    'r_sk', 'extra', 'train_inds', 'test_inds', 'labels');

end
