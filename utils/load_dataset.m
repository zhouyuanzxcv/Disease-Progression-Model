function dataset = load_dataset(dataset_name, train_inds, test_inds, options)
%LOAD_DATASET Summary of this function goes here
%   Detailed explanation goes here
dataset = [];

% preprocess_file = 'normalization_params_median_mu.mat';
preprocess_file = [];
options = insert_param_when_absent(options, 'remove_scan_outliers', 0);
options = insert_param_when_absent(options, 'use_minus_one', 1);
options = insert_param_when_absent(options, 'use_region_mean', 1);


switch dataset_name
    case {'processed_images','processed_images_flipped',...
            'processed_images_flipped_shifted', ...
            'processed_images_flipped_shifted_dilate_1', ...
            'processed_images_flipped_shifted_dilate_1_2021', ...
            'processed_images_flipped_shifted_dilate_2', ...
            'processed_images_flipped_shifted_dilate_3', ...
            'processed_images_using_masks_flipped_registered', ...
            'processed_images_using_masks_flipped_registered_dilate_1'}
        if 0
            processed_data_dir = ['../data', dataset_name(length('processed_images')+1:end)];

            [all_image_data, all_elapsed_days, subject_ids, labels, extra] = ...
                load_processed_data(processed_data_dir, preprocess_file, options);
    %         for i = 1:length(all_image_data)
    %             all_image_data{i}(2,:) = all_image_data{i}(2,:) - 0.1;
    %             all_image_data{i}(4,:) = all_image_data{i}(4,:) + 0.1;
    %         end
            
            save(['../mlds/data/', dataset_name, '_data.mat'],'all_image_data', ...
                'all_elapsed_days','subject_ids','labels','extra');
        else
            load(['../mlds/data/', dataset_name, '_data.mat']);
        end
        dataset.first_dates = extra.first_dates;
        dataset.ages = extra.ages;
        dataset.genders = extra.genders; % male is 1, female is 0
        
        % check gender consistency
        sexes = {};
        for i = 1:length(extra.all_flags)
            sexes{i,1} = extra.all_flags{i}(end).sex;
            sexes{i,2} = dataset.genders(i);
        end        
    case {'processed_images_dilate_1', 'processed_images_dilate_2', ...
            'processed_images_erode_1', 'processed_images_erode_2', ...
            }
        type = dataset_name;
        type([1:17]) = [];
        savefile_name = ['./data/', 'processed_images_data_', type,'.mat'];
        if 1
            processed_data_dir = ['../data_',type];

            [all_image_data, all_elapsed_days, subject_ids, labels, extra] = ...
                load_processed_data(processed_data_dir, preprocess_file, options);
            
            save(savefile_name,'all_image_data', ...
                'all_elapsed_days','subject_ids','labels','extra');
        else
            load(savefile_name);
        end
        dataset.first_dates = extra.first_dates;
        dataset.ages = extra.ages;
        dataset.genders = extra.genders;
    case 'csv_file'
        [all_image_data, all_elapsed_days, subject_ids, labels] = load_data_from_csv();
    case 'csv_file_363'
        [all_image_data, all_elapsed_days, subject_ids, labels] = load_data_from_csv();
        S = load('./data/processed_images_data.mat');
        ids_selected = S.subject_ids;
        inds = [];
        for i = 1:length(ids_selected)
            ind = find(strcmp(ids_selected{i}, subject_ids));
            if ~isempty(ind)
                inds(end+1) = ind;
            end
        end
        all_image_data = all_image_data(inds);
        all_elapsed_days = all_elapsed_days(inds);
        subject_ids = subject_ids(inds);
    otherwise
        S = load(dataset_name);
        all_image_data = S.all_image_data;
        all_elapsed_days = S.all_elapsed_days;
        subject_ids = S.subject_ids;
        labels = S.labels;
        S = rmfield(S, {'all_image_data','all_elapsed_days','subject_ids','labels'});
        dataset.extra = S;
        if isfield(S, 'class_inds')
            dataset.class_inds = S.class_inds;
        end
        if isfield(S, 'first_dates')
            dataset.first_dates = S.first_dates;
        end
        if isfield(S, 'extra')
            if isfield(S.extra, 'ages')
                dataset.ages = S.extra.ages;
            end
            if isfield(S.extra, 'genders')
                dataset.genders = S.extra.genders;
            end
            if isfield(S.extra, 'first_dates')
                dataset.first_dates = S.extra.first_dates;
            end
        end
end

% change (LC, RC, LP, RP) to (LC, LP, RP, RC)
% reorder = [1 3 4 2]; 
if ~isempty(labels) && strcmp(labels{1},'LC')
    [labels, reorder] = reorder_labels(labels, {'LC','LP','RP','RC'});
    for i = 1:length(all_image_data)
        all_image_data{i} = all_image_data{i}(reorder,:);
    end
end

%% split into training and testing
N = length(all_image_data);
[train_inds1, test_inds1] = random_split_indices(N, 0.8);
if isempty(train_inds)
    train_inds = 1:N;
end
if isempty(test_inds)
    test_inds = test_inds1;
end

train_data = all_image_data(train_inds);
train_ts = all_elapsed_days(train_inds);
train_ids = subject_ids(train_inds);

test_data = all_image_data(test_inds);
test_ts = all_elapsed_days(test_inds);
test_ids = subject_ids(test_inds);

% show_seq_quiver_plot(train_data(1:5:end), train_ids(1:5:end), labels);

dataset.train_data = train_data;
dataset.train_ts = train_ts;
dataset.train_ids = train_ids;
dataset.test_data = test_data;
dataset.test_ts = test_ts;
dataset.test_ids = test_ids;

dataset.train_inds = train_inds;
dataset.test_inds = test_inds;

dataset.labels = labels;

end

