function [all_image_data, all_elapsed_days, subject_ids, labels, extra] = ...
    load_processed_data(processed_data_dir, preprocess_file, options)
%LOAD_PROCESSED_DATA Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    options = [];
end

remove_scan_outliers = parse_param(options, 'remove_scan_outliers', 0);
use_minus_one = parse_param(options, 'use_minus_one', 0);
use_region_mean = parse_param(options,'use_region_mean',0);

show_only_positive_images = 0;
show_only_negative_images = 0;
use_only_first_scan = 0;

if ~isempty(preprocess_file)
    load(preprocess_file);
end

filelist = dir(processed_data_dir);
is_dir = [filelist.isdir];
filelist(is_dir) = [];

num_subjects = length(filelist);
all_elapsed_days = {};
all_image_data = {};
subject_ids = {};
all_flags = {};

list_positive_subjects = [];
list_positive_seq_id = [];
list_negative_subjects = [];
list_negative_seq_id = [];

load(['../masks/masks_for_',processed_data_dir(4:end)]);
occMask = mask.Occ;
LSmask = mask.LS;
RSmask = mask.RS;
LCmask = mask.LC;
RCmask = mask.RC;
LPmask = mask.LP;
RPmask = mask.RP;

[region_sizes,labels] = get_region_num_of_voxels(mask);

for i = 1:num_subjects
    subject_id = filelist(i).name(1:end-4);
    load(fullfile(processed_data_dir, [subject_id, '.mat']));
    
    % store subject id and seq. id with positive or negative values
    neg_seq_id = find([flags.isNegative] == 1);
    pos_seq_id = find([flags.isNegative] ~= 1);
    list_positive_subjects = [list_positive_subjects, ...
        repmat(str2double(subject_id), [1 length(pos_seq_id)])];
    list_positive_seq_id = [list_positive_seq_id, pos_seq_id];
    list_negative_subjects = [list_negative_subjects, ...
        repmat(str2double(subject_id), [1 length(neg_seq_id)])];
    list_negative_seq_id = [list_negative_seq_id, neg_seq_id];
    
    % filter the images based on some criterion
    if show_only_positive_images 
        neg_idx = find([flags.isNegative] == 1);
        elapsedDays(neg_idx) = [];
        imgData(:,neg_idx) = [];
    end
    if show_only_negative_images
        neg_idx = find([flags.isNegative] == 1);
        elapsedDays = elapsedDays(neg_idx);
        imgData = imgData(:,neg_idx);
    end
    if use_only_first_scan
        elapsedDays = elapsedDays(1);
        imgData = imgData(:,1);
    end
    
    % save data
    first_dates(i) = firstDate;
    all_elapsed_days{i} = elapsedDays;
    if use_minus_one
        imgData(1:6,:) = imgData(1:6,:) - region_sizes(1:6)' * ones(1,size(imgData,2));
    end
    if use_region_mean
        imgData(1:6,:) = imgData(1:6,:) ./ repmat(region_sizes(1:6)', [1,size(imgData,2)]);        
    end
    if ~isempty(preprocess_file) 
        imgData(1:6,:) = imgData(1:6,:) .* repmat(imgData(7,:), [6 1]);
        
        imgData = imgData .* repmat(normalization_params{i}(1,:), [7 1]);
        
        intercepts = normalization_params{i}(2,:);
        imgData(7,:) = imgData(7,:) + intercepts;
        
        num_voxels = [length(find(LCmask>0.5)), length(find(RCmask>0.5)),...
            length(find(LPmask>0.5)), length(find(RPmask>0.5)), ...
            length(find(LSmask>0.5)), length(find(RSmask>0.5))];
        offsets = num_voxels' * intercepts;
        imgData(1:6,:) = imgData(1:6,:) + offsets;
        
        imgData(1:6,:) = imgData(1:6,:) ./ repmat(imgData(7,:), [6 1]);
    end
    all_image_data{i} = imgData;
    subject_ids{i} = subject_id;
    all_flags{i} = flags;
end

extra.list_positive_subjects = list_positive_subjects;
extra.list_positive_seq_id = list_positive_seq_id;
extra.list_negative_subjects = list_negative_subjects;
extra.list_negative_seq_id = list_negative_seq_id;
extra.first_dates = first_dates;
extra.all_flags = all_flags;
extra.mask = mask;

if remove_scan_outliers
    [all_image_data, all_elapsed_days, subject_ids] = ...
        remove_scan_outliers_internal(all_image_data, ...
        all_elapsed_days, subject_ids);
end

% manually remove outliers. 
% after recreating the data with previous data removed, these two outliers
% disappeared
if 0 % removed 3360 and 3605
inds = [];
for i = 1:length(all_image_data)
    img_data = all_image_data{i}(1:4,:);
    if max(img_data(:)) > 4
        inds = [inds, i];
        subject_ids{i}
    end
end
all_image_data(inds) = [];
all_elapsed_days(inds) = [];
subject_ids(inds) = [];
end

S = load('../outliers/compiledAgeGender.mat');
for i = 1:length(subject_ids)
    id = str2num(subject_ids{i});
    ind = find(S.subjectId == id);
    ages(i) = S.subjectAge(ind);
    genders(i) = S.subjectGender(ind);
end
extra.ages = ages;
extra.genders = genders;

end

function [all_image_data1, all_elapsed_days1, subject_ids1] = ...
    remove_scan_outliers_internal(all_image_data, all_elapsed_days, subject_ids)
fcn_is_greater_than = @(x1,x2) all(x1 >= x2);

is_valid = false(1, length(all_image_data));

for i = 1:length(all_image_data)
    data = all_image_data{i}(1:4,:);
    [len,inds] = find_longest_decreasing_subsequence(data, fcn_is_greater_than);
    all_image_data1{i} = all_image_data{i}(:,inds);
    all_elapsed_days1{i} = all_elapsed_days{i}(:,inds);

    if len > 1
        is_valid(i) = true;
    end
end

all_image_data1 = all_image_data1(is_valid);
all_elapsed_days1 = all_elapsed_days1(is_valid);
subject_ids1 = subject_ids(is_valid);

disp(['Original number of scans: ', num2str(get_total_num_scans(all_image_data))]);
disp(['Number of scans after removing outliers: ', num2str( ...
    get_total_num_scans(all_image_data1))]);

end

function total_size = get_total_num_scans(all_image_data)
all_sizes = [];
for i = 1:length(all_image_data)
    all_sizes(i) = size(all_image_data{i},2);
end
total_size = sum(all_sizes);
end

function [max_len,max_inds] = find_longest_decreasing_subsequence( ...
    data, fcn_is_greater_than)
% use dynamic programming to find the LDS
len = size(data, 2);

inds = zeros(1,len); % stores the index linked to the previous 
max_lens = zeros(1,len); % stores the maximal length of subsequences that end at k

for k = 1:len
    previous = [];
    previous_inds = [];
    for j = 1:k-1
        if fcn_is_greater_than(data(:,j), data(:,k)) 
            previous = [previous, max_lens(j)];
            previous_inds = [previous_inds, j];
        end
    end
    if isempty(previous)
        max_lens(k) = 1;
        inds(k) = -1;
    else
        [M,I] = max(previous);
        max_lens(k) = M + 1;
        inds(k) = previous_inds(I);
    end
end


[max_len,I] = max(max_lens);
max_inds = I;
while inds(I) ~= -1
    I = inds(I);
    max_inds = [I, max_inds];
end
    
end