function [outputArg1,outputArg2] = create_synthetic_data(inputArg1,inputArg2)
%CREATE_SYNTHETIC_DATA Summary of this function goes here
%   Detailed explanation goes here

% result_file = 'analysis/result_processed_images_multiple_Gibbs_finite_t_3.mat';

result_file = ['analysis/result_processed_images_flipped_shifted_dilate_1',...
    '_multiple_Gibbs_finite_t_3.mat'];


load(result_file);
K = 3;

dataset = load_dataset(dataset_name, train_inds, test_inds, []);
assert(all(strcmp(labels, dataset.labels)));

train_data = dataset.train_data;     
train_ts = dataset.train_ts;         
train_ids = dataset.train_ids;       
test_data = dataset.test_data;       
test_ts = dataset.test_ts;           
test_ids = dataset.test_ids;         
train_inds = dataset.train_inds;     
test_inds = dataset.test_inds;       
labels = dataset.labels;       
D = size(train_data{1},1);

all_elapsed_days = train_ts;
subject_ids = train_ids;


for noise_level = (3:-1:1)
    model1 = model;
    % noise_level = 3 means small noise, 1 means large noise
    model1.sigmas = model1.sigmas * 0.1^(noise_level-1);
    model1.pis = ones(1,K)*(1/K);
    all_image_data = {};
    class_inds = [];
    
    for i = 1:length(train_data)
        all_image_data{i} = train_data{i}(:,1);
        
        class_idx = find(mnrnd(1,model1.pis));
        class_inds(i) = class_idx;
        
        for j = 1:length(train_ts{i}) - 1
            delta_t = (all_elapsed_days{i}(j+1) - all_elapsed_days{i}(j))/365;
%             noise = randn(D,1) * model1.sigmas(class_idx) * delta_t;
            noise = trnd(model1.nus, [D,1]) * model1.sigmas(class_idx) * delta_t;
            A = (eye(D) + delta_t * model1.As(:,:,class_idx)); 
            all_image_data{i}(:,j+1) = A * all_image_data{i}(:,j) + noise; 
        end
    end
    save(['synthetic_noise_',num2str(noise_level),'.mat'], 'class_inds', ...
        'all_image_data', 'all_elapsed_days', 'subject_ids', 'labels');
end

end

