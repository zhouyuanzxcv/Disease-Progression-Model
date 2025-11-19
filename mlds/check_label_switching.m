function perc = check_label_switching(extra)
%CHECK_LABEL_SWITCHING Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    result_file = 'analysis/result_processed_images_Gibbs_finite_t_3.mat';
    load(result_file);

    dataset = load_dataset(dataset_name, train_inds, test_inds, []);
    assert(all(strcmp(labels, dataset.labels)));

    % train_data = dataset.train_data;     
    % train_ts = dataset.train_ts;         
    % train_ids = dataset.train_ids;       
    % test_data = dataset.test_data;       
    % test_ts = dataset.test_ts;           
    % test_ids = dataset.test_ids;         
    % train_inds = dataset.train_inds;     
    % test_inds = dataset.test_inds;       
    % labels = dataset.labels;      
end

N = size(extra.zs_all, 2);

zs = extra.zs_all(extra.keep_ind:end,:);
nll = extra.nll(extra.keep_ind:end);

num_iter = size(zs,1);
K = max(zs(:));

% change to one-of-K encoding
zs1 = zeros(num_iter,N,K);

for k = 1:K
    zs1(:,:,k) = double(zs == k);
end

% take the z closest to the MAP as template
[min_nll,min_ind] = min(nll);
z_template = squeeze(zs1(min_ind,:,:));

% find best permutation
ps = [];
for iter = 1:num_iter
    z = squeeze(zs1(iter,:,:)); % N x K
    p = find_best_permutation(z_template, z, K);
    ps = [ps;p];
end

% count how many are identity permutation
perc = length(find(all(ps == repmat((1:K),[num_iter 1]), 2))) / num_iter;

end

function p = find_best_permutation(z_template, z, K)
P1 = perms([1:K]);
P1 = flipud(P1); % make identity permutation the first

for i = 1:size(P1,1)
    err(i) = sum(sum(abs(z(:,P1(i,:)) - z_template).^2));
end

[~,ind] = min(err);
p = P1(ind,:);

end