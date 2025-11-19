function [A_ks,ind] = init_As(X, V, delta_tau, K, options)
if nargin < 5
    options = [];
end

A_all = estimate_candidate_As(X, delta_tau);

% A_std = nanstd(A_all); % 1 by D
% 
% % randomly sample A from the mean and std. It does not work well.
% for k = 1:K
%     tmp = randn(1,D);
%     tmp = tmp .* A_std + mean_A;
%     A_ks(:,:,k) = diag(tmp);
% end

A_all_wo_nan = A_all;
A_all_wo_nan(any(isnan(A_all), 2),:) = [];

init_example_mode = parse_param(options, 'init_example_mode', 'random');
if strcmp(init_example_mode, 'random')
    sel_inds = randperm(size(A_all_wo_nan,1), K); 
    C = A_all_wo_nan(sel_inds,:);
    ind = knnsearch(C, A_all);
elseif strcmp(init_example_mode, 'kmeans')
    [ind, C] = kmeans(A_all_wo_nan, K);
    ind = knnsearch(C, A_all);
end

for k = 1:K
    A_ks(:,:,k) = diag(C(k,:));
end

init_As_random_noise = parse_param(options, 'init_As_random_noise', 0);
if init_As_random_noise
    A_ks = add_random_noise_to_A(A_ks);
end

end

function A_all = estimate_candidate_As(X, delta_tau)
% restore delta_t to original t
max_taus = [];
for i = 1:length(delta_tau)
    tau{i} = cumsum([0,delta_tau{i}]);
    max_taus(end+1) = tau{i}(end);
end

tq = (0:floor(median(max_taus)));

% sample to get x at the same time 
D = size(X{1}, 1);
interp_x = []; % N by T by D
for i = 1:length(X)
    resampled = interp1(tau{i}', X{i}', tq', 'linear', NaN);
    resampled = reshape(resampled, [1 length(tq) D]); 
    interp_x = cat(1, interp_x, resampled);
end

% take the median curve
% mean_x = reshape(nanmedian(interp_x, 1), [length(tq) D]); % T by D
% figure, plot(mean_x);

% since delta_t from tq is always 1, the speed will be just the distance
% mean_A = (mean_x(2:end,:) - mean_x(1:end-1,:)) ./ mean_x(1:end-1,:);
% mean_A = mean(mean_A, 1); % 1 by D

A_all = (interp_x(:,2:end,:) - interp_x(:,1:end-1,:)) ./ interp_x(:,1:end-1,:);
A_all = squeeze(nanmean(A_all, 2)); % N by D
end

function A_ks = add_random_noise_to_A(A_ks)
[rows,cols,K] = size(A_ks);
for k = 1:K
    Ak_max = max(max(abs(A_ks(:,:,k))));
    noise = randn(rows,cols) * (Ak_max/10);
    A_ks(:,:,k) = A_ks(:,:,k) + noise;
end

end