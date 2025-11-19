function X_pred = predict_LDS(X, ts, start_ind, model, extra)
%PREDICT_LDS Summary of this function goes here
%   Detailed explanation goes here
if start_ind < 2
    error('Prediction should start at least at the second scan');
end

if isfield(model, 'nus') && ~isempty(extra) % use t prediction
    X_pred = predict_LDS_gibbs_t(X, ts, start_ind, model, extra);
else
    if isfield(extra, 'cluster_sizes_all')
        model = filter_zero_clusters(model, extra);
    end

    X_pred = predict_LDS_EM(X, ts, start_ind, model);
end


end

function X_pred = predict_LDS_EM(X, ts, start_ind, model)
As = model.As;
sigmas = model.sigmas;
pis = model.pis;
mus_first = model.mus_first;
Sigmas_first = model.Sigmas_first;

[V, delta_tau] = calculate_velocity(X, ts);

K = size(As,3);
D = size(X{1},1);

Ts = zeros(1, length(X));
for i = 1:length(X)
    Ts(i) = size(X{i},2);
end
T_max = max(Ts);

% initialize the prediction
for i = 1:length(X)
    X_pred{i} = X{i}(:,1:start_ind-1);
end

for t = start_ind:T_max
    % select the indices that have length >= t
    inds = [];
    for i = 1:length(X)
        if size(X{i},2) >= t
            inds = [inds, i];
        end
    end
    
    % calculate the posterior
    X1 = {};
    ts1 = {};
    for i = 1:length(inds)
        X1{i} = X{inds(i)}(:,1:t-1);
        ts1{i} = ts{inds(i)}(:,1:t-1);
    end
    log_X_theta = calc_logprob_X_given_theta(X1, ts1, As, ...
        sigmas, mus_first, Sigmas_first);
    r_sk = calc_gamma_E_step(log_X_theta, pis);
    
    % predict
    movement = zeros(length(inds), K, D);
    for k = 1:K
        for i = 1:length(inds)
             movement(i,k,:) = (eye(D) + delta_tau{inds(i)}(t-1) * As(:,:,k)) * X1{i}(:,t-1);
        end
    end
    X_new = sum(repmat(r_sk, [1 1 D]) .* movement, 2);
    
    % assign to X_pred
    for i = 1:length(inds)
        X_pred{inds(i)}(:,t) = squeeze(X_new(i,:,:));
    end
    
end
end