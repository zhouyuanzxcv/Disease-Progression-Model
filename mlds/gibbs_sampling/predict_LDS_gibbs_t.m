function X_pred = predict_LDS_gibbs_t(X, ts, start_ind, model, extra)
%PREDICT_LDS Summary of this function goes here
%   Detailed explanation goes here
if start_ind < 2
    error('Prediction should start at least at the second scan');
end

As = model.As;
sigmas = model.sigmas;
nus = model.nus;
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
X_pred = X;
for i = 1:length(X)
    X_pred{i}(:,start_ind:end) = 0;
end

if 0
    As_all = As;
    sigmas_all = sigmas;
    nus_all = nus;
    pis_all = pis;
    cluster_sizes_all = extra.cluster_sizes_all(end,:);
else
    As_all = extra.As_all(:,:,:,extra.keep_ind:end);
    sigmas_all = sqrt(extra.sigma2s_all(extra.keep_ind:end,:));
    nus_all = extra.nus_all(extra.keep_ind:end,:);
    pis_all = extra.pis_all(extra.keep_ind:end,:);
    cluster_sizes_all = extra.cluster_sizes_all(extra.keep_ind:end,:);
end

if size(nus_all,2) == 1
    nus_all = repmat(nus_all, [1 K]);
end

for iter = 1:size(pis_all,1)
    cluster_sizes = cluster_sizes_all(iter, :);
    As = As_all(:, :, cluster_sizes~=0, iter);
    sigmas = sigmas_all(iter, cluster_sizes~=0);
    nus = nus_all(iter, cluster_sizes~=0);
    pis = pis_all(iter, cluster_sizes~=0);
    K1 = size(As,3);
    
    for t = start_ind:T_max
        % select the indices that have length >= t
        inds = find(Ts >= t);

        % calculate the posterior
        X1 = X(inds);
        ts1 = ts(inds);
        for i = 1:length(inds)
            X1{i} = X1{i}(:,1:t-1);
            ts1{i} = ts1{i}(:,1:t-1);
        end
        log_X_theta = calc_logprob_X_given_theta_t(X1, ts1, As, ...
            sigmas, nus, mus_first, Sigmas_first);
        r_sk = calc_gamma_E_step(log_X_theta, pis);

        % predict
        movement = zeros(length(inds), K1, D);
        for k = 1:K1
            for i = 1:length(inds)
                 movement(i,k,:) = (eye(D) + delta_tau{inds(i)}(t-1) * As(:,:,k)) * X1{i}(:,t-1);
            end
        end
        X_new = sum(repmat(r_sk, [1 1 D]) .* movement, 2);

        % assign to X_pred
        for i = 1:length(inds)
            X_pred{inds(i)}(:,t) = X_pred{inds(i)}(:,t) + squeeze(X_new(i,:,:));
        end
    end
end

for i = 1:length(X_pred)
    X_pred{i}(:,start_ind:end) = X_pred{i}(:,start_ind:end) / size(pis_all,1);
end

end

