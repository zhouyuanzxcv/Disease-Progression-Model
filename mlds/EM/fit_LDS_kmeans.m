function [model, r_ik, extra] = fit_LDS_kmeans(X, ts, K, options)
%FIT_LDS_KMEANS Summary of this function goes here
%   Detailed explanation goes here
[V, delta_tau] = calculate_velocity(X, ts);
for i = 1:length(X)
    X{i}(:,end) = [];
    TiMinus1(i) = size(V{i},2);
end
D = size(X{1},1);
N = length(X);

num_iter = parse_param(options,'num_iter',300);

if 1
    A = fit_A(V, X);
    error = calc_errors(V,X,A);
    [class_inds, c] = kmeans(error',K,'Start','cluster');
else
    class_inds = ceil(rand(N,1)*K);
end

As = [];
As_old = zeros(D,D,K);

As_diff = [];
for iter = 1:num_iter
    % Update A
    for k = 1:K
        V1 = V(class_inds == k);
        X1 = X(class_inds == k);
        As(:,:,k) = fit_A(V1,X1);
    end
    
    % Update idx
    error = [];
    for k = 1:K
        error(:,k) = calc_errors(V,X,As(:,:,k));
    end
    [~,class_inds] = min(error, [], 2);
    
    % check convergence
    As_diff(iter) = max(abs(As_old(:) - As(:)));
    if As_diff(iter) < 1e-4
        disp(['Converged at iteration ',num2str(iter)]);
        break;
    end
    As_old = As;
end

r_ik = one_hot_encoding(class_inds', 'to_one_hot');
r_ik = squeeze(r_ik);

pis = [];
sigma2s = [];
for k = 1:K
    inds = (class_inds == k);
    pis(k) = length(find(inds));
    
    T1 = TiMinus1(inds);
    tmp = T1 .* error(inds, k)';
    sigma2s(k) = sum(tmp) / sum(T1) / D;
end
pis = pis / N;

model = [];
model.As = As;
model.sigmas = sqrt(sigma2s);
model.pis = pis;

[mu_first, Sigma_first] = calc_first_point_statistics(X);
mus_first = repmat(mu_first, [K 1]);
Sigmas_first = repmat(Sigma_first, [1 1 K]);

model.mus_first = mus_first;
model.Sigmas_first = Sigmas_first;

extra = [];
extra.nll = As_diff;

end

function A = fit_A(V,X)
V1 = cell2mat(V);
X1 = cell2mat(X);
A = V1*X1'*inv(X1*X1');

end

function err = calc_errors(V, X, A)
for i = 1:length(V)
    err(i) = mean(sum((V{i} - A * X{i}).^2, 1));
end
end