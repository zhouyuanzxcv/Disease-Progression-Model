function [model, r_sk, extra] = fit_LDS(X, ts, K, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
t_start = tic;

use_centrosym = parse_param(options, 'use_centrosym', 1);
use_single_Sigma_for_first_point = parse_param(options, ...
    'use_single_Sigma_for_first_point', 1);
options.use_single_Sigma_for_first_point = use_single_Sigma_for_first_point;

[V, delta_tau] = calculate_velocity(X, ts);

[As, sigmas, pis, mus, Sigmas] = init_LDS(X, ts, K, options);
N = length(X);
D = size(X{1},1);
Ts = zeros(N,1);
for i = 1:length(X)
    Ts(i) = size(X{i}, 2);
end

VX1 = zeros(D,D,N);
XX1 = zeros(D,D,N);
for i = 1:length(X)
    VX1(:,:,i) = V{i} * X{i}(:,1:end-1)';
    XX1(:,:,i) = X{i}(:,1:end-1) * X{i}(:,1:end-1)';
end

X1s = zeros(N,D);
for i = 1:length(X)
    X1s(i,:) = X{i}(:,1)';
end

nll = calc_neg_loglikelihood(X, ts, As, sigmas, mus, Sigmas, pis);

for iter = 1:200
    % E step
    log_X_theta = calc_logprob_X_given_theta(X, ts, As, ...
        sigmas, mus, Sigmas);
    r_sk = calc_gamma_E_step(log_X_theta, pis);
    
    % M step
    % update A
    pis = sum(r_sk,1) / N;
    for k = 1:K
        tmp = repmat(reshape(r_sk(:,k), [1 1 N]), [D D]);
        F = sum(tmp .* VX1, 3);
        E = sum(tmp .* XX1, 3);
        if ~use_centrosym
            As(:,:,k) = F * inv(E);
        else
            As(:,:,k) = fmin_centrosym(E,F);
        end
    end
    % update sigma
    for k = 1:K
        tmp = zeros(N,1);
        for i = 1:N
            tmp(i) = sum(sum((V{i} - As(:,:,k) * X{i}(:,1:end-1)).^2, 1), 2);
        end
        tmp = sum(tmp .* r_sk(:,k));
        tmp = tmp / sum(r_sk(:,k) .* (Ts-1) * D);
        sigmas(k) = sqrt(tmp);
    end
    
    % update mus and Sigmas
    if ~use_single_Sigma_for_first_point
        for k = 1:K
            mus(k,:) = sum(X1s .* repmat(r_sk(:,k), [1 D])) / sum(r_sk(:,k));
            tmp = X1s - repmat(mus(k,:), [N 1]);
            Sigmas(:,:,k) = tmp' * (repmat(r_sk(:,k), [1 D]) .* tmp) / sum(r_sk(:,k));
        end
    end
    
    % calculate likelihood value
    nll(end+1) = calc_neg_loglikelihood(X, ts, As, sigmas, mus, Sigmas, pis);
    
    % converge
    if test_convergence(nll, 1e-6, 2)
        break;
    end
end

t_elapsed = toc(t_start);
disp(['Elapsed time is ',num2str(t_elapsed),' seconds']);

extra = [];
extra.nll = nll;

model = [];
model.As = As;
model.sigmas = sigmas;
model.pis = pis;
model.mus_first = mus;
model.Sigmas_first = Sigmas;

end


