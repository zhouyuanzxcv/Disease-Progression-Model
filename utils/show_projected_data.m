function proj_data = show_projected_data(train_data, train_ts, Us, class_inds)
%SHOW_PROJECTED_DATA Summary of this function goes here
%   Detailed explanation goes here
if isempty(Us)
    K = max(class_inds);
    for k = 1:K
        Us{k} = eye(size(train_data{1}, 1));
    end
end

K = length(Us);
M = size(Us{1},2);

figure;
for k = 1:K
    train_data_k = train_data(class_inds == k);
    train_ts_k = train_ts(class_inds == k);


    projs = {};
    for n = 1:length(train_data_k)
        for j = 1:M
            u1 = Us{k}(:,j);
            if u1(1) < 0, u1 = -u1; end

            proj1 = u1' * train_data_k{n};
            t = train_ts_k{n}/365;

            subplot(M, K, k + (j-1)*K);
            plot(t, proj1); hold on
            xlabel({'Year',['Cluster ',num2str(k)]});
            ylabel(['Projected to u_',num2str(j)]);

            projs{n}(j,:) = proj1;
        end
    end

    proj_data(class_inds == k) = projs;
    
end
% subject_ids = train_ids;

% save('eigenvector_projected_data.mat','proj_data', ...
%     'proj_ts','subject_ids','first_dates');
% saveR('eigenvector_projected_data.Rdata','proj_data_1','proj_data_4', ...
%     'proj_data_ts','subject_ids','first_dates');
end

