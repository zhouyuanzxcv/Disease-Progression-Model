function [model, r_sk, extra] = merge_models_t(models, r_sks, extras)
%MERGE_MODELS_T Summary of this function goes here
%   Detailed explanation goes here
num_runs = length(models);

fprintf('The proportions of clusters of the %d runs are: \n', num_runs);
all_pis = [];
for i = 1:num_runs
    all_pis(i,:) = models{i}.pis;
end
all_pis

for run = 1:num_runs
    nll(:,run) = extras{run}.nll;
end

if 0
    figure;
    colors = distinguishable_colors(num_runs);
    
    for i = 1:num_runs
        semilogx(nll(:,i), 'color', colors(i,:)); 
        hold on;
    end
    xlabel('Iteration number'); ylabel('NLL');
end

chosen = find_optimal_converged_chains(nll);

if isempty(chosen)
    model = [];
    r_sk = [];
    extra = [];
    return;
end

% if 1

idx = 1;

% The following code is added after running the algorithm on my 
% personal computer and finding that it tends to consider all the 
% chains converge by only analyzing NLL. The possible reason is that
% the originally processed PPMI data lead to different clustering results
% with different NLL but the flipped (dilated) data produce different clustering
% results with very similar NLL. However, on the desktop computer, it still tends
% to correctly exclude 1 different result when 4 chains converge to the same
% result. To make the code behave similarly across different machines and 
% different datasets, the following code is added after NLL-based combination
% analysis.

if length(chosen) > 2
    similarity_mat = zeros(length(chosen), length(chosen));
    for i = 1:length(chosen)
        for j = 1:length(chosen)
            [~, ci] = max(r_sks{chosen(i)}, [], 2);
            [~, cj] = max(r_sks{chosen(j)}, [], 2);
            [AR,RI,MI,HI] = rand_index(ci,cj);
            similarity_mat(i,j) = RI;
        end
    end
    similarity_chosen = sum(similarity_mat, 2);
    [~, idx] = max(similarity_chosen, [], 1);
end

% choose one of them
fprintf('The chosen chain is No. %d \n', chosen(idx));
model = models{chosen(idx)};
r_sk = r_sks{chosen(idx)};
extra = extras{chosen(idx)};
    
% else   
    % merge different results
%     [model, r_sk, extra] = merge_chosen(chosen, models, r_sks, extras);


end

% function [model, r_sk, extra] = merge_chosen(chosen, models, r_sks, extras)
% template_ind = chosen(1);
% [D,D,K] = size(models{template_ind}.As);
% template_As = reshape(models{template_ind}.As, [D*D K])';
% 
% As_all_merged = extras{template_ind}.As_all;
% sigma2s_all_merged = extras{template_ind}.sigma2s_all;
% nus_all_merged = extras{template_ind}.nus_all;
% pis_all_merged = extras{template_ind}.pis_all;
% cluster_sizes_all_merged = extras{template_ind}.cluster_sizes_all;
% zs_all_merged = extras{template_ind}.zs_all;
% 
% for model_ind = 2:length(chosen)
%     extra = extras{chosen(model_ind)};
%     model = models{chosen(model_ind)};
%     As1 = reshape(model.As, [D*D K])';
%     
%     % how to handle nan in As is not yet solved in the permute function
%     [P, As2] = permute_endmembers(template_As, As1);
%     
%     nus_all = extra.nus_all(keep_ind:end,:);
%     sigma2s_all = extra.sigma2s_all(keep_ind:end,:);
%     pis_all = extra.pis_all(keep_ind:end,:);
%     cluster_sizes_all = extra.cluster_sizes_all(keep_ind:end,:);
%     zs_all = extra.zs_all(keep_ind:end,:);
%     As_all = extra.As_all(:,:,:,keep_ind:end);
%     
%     nus_all_merged = cat(1, nus_all_merged, nus_all);
%     sigma2s_all_merged = cat(1, sigma2s_all_merged, sigma2s_all * P');
%     pis_all_merged = cat(1, pis_all_merged, pis_all * P');
%     cluster_sizes_all_merged = cat(1, cluster_sizes_all_merged, ...
%         cluster_sizes_all * P');
%     
%     zs_all1 = one_hot_encoding(zs_all, 'to_one_hot');
% 
%     for i = 1:size(pis_all,1)
%         As1 = reshape((P * reshape(extra.As_all(:,:,:,i), [D*D K])')', [D D K]);
%         As_all_merged(:,:,:,end+1) = As1;
%         zs1 = squeeze(zs_all1(i,:,:)) * P';
%         zs1 = one_hot_encoding(reshape(zs1, [1 size(zs1,1) size(zs1,2)]), 'from_one_hot');
%         zs_all_merged(end+1,:) = zs1;
%     end
% end
% 
% extra.keep_ind = 1;
% extra.As_all = As_all_merged;
% extra.sigma2s_all = sigma2s_all_merged;
% extra.nus_all = nus_all_merged;
% extra.pis_all = pis_all_merged;
% extra.cluster_sizes_all = cluster_sizes_all_merged;
% extra.zs_all = zs_all_merged;
% 
% mus_first = models{template_ind}.mus_first;
% Sigmas_first = models{template_ind}.Sigmas_first;
% [model, r_sk] = calc_model(extra, mus_first, Sigmas_first);
% 
% end
% 
% function [model, r_sk] = calc_model(extra, mus_first, Sigmas_first)
% keep_ind = extra.keep_ind;
% r_ik = calc_prob_z(extra.zs_all(keep_ind:end,:));
% 
% model = [];
% model.As = mean(extra.As_all(:,:,:,keep_ind:end), 4);
% model.sigmas = mean(sqrt(extra.sigma2s_all(keep_ind:end,:)), 1);
% model.nus = mean(extra.nus_all(keep_ind:end,:), 1);
% model.pis = mean(extra.pis_all(keep_ind:end,:), 1);
% model.mus_first = mus_first;
% model.Sigmas_first = Sigmas_first;
% 
% end

