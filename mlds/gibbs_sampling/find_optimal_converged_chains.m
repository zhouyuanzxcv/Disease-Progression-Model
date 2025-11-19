function [converged_optimal] = find_optimal_converged_chains(nll)
%FIND_OPTIMAL_CONVERGED_CHAINS Summary of this function goes here
%   Detailed explanation goes here
num_runs = size(nll,2);

converged = {};
last_means = [];
is_found = 0;
for k = (num_runs:-1:2)
    chosen = nchoosek((1:num_runs), k);
    for r = 1:size(chosen,1)
        [num_attempts, ratios, means] = assess_mixing(nll(:,chosen(r,:)));
        if ratios(end) < 1.03
            converged(end+1) = {chosen(r,:)};
            last_means(end+1) = means(end);
            if k > num_runs/2 % more than half converged
                is_found = 1;
                break;
            end
        end
    end
    if is_found
        break;
    end
end

if length(converged) > 1 % more than 1 group with converged distributions
%     [~,ind] = min(last_means); % for nll, the less, the better
    lens = cellfun(@length, converged);
    [~,ind] = max(lens);
    if length(find(lens == lens(ind))) > 1 % more than 1 group with the same maximal size
        converged_optimal = [];
    else
        converged_optimal = converged{ind};
    end
elseif length(converged) == 1
    converged_optimal = converged{1};
else
    converged_optimal = [];
end

end

