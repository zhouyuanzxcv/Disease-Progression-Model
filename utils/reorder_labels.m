function [labels_new, reorder_inds] = reorder_labels(labels, labels_ordered)
%REORDER_LABELS Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   labels: cell array say {'LC','RC','LP','RP'}
%   labels_ordered: cell array say {'LC','LP','RP','RC'}
% Output:
%   labels_new: cell array say {'LC','LP,'RP','RC'}
%   reoder_inds: reoder indices say [1 3 4 2] 
for i = 1:length(labels_ordered)
    ind = find(strcmp(labels, labels_ordered{i}));
    if isempty(ind)
        ind = i;
    end
    reorder_inds(i) = ind;
    labels_new{i} = labels{ind};
end

end

