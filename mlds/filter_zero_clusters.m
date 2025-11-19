function model = filter_zero_clusters(model, extra)
%FILTER_ZERO_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here
szs = extra.cluster_sizes_all(end,:);
assert(length(szs) == size(model.As,3));

model.As = model.As(:,:,szs ~= 0);
model.sigmas = model.sigmas(szs ~= 0);
model.pis = model.pis(szs ~= 0);
model.mus_first = model.mus_first(szs ~= 0,:);
model.Sigmas_first = model.Sigmas_first(:,:,szs ~= 0);

end

