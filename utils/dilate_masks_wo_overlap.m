function masks_dilated = dilate_masks_wo_overlap(masks_to_dilate, dilate_size)
%DILATE_MASKS_WO_OVERLAP Used in "load_dilated_masks.m". It dilates the
%masks while avoiding overlapping among the masks.
%   Usage:
%   >> masks_to_dilate = {LCmask, RCmask, LPmask, RPmask, LGPmask, RGPmask,LTHmask, RTHmask};
%   >> masks_dilated = dilate_masks_wo_overlap(masks_to_dilate, 1);

% calculate the mask of all the regions
mask_all_regions = zeros(size(masks_to_dilate{1}));
for i = 1:length(masks_to_dilate)
    mask_all_regions = mask_all_regions + double(masks_to_dilate{i} > 0.5);
end

% calculate the positions of pixels outside the masks
k = find(mask_all_regions == 0);
[r,c,s] = ind2sub(size(mask_all_regions), k);
non_mask_pixel_pos = [c,r,s];

% calculate the positions of pixels inside the masks
for i = 1:length(masks_to_dilate)
    k = find(masks_to_dilate{i} > 0.5);
    [r,c,s] = ind2sub(size(masks_to_dilate{i}), k);
    mask_pixel_poses{i} = [c,r,s];
end

% calculate the distance of pixels outside the masks to the pixels inside
for i = 1:length(mask_pixel_poses)
    d = pdist2(mask_pixel_poses{i}, non_mask_pixel_pos, 'Euclidean', 'smallest', 1);
    D(:,i) = d';
end

[M,ind] = min(D,[],2);

% separate the whole brain pixels into different regions
for i = 1:length(mask_pixel_poses)
    pos = non_mask_pixel_pos(ind == i,:);
    k = sub2ind(size(mask_all_regions), pos(:,2), pos(:,1), pos(:,3));
    mask_tmp = zeros(size(mask_all_regions));
    mask_tmp(k) = 1;
    regions_for_masks{i} = double(mask_tmp > 0.5 | masks_to_dilate{i} > 0.5);
end

if 0
    fh = figure;
    hold on;
    plot_mask(regions_for_masks{1}, [1 1 0], 0);
    plot_mask(regions_for_masks{2}, [0.8 0.8 0], 0);
    plot_mask(regions_for_masks{3}, [0 1 0], 0);
    plot_mask(regions_for_masks{4}, [0 0.8 0], 1);
end

%% dilate the masks inside the regions for each mask
for i = 1:length(masks_to_dilate)
    mask_i = dilate_mask(masks_to_dilate{i}, dilate_size);
    masks_dilated{i} = mask_i .* regions_for_masks{i};
end


end

