function [imgs_rescaled, scales, scale2s, normalizing_voxels] = ...
    rescale_by_occ_lobe(imgs, occMask, mode)
%RESCALE_BY_OCC_LOBE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    mode = 'median';
end

nImgs=size(imgs,4);
occIndex=find(occMask>0.5);

imgs_rescaled = zeros(size(imgs));
scales = zeros(1, nImgs);
scale2s = zeros(1, nImgs);

for i=1:nImgs,
    tmp=squeeze(imgs(:,:,:,i));
    occ_values = tmp(occIndex);
    normalizing_voxels_i = occ_values;
    
    occMed = median(occ_values);
    occMean = mean(occ_values);
    
    if strcmp(mode, 'median')
        tmp = tmp / occMed;
        scales(i) = occMed;
        scale2s(i) = occMean;
    elseif strcmp(mode, 'mean')
        tmp = tmp / occMean;
        scales(i) = occMean;
        scale2s(i) = occMed;
    end
    
    imgs_rescaled(:,:,:,i) = tmp;
    normalizing_voxels(:,i) = normalizing_voxels_i;
    
end

end