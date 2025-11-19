function [occMask,LSmask,RSmask,LCmask1,RCmask1,LPmask1,RPmask1, ...
    WSmask, LGPmask1, RGPmask1, LTHmask1, RTHmask1, SKmask] = ...
    load_dilated_masks(dilate_size)
%LOAD_DILATED_MASKS dilate_size negative means erosion
if nargin < 1
    dilate_size = 1;
end


[occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
    LGPmask, RGPmask, LTHmask, RTHmask, SKmask] = load_masks_original();

% load('../masks/flipped_registered_masks.mat');


masks_to_dilate = {LCmask, RCmask, LPmask, RPmask, LGPmask, RGPmask, ...
    LTHmask, RTHmask};

masks_dilated = dilate_masks_wo_overlap(masks_to_dilate, dilate_size);

LCmask1 = masks_dilated{1};
RCmask1 = masks_dilated{2};
LPmask1 = masks_dilated{3};
RPmask1 = masks_dilated{4};
LGPmask1 = masks_dilated{5};
RGPmask1 = masks_dilated{6};
LTHmask1 = masks_dilated{7};
RTHmask1 = masks_dilated{8};


if 0
    fh = figure;
    hold on;
    plot_mask(LCmask1, [1 1 0], 0);
    plot_mask(RCmask1, [0.8 0.8 0], 0);
    plot_mask(LPmask1, [0 1 0], 0);
    plot_mask(RPmask1, [0 0.8 0], 1);
end

% figure, imshow(LCmask(:,:,20));
% figure, imshow(LCmask1(:,:,20));
end
