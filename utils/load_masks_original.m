function [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
    LGPmask, RGPmask, LTHmask, RTHmask, SKmask] = load_masks_original()
%RETRIEVE_META_DATA Summary of this function goes here
%   Detailed explanation goes here
load('..\outliers\metaData.mat');
origImgSize=metaData.origImgSize;
zLimU=metaData.zLim.uLim;
zLimL=metaData.zLim.lLim;
clear metaData

LCmask = niftiread('..\masks\caudate_left.nii');
RCmask = niftiread('..\masks\caudate_right.nii');
LPmask = niftiread('..\masks\putamen_left.nii');
RPmask = niftiread('..\masks\putamen_right.nii');
LGPmask = niftiread('..\masks\leftPallidum.nii');
RGPmask = niftiread('..\masks\rightPallidum.nii');
LTHmask = niftiread('..\masks\leftThalamus.nii');
RTHmask = niftiread('..\masks\rightThalamus.nii');

% LCmask = shiftdim(LCmask, 1);
% RCmask = shiftdim(RCmask, 1);
% LPmask = shiftdim(LPmask, 1);
% RPmask = shiftdim(RPmask, 1);
LCmask = transpose_mask(LCmask);
RCmask = transpose_mask(RCmask);
LPmask = transpose_mask(LPmask);
RPmask = transpose_mask(RPmask);
LGPmask = transpose_mask(LGPmask);
RGPmask = transpose_mask(RGPmask);
LTHmask = transpose_mask(LTHmask);
RTHmask = transpose_mask(RTHmask);

LCmask = double(LCmask(:,:,zLimL:zLimU) > 0.5);
RCmask = double(RCmask(:,:,zLimL:zLimU) > 0.5);
LPmask = double(LPmask(:,:,zLimL:zLimU) > 0.5);
RPmask = double(RPmask(:,:,zLimL:zLimU) > 0.5);
LGPmask = double(LGPmask(:,:,zLimL:zLimU) > 0.5);
RGPmask = double(RGPmask(:,:,zLimL:zLimU) > 0.5);
LTHmask = double(LTHmask(:,:,zLimL:zLimU) > 0.5);
RTHmask = double(RTHmask(:,:,zLimL:zLimU) > 0.5);

% LSmask = double(LCmask > 0.5 | LPmask > 0.5);
% RSmask = double(RCmask > 0.5 | RPmask > 0.5);

% striatum masks are created in "normalization\create_striatum_masks"
load('..\masks\striatum_masks_from_mean_control_image.mat');
LSmask = LSmask(:,:,zLimL:zLimU);
RSmask = RSmask(:,:,zLimL:zLimU);
WSmask = WSmask(:,:,zLimL:zLimU);
SKmask = SKmask(:,:,zLimL:zLimU);

WSmask = force_symmetric(WSmask);
SKmask = force_symmetric(SKmask);


if 0
    LoccMask = niftiread('..\masks\occl_l_.nii');
    RoccMask = niftiread('..\masks\occr_r_.nii');
    
    LoccMask = transpose_mask(LoccMask);
    RoccMask = transpose_mask(RoccMask);
    
    LoccMask = double(LoccMask(:,:,zLimL:zLimU) > 0.5);
    RoccMask = double(RoccMask(:,:,zLimL:zLimU) > 0.5);
    
    occMask = double(LoccMask > 0.5 | RoccMask > 0.5);
else
    % Use our occipital lobe mask
    load('..\masks\occMask');
end

% load('..\masks\leftStriatum');
% LSmask=mask;
% load('..\masks\rightStriatum');
% RSmask=mask;
% load('..\masks\leftCaudate');
% LCmask=mask;
% load('..\masks\rightCaudate');
% RCmask=mask;
% load('..\masks\leftPutamen');
% LPmask=mask;
% load('..\masks\rightPutamen');
% RPmask=mask;
% clear mask

end


function WSmask = force_symmetric(WSmask)
% Since the center of the MNI atlas is at the 45th column, and the images 
% are also registered to the MNI atlas, we flip the whole striatum mask 
% along the 45th column and make it symmetric along this column.
WSmask1 = fliplr(WSmask); WSmask1 = WSmask1(:,[3:91,91,91],:);
WSmask = (WSmask > 0.5) | (WSmask1 > 0.5);
WSmask = double(WSmask);
end


function show_mask(mask)
figure;
ind = 1; 
for i = 1:5
    for j = 1:6
        if ind <= size(mask,3)
            subplot(5,6,ind); 
            imshow(mask(:,:,ind),[]); 
            ind = ind+1; 
        end
    end 
end

end