function [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
    LGPmask, RGPmask, LTHmask, RTHmask, SKmask] = load_masks()
%RETRIEVE_META_DATA This function will be used throughout the data creation
% step so changing the implementation here will affect all the steps.

if 0 % load LR flipped and registered masks
    load('../masks/flipped_registered_masks.mat');
elseif 1 % load dilated masks from the original masks
    [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask, ...
        WSmask, LGPmask, RGPmask, LTHmask, RTHmask, SKmask] = ...
        load_dilated_masks(1);
    
    WSmask = dilate_WSmask(WSmask, 2);
    SKmask = dilate_mask(SKmask, -4);

elseif 0 % load original masks
    [occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
        LGPmask, RGPmask, LTHmask, RTHmask] = load_masks_original();
else % load Hemant's masks
    load('..\masks\occMask');
    
    load('..\masks\leftStriatum');
    LSmask=mask;
    load('..\masks\rightStriatum');
    RSmask=mask;
    load('..\masks\leftCaudate');
    LCmask=mask;
    load('..\masks\rightCaudate');
    RCmask=mask;
    load('..\masks\leftPutamen');
    LPmask=mask;
    load('..\masks\rightPutamen');
    RPmask=mask;
    clear mask
end

end


function WSmask = dilate_WSmask(WSmask, dilate_size)
WSmask = dilate_mask(WSmask, dilate_size);

% Since the center of the MNI atlas is at the 45th column, and the images 
% are also registered to the MNI atlas, we flip the whole striatum mask 
% along the 45th column and make it symmetric along this column.
WSmask1 = fliplr(WSmask); WSmask1 = WSmask1(:,[3:91,91,91],:);
WSmask = (WSmask > 0.5) | (WSmask1 > 0.5);
WSmask = double(WSmask);
end

% 
% function show_mask(mask)
% figure;
% ind = 1; 
% for i = 1:5
%     for j = 1:6
%         if ind <= size(mask,3)
%             subplot(5,6,ind); 
%             imshow(mask(:,:,ind),[]); 
%             ind = ind+1; 
%         end
%     end 
% end
% 
% end