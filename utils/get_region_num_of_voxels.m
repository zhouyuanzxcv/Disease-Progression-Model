function [nums,labels] = get_region_num_of_voxels(mask)
%GET_REGION_NUM_OF_VOXELS Summary of this function goes here
%   Detailed explanation goes here
occMask = mask.Occ;
LSmask = mask.LS;
RSmask = mask.RS;
LCmask = mask.LC;
RCmask = mask.RC;
LPmask = mask.LP;
RPmask = mask.RP;

nums = [count_num(LCmask), count_num(RCmask)];
nums = [nums, count_num(LPmask), count_num(RPmask)];
nums = [nums, count_num(LSmask), count_num(RSmask)];
nums = [nums, count_num(occMask), count_num(occMask)];

labels = {'LC','RC','LP','RP','LS','RS','Scale','Scale2'};

end

function num = count_num(mask)
num = length(find(mask > 0.5));
end