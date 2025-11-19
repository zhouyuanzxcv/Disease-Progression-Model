function masks = load_masks_structure()
%LOAD_MASKS_STRUCTURE Summary of this function goes here
%   Detailed explanation goes here
[occMask,LSmask,RSmask,LCmask,RCmask,LPmask,RPmask,WSmask, ...
    LGPmask, RGPmask, LTHmask, RTHmask, SKmask] = load_masks();
masks = [];
masks.occ = occMask;
masks.LS = LSmask;
masks.RS = RSmask;
masks.LC = LCmask;
masks.RC = RCmask;
masks.LP = LPmask;
masks.RP = RPmask;
masks.WS = WSmask;
masks.LGP = LGPmask;
masks.RGP = RGPmask;
masks.LTH = LTHmask;
masks.RTH = RTHmask;
masks.SK = SKmask; % skull
end

