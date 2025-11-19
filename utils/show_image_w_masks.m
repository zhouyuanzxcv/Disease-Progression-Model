function img_color = show_image_w_masks(img_ori, list_oper_slices, list_masks,...
    list_colors, fh_image_w_mask, options)
% Example:
%   list_oper_slices{1} = [zLimL:zLimU];
%   list_masks{1} = cat(4,LCmask+RCmask,LPmask+RPmask,occMask);
%   list_colors{1} = [0.2,0,0;0,0.2,0;0,0,0.2];
%   list_styles{1} = {'contour','contour,'area'};
%   list_oper_slices{2} = [1:91];
%   list_masks{2} = 0;
%   list_colors{2} = [0.1,0.1,0];
%   list_styles{2} = {'area'}; 
%   options = [];
%   options.list_styles = list_styles;
%   options.selected_slices = 1:size(img_ori,3);

if nargin < 6
    options = [];
end

list_styles = parse_param(options,'list_styles',[]);
selected_slices = parse_param(options, 'selected_slices', 1:size(img_ori,3));
montage_size = parse_param(options, 'montage_size', [NaN NaN]);

[rows,cols,slices] = size(img_ori);

img_color = reshape(img_ori,[rows,cols,1,slices]);
img_color = repmat(img_color, [1,1,3,1]);

for i = 1:length(list_oper_slices)
    oper_slices = list_oper_slices{i};
    masks = list_masks{i};
    colors = list_colors{i};
    
    img_per_oper = img_color(:,:,:,oper_slices);

    for j = 1:size(masks,4)
        style = 'area';
        if ~isempty(list_styles) && ~isempty(list_styles{i})
            style = list_styles{i}{j};
        end
        
        mask1 = masks(:,:,:,j);
        if strcmp(style, 'area')
            % do nothing
        elseif strcmp(style, 'contour')
            [fx,fy,fz] = gradient(mask1);
            grad_mag = sqrt(fx.^2 + fy.^2 + fz.^2);
            mask1 = 2 * grad_mag / max(grad_mag(:));
        end
        img_per_oper(:,:,1,:) = squeeze(img_per_oper(:,:,1,:)) + ...
            mask1*colors(j,1);
        img_per_oper(:,:,2,:) = squeeze(img_per_oper(:,:,2,:)) + ...
            mask1*colors(j,2);
        img_per_oper(:,:,3,:) = squeeze(img_per_oper(:,:,3,:)) + ...
            mask1*colors(j,3);
    end

    img_color(:,:,:,oper_slices) = img_per_oper;
end

ax_h = findobj(fh_image_w_mask,'type','axes');
if isempty(ax_h)
    figure(fh_image_w_mask);
    montage(img_color, 'Indices', selected_slices, 'Size', montage_size);
else
    montage(img_color, 'Indices', selected_slices, 'size', montage_size, 'parent',ax_h);
end

end

