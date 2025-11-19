function show_basis_on_image(v, sizes, mask_inds, options)
%SHOW_BASIS_ON_IMAGE Summary of this function goes here
%   Detailed explanation goes here
%
% The difference between this function and "show_image_w_mask.m" is that
% the basis (v) can take arbitrary values while masks take range [0,1].

if nargin < 4
    options = [];
end

X = parse_param(options, 'X', []);
v_range = parse_param(options, 'v_range', [min(v), max(v)]);
fh_image_w_mask = parse_param(options, 'fh_image_w_mask', []);
masks = parse_param(options, 'masks', []);
selected_slices = parse_param(options, 'selected_slices', 1:sizes(3));
montage_size = parse_param(options, 'montage_size', [NaN NaN]);
roi = parse_param(options, 'roi', []);
colorbar_tick_num = parse_param(options, 'colorbar_tick_num', 5);

if isempty(fh_image_w_mask)
    fh_image_w_mask = figure;
end

if isempty(roi)
    roi.xs = [1:sizes(2)];
    roi.ys = [1:sizes(1)];
end

% mask_brightness = 0.2; % default value

img = zeros(sizes);
v_img = zeros(sizes);

if ~isempty(X)
    img(mask_inds) = X;
    img = img / max(max(img(:)), 1e-9); % the maximum is 1 
    % mask_brightness = 1 - max(img(:));
else
    img = get_background_image();
end

v_img(mask_inds) = v;
v_img_abs = abs(v_img);
v_img_n = v_img_abs / max(abs(v_range));


[rows,cols,slices] = size(img);

img_color = reshape(img,[rows,cols,1,slices]);
img_color = repmat(img_color, [1,1,3,1]);


% impose v on the image
mask_r = zeros(rows,cols,slices);
mask_g = zeros(rows,cols,slices);
mask_b = zeros(rows,cols,slices);
mask_r(v_img > 0) = 1;
mask_b(v_img < 0) = 1;

% superimpose the coefficients of v on the image
img_per_oper = img_color(:,:,:,:);

img_per_oper(:,:,1,:) = squeeze(img_per_oper(:,:,1,:)) .* (1-v_img_n) + mask_r .* v_img_n;
img_per_oper(:,:,2,:) = squeeze(img_per_oper(:,:,2,:)) .* (1-v_img_n) + mask_g .* v_img_n;
img_per_oper(:,:,3,:) = squeeze(img_per_oper(:,:,3,:)) .* (1-v_img_n) + mask_b .* v_img_n;

if ~isempty(masks)
% if 0
    % superimpose the masks of caudate, putamen, and occipital lobe on the
    % image
    colors = get_mask_colors();
    masks1 = cat(4, masks.LC + masks.RC, masks.LP + masks.RP, masks.WS);
    colors1 = 0.6*[colors.caudate; colors.putamen; colors.striatum];
    
    for k = 1:size(masks1,4)
        mask1 = masks1(:,:,:,k);
        [fx,fy,fz] = gradient(mask1);
        grad_mag = sqrt(fx.^2 + fy.^2 + fz.^2);
        mask1 = grad_mag / max(grad_mag(:));
        img_per_oper(:,:,1,:) = squeeze(img_per_oper(:,:,1,:)) + ...
            mask1 * colors1(k,1);
        img_per_oper(:,:,2,:) = squeeze(img_per_oper(:,:,2,:)) + ...
            mask1 * colors1(k,2);
        img_per_oper(:,:,3,:) = squeeze(img_per_oper(:,:,3,:)) + ...
            mask1 * colors1(k,3);
    end
end

img_color(:,:,:,:) = img_per_oper;

% show image with mask
ax_h = findobj(fh_image_w_mask,'type','axes');
img_color1 = img_color(:,:,:,:);

if isempty(ax_h)
    figure(fh_image_w_mask);
    ax_h = gca;
end
    
montage(img_color1(roi.ys,roi.xs,:,:), 'Indices', selected_slices, 'Size', montage_size, 'parent', ax_h);

%% show colorbar
cmap = zeros(64, 3);
v_range_n = v_range / max(abs(v_range));
v_interp = linspace(v_range_n(1), v_range_n(2), 64)';
v_ind_p = find(v_interp>0);
v_ind_m = find(v_interp<0);

% bg_value = mean(img(mask_inds));
alpha_p = v_interp(v_ind_p);
cmap(v_ind_p,:) = alpha_p * [1 0 0] + (1-alpha_p) * [1 1 1];
alpha_m = -v_interp(v_ind_m);
cmap(v_ind_m,:) = alpha_m * [0 0 1] + (1-alpha_m) * [1 1 1];

% if v_range = [-0.11, 0.09], prec = 2, i.e. the number of digits required after 0
% prec = ceil(log10(1 / max(abs(v_range)))) + 1;
% if prec < 0, prec = 0; end
% prec_format = ['%.',num2str(prec),'f'];

num_labels = colorbar_tick_num; % approximately 5 labels
tick_interval = (v_range(2) - v_range(1)) / (num_labels - 1);
prec1 = -floor(log10(tick_interval));
prec_format = ['%.',num2str(prec1),'f'];
tick_interval = floor(tick_interval * (10^prec1)) * (10^(-prec1));

ticklabels = {};
ticks_actual = [];
tick_display = sign(v_range(1)) * floor(abs(v_range(1)) / tick_interval) * tick_interval;
while tick_display < v_range(2)
    ticklabels{end+1} = num2str(tick_display, prec_format);
    tick_actual = (tick_display - v_range(1)) / (v_range(2) - v_range(1));
    ticks_actual(end+1) = tick_actual;
    tick_display = tick_display + tick_interval;
end

cb = colorbar('eastoutside','ticks',ticks_actual,'TickLabels',ticklabels); % 'TickLabels'
colormap(cb, cmap);


end

function img = get_background_image()
load('../high_dimension/data/mean_control.mat');
img = img_mean;

load('../outliers/metaData');
zLimU=metaData.zLim.uLim;
zLimL=metaData.zLim.lLim;
clear metaData

img = img(:,:,zLimL:zLimU);
img = img / max(img(:)) * 0.8;

end
