function mask1 = dilate_mask(mask, dilate_size)
if dilate_size > 0
    se = strel('sphere', dilate_size);
    mask1 = imdilate(mask, se);
else
    se = strel('sphere', -dilate_size);
    mask1 = imerode(mask, se);
end

sigma = 0.5;
mask1 = imgaussfilt3(mask1, sigma, 'filtersize', 3, ...
    'filterdomain', 'spatial');
    
end


