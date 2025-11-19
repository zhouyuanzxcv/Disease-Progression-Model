function mask1 = transpose_mask(mask)
% transpose mni mask
for i = 1:size(mask,3)
    mask1(:,:,i) = flipud(mask(:,:,i)');
end
end


