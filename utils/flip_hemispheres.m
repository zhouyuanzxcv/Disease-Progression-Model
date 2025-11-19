function [Y, flip_flag] = flip_hemispheres(train_data)
%FLIP_HEMISPHERES Flip hemispheres s.t. the left hemisphere is always
%brighter than the right one.
dim = size(train_data{1},1);

flip_flag = zeros(size(train_data));
% flip left and right if left is darker than right
for i = 1:length(train_data)
    Xi1 = mean(train_data{i},2);
    if mean(Xi1(1:dim/2)) < mean(Xi1(dim/2+1:end))
        train_data{i} = flipud(train_data{i});
        flip_flag(i) = 1;
    end
end

Y = train_data;

end

