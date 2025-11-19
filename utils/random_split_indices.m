function [indices1,indices2] = random_split_indices(N, ratio)
%RANDOM_SPLIT_INDICES Summary of this function goes here
%   Detailed explanation goes here
p = randperm(N);
split_ind = round(N*ratio);
indices1 = p(1:split_ind);
indices2 = p(split_ind+1:end);
end

