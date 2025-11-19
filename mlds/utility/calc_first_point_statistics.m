function [mu, Sigma] = calc_first_point_statistics(X)
D = size(X{1},1);
mu = zeros(1,D);
X1 = [];
for i = 1:length(X)
    x1 = X{i}(:,1);
    X1 = cat(1, X1, x1');
end
Sigma = (X1'*X1) / length(X);

end

