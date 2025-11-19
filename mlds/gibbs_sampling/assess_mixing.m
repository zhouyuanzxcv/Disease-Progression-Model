function [num_attempts, ratios, means] = assess_mixing(X)
[num_iter,num_chains] = size(X);
num_attempts = 200:100:num_iter;
ratios = [];
means = [];

for i = num_attempts
    % discard the first half and keep the second half
    X1 = X(round(i/2)+1:i,:);
    
    % break the valid samples into two chains
    num = size(X1,1);
    X2 = [X1(1:num/2,:),X1(num/2+1:end,:)];
    
    X2(:,any(isnan(X2), 1)) = [];
    [n,m] = size(X2);
    
    % calculate between sequence variance
    mean_per_chain = mean(X2, 1);
    mean_all = mean(mean_per_chain);
    B = sum((mean_per_chain - mean_all).^2)*(n/(m-1));
    
    % calculate within sequence variance
    s2 = sum((X2 - repmat(mean_per_chain, [n 1])).^2, 1) / (n-1);
    W = sum(s2) / m;
    
    % calculate marginal posterior variance
    vars = (n-1)/n*W + B/n;
    
    ratios(end+1) = sqrt(vars / W);
    means(end+1) = mean_all;
end


end

