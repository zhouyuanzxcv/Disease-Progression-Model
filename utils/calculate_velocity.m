function [V, delta_tau, X_base] = calculate_velocity(X, ts)
% X - cell array of data, each has a size of number of voxels x number of time points.
% ts - cell array of time stamps (unit: day), will be converted to year.
V = {};
delta_tau = {};
X_base = {};
for i = 1:length(X)
    if size(X{i},2) == 1
        delta_tau{i} = [];
        V{i} = [];
        X_base{i} = [];
    else
        delta_tau{i} = (ts{i}(2:end) - ts{i}(1:end-1)) / 365;
        dist = X{i}(:,2:end) - X{i}(:,1:end-1);
        delta_t = repmat(delta_tau{i}, [size(dist,1) 1]);
        V{i} = dist ./ delta_t;
        X_base{i} = X{i}(:,1:end-1);
    end
end
end
