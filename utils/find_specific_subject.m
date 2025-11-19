function subjects_sel = find_specific_subject(train_data, inds, data_interest, train_ids)
X = cell2mat(train_data);
mag = median(abs(X), 2);

diffs = [];
subjects_sel = [];
for i = 1:length(train_data)
    Ti = size(train_data{i},2);
    diff = abs(train_data{i}(inds,:) - repmat(data_interest, [1 Ti]));
    thresh = repmat(1e-3 * mag(inds), [1 Ti]);
    [M,ind] = min(mean(diff, 1));
    if any(all(diff - thresh < 0, 1))
        disp(['Find ',train_ids{i},' at time ',num2str(ind)]);
        subjects_sel(end+1) = str2num(train_ids{i});
    end
    diffs(i) = M;
end

[M,ind] = min(diffs);
disp(['Min diff: ', num2str(M),', Subject id: ',train_ids{ind}]);

end