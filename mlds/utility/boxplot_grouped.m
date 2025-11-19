function hs = boxplot_grouped(data, colors)
% plot boxplot for grouped data.
% data - G X C cell array. C is the number of category (e.g. train/test).
%        G is the number of groups (e.g. folder 1,2,...). Each cell of the
%        array is a column vector containing data.
% colors - 1 X C cell array containing colors.

hs = {};
for i = 1:size(data,2)
    X = [];
    g = [];
    for j = 1:size(data,1)
        X = [X; data{j,i}];
        g = [g; repmat(j, length(data{j,i}), 1)];
    end
    position = [1:size(data,1)] + (i-1)*0.3;
    box_i = boxplot(X,g,'colors',colors{i},'positions',position,'width',0.18,...
        'plotstyle','compact','labelorientation','horizontal', ...
        'outliersize', 2); 
    hold on
    hs = cat(2, hs, {box_i});
end

% xlabels = strsplit(num2str([1:size(data,1)]));
% set(gca,'XTickLabel',xlabels)  % Erase xlabels   

end
