function show_seq_quiver_in_figure(train_data, train_ts, train_ids, ...
    labels, r_sk, model)

if 1
    % 2 by 4 subplots where the first column is the entire data, the second
    % to the last columns are the 3 clusters.
    show_whole_data_and_clusters_together(train_data, train_ts, ...
        train_ids, labels, r_sk, model);
else
    % One figure showing 2 by 1 subplots of the entire data, one figure
    % showing 2 by 3 subplots of the 3 clusters.
    show_whole_data_and_clusters_separately(train_data, train_ts, ...
        train_ids, labels, r_sk, model);
end

end


function show_whole_data_and_clusters_separately(train_data, train_ts,  ...
    train_ids, labels, r_sk, model)

plot_rows = 2;
plot_cols = 1;
options = [];
options.create_figure = 0;
options.subplot_inds = {{plot_rows,plot_cols,1},{plot_rows,plot_cols,2}};
options.all_elapsed_days = train_ts;
options.relative_arrow_size = 1/70;

figure;

%% plot whole data
show_seq_quiver_plot(train_data, train_ids, labels, options);
subplot(plot_rows,plot_cols,1); xlabel({'LC','(a) LC vs RC for all subjects'});
subplot(plot_rows,plot_cols,2); xlabel({'LP','(b) LP vs RP for all subjects'});

create_time_series_colorbar();
set(gcf,'position', [0.448, 0.010, 0.565, 0.846]);

[M,class_inds] = max(r_sk, [], 2);

% move subplot 2 higher to save space
hs = findobj(gcf,'type','axes');
pos2 = get(hs(2),'position');
pos2(2) = pos2(2) * 0.9;
set(hs(2),'position',pos2);

%% plot the 3 clusters
figure;

plot_cols = 3;

subplot_inds{1} = {{plot_rows,plot_cols,1},{plot_rows,plot_cols,4}};
subplot_inds{2} = {{plot_rows,plot_cols,2},{plot_rows,plot_cols,5}};
subplot_inds{3} = {{plot_rows,plot_cols,3},{plot_rows,plot_cols,6}};

%     data_interest{1} = [1.9501; 1.3787];
%     data_interest{2} = [1.601; 1.0843];
%     data_interest{3} = [0.94443; 1.4387];

data_interest{1} = [1; 0.75; 1; 1.5];
data_interest{2} = [1; 0.75; 1; 1.5];
data_interest{3} = [1; 0.75; 1; 1.5];

subplot_start_char = 97; % start at '(a)'

options.relative_arrow_size = 1/70;

for k = [1,2,3]
    train_data_1 = train_data(class_inds == k);
    train_ts_1 = train_ts(class_inds == k);
    train_ids_1 = train_ids(class_inds == k);


    options.subplot_inds = subplot_inds{k};
    options.all_elapsed_days = train_ts_1;
    show_seq_quiver_plot(train_data_1, train_ids_1, labels, options);
    for i = 1:length(subplot_inds{k})
        subplot(subplot_inds{k}{i}{:});
        set_axes_range(train_data);
    end

    if 0 % find the subject with closest time point to the specified one
        sub_idx = find_specific_subject(train_data_1, [1,4], data_interest{k});
        data_sel = train_data_1(sub_idx);
    else % directly use the specified time point
        data_sel = {data_interest{k}};
    end

    for i = 1:length(data_sel)
        data_sel{size(data_interest{k},2) + i} = flipud(data_sel{i});
    end
    for i = 1:length(data_sel)
        ts_sel{i} = 20*365;
    end
    if ~isempty(model)
        show_solution_curve_and_xlabel(model.As(:,:,k), options, k, ...
            data_sel, ts_sel, subplot_start_char);
    end
    
    if k == 1
        create_time_series_colorbar();
    end
end

create_model_trajectory_colorbar();

set(gcf,'position', [0.17, 0.13, 0.71, 0.80]);

end


function show_whole_data_and_clusters_together(train_data, train_ts,  ...
    train_ids, labels, r_sk, model)

plot_rows = 2;
plot_cols = 4;
options = [];
options.create_figure = 0;
options.subplot_inds = {{plot_rows,plot_cols,1},{plot_rows,plot_cols,5}};
options.all_elapsed_days = train_ts;

figure;

%% plot whole data
show_seq_quiver_plot(train_data, train_ids, labels, options);
subplot(plot_rows,plot_cols,1); xlabel({'LC','LC vs RC for all subjects'});
subplot(plot_rows,plot_cols,5); xlabel({'LP','LP vs RP for all subjects'});

create_time_series_colorbar();

[M,class_inds] = max(r_sk, [], 2);

%% plot the 3 clusters
subplot_inds{1} = {{plot_rows,plot_cols,2},{plot_rows,plot_cols,6}};
subplot_inds{2} = {{plot_rows,plot_cols,3},{plot_rows,plot_cols,7}};
subplot_inds{3} = {{plot_rows,plot_cols,4},{plot_rows,plot_cols,8}};

%     data_interest{1} = [1.9501; 1.3787];
%     data_interest{2} = [1.601; 1.0843];
%     data_interest{3} = [0.94443; 1.4387];

data_interest{1} = [1; 0.75; 1; 1.5];
data_interest{2} = [1; 0.75; 1; 1.5];
data_interest{3} = [1; 0.75; 1; 1.5];

subplot_start_char = 99; % start at '(c)'
subplot_start_char = [];

for k = [1,2,3]
    train_data_1 = train_data(class_inds == k);
    train_ts_1 = train_ts(class_inds == k);
    train_ids_1 = train_ids(class_inds == k);


    options.subplot_inds = subplot_inds{k};
    options.all_elapsed_days = train_ts_1;
    show_seq_quiver_plot(train_data_1, train_ids_1, labels, options);
    for i = 1:length(subplot_inds{k})
        subplot(subplot_inds{k}{i}{:});
        set_axes_range(train_data);
    end

    if 0 % find the subject with closest time point to the specified one
        sub_idx = find_specific_subject(train_data_1, [1,4], data_interest{k});
        data_sel = train_data_1(sub_idx);
    else % directly use the specified time point
        data_sel = {data_interest{k}};
    end

    for i = 1:length(data_sel)
        data_sel{size(data_interest{k},2) + i} = flipud(data_sel{i});
    end
    for i = 1:length(data_sel)
        ts_sel{i} = 20*365;
    end
    show_solution_curve_and_xlabel(model.As(:,:,k), options, k, data_sel, ts_sel, ...
        subplot_start_char);
end

create_model_trajectory_colorbar();

set(gcf,'position', [0.0805    0.1986    0.7974    0.6083]);

end

%% Supporting functions

function create_time_series_colorbar()
ticks = [0,0.2,0.4,0.6,0.8,1];
ticklabels = {'Baseline','1 yr','2 yrs', '3 yrs', '4 yrs', '5 yrs'};
cb = colorbar(gca,'east','Direction','reverse','ticks',ticks,'TickLabels', ...
    ticklabels,'AxisLocation','in');
colormap(cb, 'parula');
end

function create_model_trajectory_colorbar()
ticks = [0,0.2,0.4,0.6,0.8,1];
ticklabels = {'Baseline','4 yrs','8 yrs', '12 yrs', '16 yrs', '20 yrs'};
cb = colorbar(gca,'east','Direction','reverse','ticks',ticks,'TickLabels',ticklabels);
colormap(cb, 'jet');
end

function sub_idx = find_specific_subject(train_data, inds, data_interest)
sub_idx = [];
for k = 1:size(data_interest,2)
    dist = [];
    for i = 1:length(train_data)
        dist(i) = min(mean(abs(train_data{i}(inds,:) - repmat(data_interest(:,k), ...
            [1 size(train_data{i},2)]))));
    end
    [~,min_i] = min(dist);
    sub_idx(end+1) = min_i;
end
end

function show_solution_curve_and_xlabel(A, options, cluster_idx, data, ts, subplot_start_char)
% calculate solution curve
[xts,ts] = calc_LDS_curve(data, ts, A);

cmap = colormap('jet');
% plot LC vs RC
subplot(options.subplot_inds{1}{:});
for i = 1:length(xts)
    plot(xts{i}(1,:),xts{i}(4,:),'-k','linewidth',3);
    plot_line_gradient_color(xts{i}(1,:),xts{i}(4,:),'-',1.5,cmap);
end

if isempty(subplot_start_char)
    subplot_char = [];
else
    subplot_char = char(subplot_start_char-1+cluster_idx*2-1);
    subplot_char = ['(',subplot_char,') '];
end

xlabel({'LC',[subplot_char, 'LC vs RC for subtype ',num2str(cluster_idx)]});

% plot LP vs RP
subplot(options.subplot_inds{2}{:});
for i = 1:length(xts)
    plot(xts{i}(2,:),xts{i}(3,:),'-k','linewidth',3);
    plot_line_gradient_color(xts{i}(2,:),xts{i}(3,:),'-',1.5,cmap);
end

if isempty(subplot_start_char)
    subplot_char = [];
else
    subplot_char = char(subplot_start_char-1+cluster_idx*2);
    subplot_char = ['(',subplot_char,') '];
end

xlabel({'LP',[subplot_char, 'LP vs RP for subtype ',num2str(cluster_idx)]});

end

function plot_line_gradient_color(x,y,style,linewidth,colormap)
npts = length(x);
interval = round(npts / size(colormap,1));
for i = 1:interval:npts-1   
    plot(x(i:min(i+interval, npts)),y(i:min(i+interval, npts)), ...
        style,'linewidth',linewidth, ...
        'color',colormap(ceil(i/npts*size(colormap,1)),:));
end
end



function set_axes_range(all_image_data)
X = cell2mat(all_image_data);
data_max = max(X(:));
xlim([0 data_max]);
ylim([0 data_max]);

end

