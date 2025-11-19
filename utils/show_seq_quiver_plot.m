function show_seq_quiver_plot(all_image_data, subject_ids, labels, options)
if nargin < 4
    options = [];
end

values_for_probe = {};

if ~isempty(subject_ids)
    for i = 1:4
        values = [];
        subject_ids_expanded = {};
        for n = 1:length(all_image_data)
            values = [values, all_image_data{n}(i,:)];
            subject_ids_expanded = [subject_ids_expanded,  ...
                repmat(subject_ids(n), [1, size(all_image_data{n},2)])];
        end
        
        values_for_probe{i} = values;
    end
end

if ~isempty(labels) 
    if strcmp(labels{1}, 'LC') && strcmp(labels{2}, 'LP')
        % Since all_image_data are ordered in (LC, LP, RP, RC), the (LC, RC, LP,
        % RP) order is (1, 4, 2, 3)
        labels_expected = parse_param(options, 'labels_expected', {'LC','RC','LP','RP'});
        [labels_expected, order_inds] = reorder_labels(labels, labels_expected);
        % order_inds = [1 4 2 3]; 
    else
        labels_expected = labels;
        order_inds = (1:size(all_image_data{1},1));
    end
else
    for j = 1:size(all_image_data{1},1)
        labels{j} = ['Dim ',num2str(j)];
    end
    order_inds = (1:size(all_image_data{1},1));
    labels_expected = labels;
end

classes = parse_param(options,'classes',[]);
predicted_data = parse_param(options,'predicted_data',[]);
options.order_inds = order_inds;
options.labels_expected = labels_expected;

if ~isempty(classes)
    show_time_sq_as_quiver_with_classes(all_image_data, options);
elseif ~isempty(predicted_data)
    show_time_sq_as_quiver_with_prediction(all_image_data, options);
else
    fh_quiver = show_time_sq_as_quiver(all_image_data, order_inds, options);
    dcm_obj = datacursormode(fh_quiver);
    dcm_quiver_c = @(obj, event) show_extra_for_quiver(obj, event, ...
        values_for_probe{1}, values_for_probe{2}, subject_ids_expanded);
    dcm_obj.UpdateFcn = dcm_quiver_c;
end
end

function fh_quiver = show_time_sq_as_quiver(all_image_data, order_inds, options)

create_figure = parse_param(options, 'create_figure', 1);
subplot_inds = parse_param(options, 'subplot_inds', {{1,2,1},{1,2,2}});

if create_figure
    fh_quiver = figure;
else
    fh_quiver = gcf;
end

labels_expected = options.labels_expected;
options.original_all_image_data = all_image_data;

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.1, 0.9, 0.8])
subplot(subplot_inds{1}{:}); 
show_time_sq_as_quiver_impl(all_image_data, order_inds(1), order_inds(2), options); % LC vs. RC
xlabel(labels_expected{1}); ylabel(labels_expected{2});

subplot(subplot_inds{2}{:}); 
show_time_sq_as_quiver_impl(all_image_data, order_inds(3), order_inds(4), options); % LP vs. RP
xlabel(labels_expected{3}); ylabel(labels_expected{4});

end


function set_axes_range(all_image_data, options)
axis_mode = parse_param(options,'axis_mode',1);
if axis_mode == 1
    axis equal;
end

[x_range, y_range] = get_axes_range(all_image_data, options);
xlim([x_range(1) x_range(2)]);
ylim([y_range(1) y_range(2)]);

end

function [x_range,y_range] = get_axes_range(all_image_data, options)
axis_mode = parse_param(options,'axis_mode',1);
if axis_mode == 1
    X = cell2mat(all_image_data);
    data_max = max(X(:));
    data_min = 0;
    range = [data_min data_max];
    x_range = range;
    y_range = range;
elseif axis_mode == 2
    x_data = [];
    y_data = [];
    for i = 1:length(all_image_data)
        x_data = [x_data, all_image_data{i}(options.x_ind,:)];
        y_data = [y_data, all_image_data{i}(options.y_ind,:)];
    end
    x_range = [min(x_data), max(x_data)];
    y_range = [min(y_data), max(y_data)];
end
end

function draw_ref_line()
x = xlim;
y = ylim;
xmin = min(x(2),y(2));
plot([0 xmin],[0 xmin],'k--','linewidth',0.5);
end

function qhs = show_time_sq_as_quiver_impl(all_image_data, x_ind, y_ind, options)
qhs = [];

if 1
if isfield(options, 'all_elapsed_days')
    all_elapsed_days = options.all_elapsed_days;
    if isfield(options, 'max_elapsed_days')
        max_elapsed_days = options.max_elapsed_days;
    else
        max_elapsed_days = 5*365;
    end

    qhs = show_time_sq_as_quiver_impl_cmap(all_image_data, ...
        all_elapsed_days, x_ind, y_ind, max_elapsed_days, options);
    return;
end
end


num_subjects = length(all_image_data);
for i = 1:length(all_image_data)
    Ts(i) = size(all_image_data{i},2);
end
T_max = max(Ts);

for t = 1:T_max
    lc = [];
    rc = [];
    for i = 1:num_subjects
        if size(all_image_data{i},2) >= t+1
            lc = cat(1, lc, all_image_data{i}(x_ind,t:t+1));
            rc = cat(1, rc, all_image_data{i}(y_ind,t:t+1));
        end
    end
    if ~isempty(lc) && ~isempty(rc)
        hold on;
        qhs(t) = quiver(lc(:,1), rc(:,1), lc(:,2)-lc(:,1), rc(:,2)-rc(:,1), 0, ...
            'AutoScale','off','linewidth',1);
        set_axes_range(all_image_data, options);
    end
end
draw_ref_line();

end

function qhs = show_time_sq_as_quiver_impl_cmap(all_image_data, ...
    all_elapsed_days, x_ind, y_ind, max_t, options)
qhs = [];

cmap = colormap('parula');
num_colors = size(cmap,1);

num_subjects = length(all_image_data);

hs = {};

options.x_ind = x_ind;
options.y_ind = y_ind;
[x_range, y_range] = get_axes_range(options.original_all_image_data, options);

range = sqrt((x_range(2) - x_range(1)) * (y_range(2) - y_range(1)));

relative_arrow_size = parse_param(options, 'relative_arrow_size', 1/50);
arrow_size = range * relative_arrow_size;

arrow_angle = pi/10;

hold on;
for i = 1:num_subjects
    image_data = all_image_data{i};
    for j = 1:size(image_data,2)-1
        x = image_data(x_ind,j:j+1);
        y = image_data(y_ind,j:j+1);
        angle = atan2(y(2)-y(1), x(2)-x(1)) + pi;
        x1 = arrow_size * cos(angle + arrow_angle) + x(2);
        y1 = arrow_size * sin(angle + arrow_angle) + y(2);
        x2 = arrow_size * cos(angle - arrow_angle) + x(2);
        y2 = arrow_size * sin(angle - arrow_angle) + y(2);
        x_w_arrow = [x,x1,x2,x(2)];
        y_w_arrow = [y,y1,y2,y(2)];
        h = plot(x_w_arrow, y_w_arrow, 'linewidth', 1);

        hs{end+1} = h;
    end
end

drawnow

idx = 1;
for i = 1:num_subjects
    ts = all_elapsed_days{i};
    for j = 1:size(ts,2)-1
        start_ind = max(ts(j)/max_t*num_colors, 1);
        end_ind = min(ts(j+1)/max_t*num_colors, num_colors);
        cd = uint8(255*cat(1, cmap(round([start_ind,end_ind]),:)', ones(1,2)));
        cd = [cd, repmat(cd(:,2),1,3)];
        
        set(hs{idx}.Edge, 'ColorBinding', 'interpolated', 'ColorData', cd);
        idx = idx + 1;
    end
end


set_axes_range(options.original_all_image_data, options);

show_ref_line = parse_param(options, 'show_ref_line', 1);
if show_ref_line
    draw_ref_line();
end

end


function fh_quiver = show_time_sq_as_quiver_with_classes(all_image_data, options)
fh_quiver = figure;
set(fh_quiver,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);

classes = parse_param(options,'classes',[]);
weights = parse_param(options,'weights',[]);
all_elapsed_days = parse_param(options,'all_elapsed_days',[]);

class_ids = sort(unique(classes));
order_inds = options.order_inds;

K = length(class_ids);
colors = distinguishable_colors(K);
legend_names = {};
for k = 1:K
    legend_names{k} = ['Subtype ',num2str(class_ids(k))];
end

M = size(all_image_data{1}, 1);
x_inds = order_inds(1:2:M);
y_inds = order_inds(2:2:M);
% x_inds = [order_inds(1), order_inds(3)];
% y_inds = [order_inds(2), order_inds(4)];

labels_expected = options.labels_expected;
x_strs = labels_expected((1:2:M));
y_strs = labels_expected((2:2:M));

% if ~isempty(weights)
%     % normalize weights such that the maximum is one
%     w_max = max(cell2mat(weights));
%     lambda_fcn = @(x) x / w_max;
%     weights = cellfun(lambda_fcn, weights, 'UniformOutput', 0);
% end
options.original_all_image_data = all_image_data;

for i = 1:M/2
    for k = 1:K
        subplot(M/2, K, k+(i-1)*K);

        all_image_data1 = all_image_data(classes == class_ids(k));
        if ~isempty(all_elapsed_days)
            all_elapsed_days1 = all_elapsed_days(classes == class_ids(k));
            options.all_elapsed_days = all_elapsed_days1;
            
            qhs = show_time_sq_as_quiver_impl(all_image_data1, x_inds(i), y_inds(i), options);
        elseif ~isempty(weights)
            weights1 = weights(classes == class_ids(k));
            cmap = colormap('jet');
            qhs = show_time_sq_as_quiver_impl_with_weights(all_image_data1, ...
                x_inds(i), y_inds(i), weights1, cmap, options);
        else
            qhs = show_time_sq_as_quiver_impl(all_image_data1, x_inds(i), y_inds(i), options);
        end
%         for i = 1:length(qhs)
%             set(qhs(i), 'color', colors(k,:));
%         end
        xlabel_str = {x_strs{i}, legend_names{k}};
        ylabel_str = y_strs{i};
        xlabel(xlabel_str);
        ylabel(ylabel_str);
        
        options.x_ind = x_inds(i);
        options.y_ind = y_inds(i);
        set_axes_range(all_image_data, options);
    end
end

if ~isempty(all_elapsed_days)
    show_colorbar(gca);
elseif ~isempty(weights)
    colorbar('Location','eastoutside');
end

end 

function show_colorbar(axe_h)
ticks = [0,0.2,0.4,0.6,0.8,1];
ticklabels = {'Baseline','1 yr','2 yrs', '3 yrs', '4 yrs', '5 yrs'};
cb = colorbar(axe_h,'east','Direction','reverse','ticks',ticks,'TickLabels', ...
    ticklabels,'AxisLocation','in');
colormap(cb, 'parula');

end

function output_txt = show_extra_for_quiver(obj, event, values1, ...
    values2, subject_ids_expanded)
txt_pos = [num2str(event.Position(1)), ', ', num2str(event.Position(2))];

xdata = event.Position(1);
ydata = event.Position(2);

n = find(xdata == values1 & ydata == values2, 1);

if ~isempty(n)
    id = subject_ids_expanded{n};
    idxs = find(strcmp(id, subject_ids_expanded));
    k = find(idxs == n);
    extra_info = ['Sub. ID: ', id, '; Seq. ID: ', num2str(k)];
else
    extra_info = '';
end

output_txt = {txt_pos, extra_info};

end

%% show quivers with prediction
function fh_quiver = show_time_sq_as_quiver_with_prediction(all_image_data, options)
predicted_data = options.predicted_data;
order_inds = options.order_inds;
hide_original_in_prediction = parse_param(options, 'hide_original_in_prediction', 0);
create_figure = parse_param(options, 'create_figure', 1);
subplot_inds = parse_param(options, 'subplot_inds', {{1,2,1},{1,2,2}});
if create_figure
    fh_quiver = figure;    
end

labels_expected = options.labels_expected;

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.1, 0.9, 0.8])
subplot(subplot_inds{1}{:}); 
if ~hide_original_in_prediction
    qhs = show_time_sq_as_quiver_impl(all_image_data, order_inds(1), ...
        order_inds(2), options); % LC vs. RC
    set(qhs, 'color', 'blue');
end
% stack the predicted data
qhs = show_time_sq_as_quiver_prediction_impl(all_image_data, ...
    predicted_data, order_inds(1), order_inds(2), options); % LC vs. RC
set(qhs, 'color', 'red');
xlabel(labels_expected{1}); ylabel(labels_expected{2});

subplot(subplot_inds{2}{:}); 
if ~hide_original_in_prediction
    qhs = show_time_sq_as_quiver_impl(all_image_data, order_inds(3), ...
        order_inds(4), options); % LP vs. RP
    set(qhs, 'color', 'blue');
end
% stack the predicted data
qhs = show_time_sq_as_quiver_prediction_impl(all_image_data, ...
    predicted_data, order_inds(3), order_inds(4), options); % LP vs. RP
set(qhs, 'color', 'red');
xlabel(labels_expected{3}); ylabel(labels_expected{4});

end

function qhs = show_time_sq_as_quiver_prediction_impl(all_image_data, ...
    predicted_data, x_ind, y_ind, options)
qhs = [];

num_subjects = length(all_image_data);
for i = 1:length(all_image_data)
    Ts(i) = size(all_image_data{i},2);
end
T_max = max(Ts);

for t = 1:T_max
    lc = [];
    rc = [];
    for i = 1:num_subjects
        if size(all_image_data{i},2) >= t+1
            lc1 = all_image_data{i}(x_ind,t);
            lc1 = [lc1, predicted_data{i}(x_ind,t+1)];
            lc = cat(1, lc, lc1);
            
            rc1 = all_image_data{i}(y_ind,t);
            rc1 = [rc1, predicted_data{i}(y_ind,t+1)];
            rc = cat(1, rc, rc1);
        end
    end
    if ~isempty(lc) && ~isempty(rc)
        hold on;
        qhs(t) = quiver(lc(:,1), rc(:,1), lc(:,2)-lc(:,1), rc(:,2)-rc(:,1), 0, ...
            'AutoScale','off','linewidth',1);
        set_axes_range(all_image_data, options);
    end
end

end


%% show quivers with weights 
function qhs = show_time_sq_as_quiver_impl_with_weights(all_image_data, ...
    x_ind, y_ind, weights, cmap, options)
qhs = [];

num_subjects = length(all_image_data);
for i = 1:length(all_image_data)
    Ts(i) = size(all_image_data{i},2);
end
T_max = max(Ts);

num_colors = size(cmap,1);

for t = 1 : T_max - 1
    lc = [];
    rc = [];
    color_inds = [];
    for i = 1:num_subjects
        if size(all_image_data{i},2) >= t+1
            lc = cat(1, lc, all_image_data{i}(x_ind,t:t+1));
            rc = cat(1, rc, all_image_data{i}(y_ind,t:t+1));
            color_inds = cat(1, color_inds, round(weights{i}(t) * num_colors));
        end
    end
    color_inds(color_inds > num_colors) = num_colors;
    
    if ~isempty(lc) && ~isempty(rc)
        hold on;
        for i = 1:size(lc,1)
            qhs(end+1) = quiver(lc(i,1), rc(i,1), lc(i,2)-lc(i,1), rc(i,2)-rc(i,1), 0, ...
                'AutoScale','off','linewidth',1,'MaxHeadSize',0.5, ...
                'color', cmap(color_inds(i),:));
        end
        set_axes_range(all_image_data, options);
    end
end
draw_ref_line();

end