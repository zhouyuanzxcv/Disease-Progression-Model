function [all_image_data, all_elapsed_days, subject_ids, labels] = ...
    load_data_from_csv()
%LOAD_DATA_FROM_CSV Summary of this function goes here
%   Detailed explanation goes here
labels = {'LC','RC','LP','RP'};
M = read_mixed_csv('./data/DATScan_Analysis.csv',',');
M1 = M(2:end,:); % remove header
patnos = M1(1:end,1);

% sort patient no.
patnos = cellarray2num(patnos);
[patnos_sorted, ind] = sort(patnos);
M2 = M1(ind,:); % sorted rows

% group image data
all_image_data = {};
all_elapsed_days = {};
subject_ids = {};

image_data = [];
last_patno = '';

for i = 1:size(M2,1)
    current_patno = M2{i,1}(2:end-1);
    if ~strcmp(last_patno, current_patno)
        % save old patient
        if ~isempty(image_data)
            all_image_data{end+1} = image_data;
            all_elapsed_days{end+1} = elapsed_days;
            subject_ids{end+1} = last_patno;
        end
        
        % start new patient
        image_data = read_ratios(M2(i,3:6));
        elapsed_days = read_elapsed_days(M2(i,2));
    else
        image_data = [image_data, read_ratios(M2(i,3:6))];
        elapsed_days = [elapsed_days, read_elapsed_days(M2(i,2))];
    end
    
    last_patno = current_patno;
end

% sort elapsed time for each subject
for i = 1:length(all_image_data)
    % remove data with unknown date
    elapsed_days = all_elapsed_days{i};
    image_data = all_image_data{i};
    ind_remove = find(isnan(elapsed_days));
    elapsed_days(ind_remove) = [];
    image_data(:,ind_remove) = [];
    
    % sort time
    [elapsed_days1, ind] = sort(elapsed_days);
    image_data1 = image_data(:,ind);
    if isempty(elapsed_days1)
        all_elapsed_days{i} = elapsed_days1;
    else
        all_elapsed_days{i} = elapsed_days1 - elapsed_days1(1);
    end
    all_image_data{i} = image_data1;
end

% remove sequences with length 1
remove_ind = [];
for i = 1:length(all_image_data)
    if size(all_image_data{i}, 2) <= 1
        remove_ind(end+1) = i;
    end
end

all_image_data(remove_ind) = [];
all_elapsed_days(remove_ind) = [];
subject_ids(remove_ind) = [];

end

function year = convert_event_id_to_year(event_id)
% get year stamp from the ppmi website
% http://www.ppmi-info.org/wp-content/uploads/2010/10/PPMI-Schedule-of-Activities.pdf
switch event_id
    case 'SC'
        year = 0;
    case 'V04'
        year = 1;
    case 'V06'
        year = 2;
    case 'V08'
        year = 3;
    case 'V10'
        year = 4;
    case {'U01','U02','ST'} % unscheduled visit or symptomatic therapy
        year = NaN; 
    otherwise
        fprintf('unsupported event id %s\n', event_id);
        year = NaN;
end

end

function image_data = read_ratios(row)
image_data = cellarray2num(row);
image_data = image_data([2,1,4,3]);
image_data = image_data';

end

function days = read_elapsed_days(cell)
event_id = cell{1}(2:end-1);
year = convert_event_id_to_year(event_id);
days = year * 365;

end

function Y = cellarray2num(X)
Y = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Y(i,j) = str2num(X{i,j}(2:end-1)); % removing leading and trailing "
    end
end
end