function readings = importTHLReadings( input_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dd = dir(fullfile(input_dir, '*.csv'));
if (isempty(dd))
    warning ('thl:emptydir', ['no csv files found in ' input_dir]);
    readings = [];
    return;
end
fn = {};
for j = 1:length(dd)
    data{j} = importTHLData(fullfile(input_dir, dd(j).name)); %#ok<*AGROW>
    fn = union(fn, fieldnames(data{j}));
end

fn = setdiff(fn, {'data', 'textdata'});
for j = 1:length(data)
    for k = 1:length(fn)
        if (~any(cellfun(@(s) ~isempty(strfind(s, fn{k})), fieldnames(data{j}))))
            warning ('thl:notfound', ['field ' fn{k} ' not present in file: ' dd(j).name]);
        end
    end
end

for j = 1:length(data)
    fn = intersect(fn, fieldnames(data{j}));
end
for j = 1:length(data)
    data{j} = rmfield(data{j}, setdiff(fieldnames(data{j}),fn));
end
data = cell2mat(data);



for j = 1:length(fn)
    readings.(fn{j}) = [data.(fn{j})];
end
if (~isfield(readings, 'time'))
    warning ('thl:notime', 'lost time information');
    return;
end

[~,I] = sort(readings.time, 'ascend');
for j = 1:length(fn)
    readings.(fn{j}) = readings.(fn{j})(I);
end

