function datastruct = importdata2( fname )
%function datastruct = importdata2( fname )
%   imports data using impordata then adds fields using colheaders

datastruct = importdata(fname);
if (~isfield(datastruct, 'colheaders'))
    return;
end

%check for rows with only valid number = frameNum and eliminate them
%these rows have no metadata

whichcol = find(strcmpi (datastruct.colheaders, 'frameNum'), 1, 'first');
inds = setdiff(1:length(datastruct.colheaders), whichcol);
valid = any(datastruct.data(:,inds), 2);
datastruct.data = datastruct.data(valid,:);

for j = 1:length(datastruct.colheaders)
    ch = datastruct.colheaders{j};
    ch(isspace(ch)) = '_';
     datastruct.(ch) = datastruct.data(:,j);
end

end

