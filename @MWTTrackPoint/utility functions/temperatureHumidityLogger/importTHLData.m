
function data = importTHLData (fn)
    
    data = importdata(fn);
    colheaders = data.textdata(1,:);
    colheaders = regexp(colheaders, '(\w)+', 'match', 'once');
    colheaders = lower(colheaders);
    ind = find(strcmpi(colheaders, 'time'));
    data.(colheaders{ind}) = reshape(datenum(data.textdata(2:end,ind)),1,[]);
    colheaders = colheaders(~(strcmpi(colheaders, 'time')));
    for j = 1:length(colheaders)
        data.(colheaders{j}) = reshape(data.data(:,j),1,[]);
    end
    