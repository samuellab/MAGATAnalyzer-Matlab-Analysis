function ind = indexCrossingBox (ptlist, index, range, position)
%function ind = indexCrossingBox (ptlist, index, range, position)
%
%finds either the beginning of the track given in ptlist (position = start)
%or the end of the track that goes through ptlist(index) and remains
%entirely within the box defined by ptlist(index) +/- range/2 
%if range is 1D, it is used for both

st = strcmpi(position, 'start');
if (st)
    ptlist = ptlist(:,1:index);
else 
    ptlist = ptlist(:,index:end);
end
range = abs(range);
if size(range,1) == 1
    range = range';
end
sz = size(ptlist);
sz(1) = 1;
if (size(range == 1))
    range = repmat (range, size(ptlist));
else
    range = repmat (range, sz);
end
center = repmat(ptlist(:,index), sz);
inrange = all(ptlist < (center + range/2), 1) & all(ptlist > (center - range/2), 1);
if (st)
    ind = find(~inrange, 1, 'last');
    if (isempty(ind))
        ind = 1;
    end
else
    ind = find(~inrange, 1, 'first');
    if (isempty(ind))
        ind = length(ptlist);
    else
        ind = index + ind;
    end
end

