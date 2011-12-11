function qv = getDerivedQuantity (tp, field, varargin)
%function qv = getDerivedQuantity (tp, field, varargin)
%
%gets derived quantity field from tp.track, subject to the following
%parameter/value pairs
%'position', pos
%pos = 'first' or 'start' -- value at tp.startInd
%pos = 'last' or 'end' -- value at tp.endInd
%pos = fieldname -- value at tp.fieldname (e.g. fieldname = startInd,
%                       centralInd)
%pos = 'atMax' -- value at tp.maxInd : only valid for head swings
%pos = 'mean' -- mean of all values
%pos = 'all' -- all values, this is the default behavior
%pos = 'min' or 'minimum' -- minimum of all values
%pos = 'max' or 'maximum' -- maximum of all values
%pos = 'maxabs' or 'minabs' -- value that is farthest or closest from 0
%
%'reltpfield', rtpf
%relative track part field:  calls tp.(rtpf).getDerivedQuantity(field,
%varargin(:))
%e.g. 'reltpfield', 'prevRun', 'pos', 'end' will get the last point from
%the previous run
%
%'trimpct', f fraction btwn 0 and 0.5;  only take the part of the trackpart
%       between f and 1-f
%'trimpts', npts;  only take the part of the track between npts and
%           end-npts; if this would be no points, then the return is empty
%'operation', op
%applied to quantity before position is applied (i.e. before taking
%min, max, mean, etc.)
%newquantity = op(quantity) 
%examples -- 'operation', 'abs' or 'operation', @abs

if (length(tp) > 1)
    qv = [];
    for j = 1:length(tp)
        
        qv = [qv tp(j).getDerivedQuantity(field, varargin{:})];
        
    end
    return;
end

reltpfield = [];
varargin = assignApplicable(varargin);
if (~isempty(reltpfield))
    qv = tp.(reltpfield).getDerivedQuantity(field, varargin{:});
    return;
end

position = 'all';
operation = [];
trimpct = -1;
trimpts = -1;
varargin = assignApplicable (varargin);

if (tp.startInd > tp.endInd)
%    warning ('TP:GDQ', 'trackpart indices are out of order!');
    temp = tp.startInd;
    tp.startInd = tp.endInd;
    tp.endInd = temp;
    tp.inds = tp.startInd:tp.endInd;
end
if (iscell(position))
    qv = tp.(position{1}).getDerivedQuantity(field, 'position', position{2});
    return;
end

switch (lower(position))
    case {'first','start'}
        inds = tp.startInd;
    case {'last', 'end'}
        inds = tp.endInd;
    case {'atmax'}
        inds = tp.maxInd;
    otherwise
        if (any(strcmp(properties(tp), position)))
            inds = tp.(position);
        else
            inds = tp.inds;
            %inds = tp.startInd:tp.endInd;
        end
end
if (trimpct > 0 && trimpts > 0)
    warning('TP:QGD', 'trackpart.getDerivedQuantity: you cannot ask me to trim both by percent and by points');
end
if (trimpct > 0)
    if (trimpct > 0.5)
        warning('TP:QGD', 'trackpart.getDerivedQuantity: trimpct must be < 0.5');
    else
        npts = tp.endInd - tp.startInd;
        inds = floor(tp.startInd + trimpct * npts) : ceil(tp.endInd - trimpct*npts);
    end
end
if (trimpts > 0)
    inds = (tp.startInd + trimpts):(tp.endInd - trimpts);
end

if (~isempty(tp.track) && ~isempty(inds))
    tr = tp.track;
    tr.calculateDerivedQuantity({'eti', field});
    inds = inds(inds >= 1 & inds <= length(tp.track.dq.eti));
    qv = tr.getDerivedQuantity(field, false, inds);
else
%    warning('TP:GDQ', 'trackpart.getDerivedQuantity: trackpart has no track or empty indices');
 %   tp
%    tp.track
    qv = [];
    return;
end
if (~isempty(operation))
    if (ischar(operation))
        operation = str2func(operation);
    end
    if (isa (operation, 'function_handle'))
        qv = operation(qv);
    end
end

switch (lower(position))
    case 'mean'
        qv = mean(qv, 2);
    case {'min', 'minimum'}
        qv = min(qv,[], 2);
    case {'max', 'maximum'}
        qv = max(qv,[],2);
    case {'maxabs', 'minabs'}
        if (size(qv,1) > 1)
            qv2 = sum(qv.^2);
        else
            qv2 = abs(qv);
        end
        if (strcmpi(position, 'maxabs'))
            [~,I] = max(qv2);
        else
            [~,I] = min(qv2);
        end
        qv = qv(:,I);
end
end

        
