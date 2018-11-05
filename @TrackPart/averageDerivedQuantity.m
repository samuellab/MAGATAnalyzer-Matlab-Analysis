function [qv, datamatrix] = averageDerivedQuantity (tp, field, centerpos, offsetinds, varargin)
%function [qv, datamatrix] = averageDerivedQuantity (tp, field, centerpos, offsetinds, varargin)
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
%'posoffset', m
% gets the value at the index offset by m points from position
% (just adds m to inds); m can be positive or negative
% if posoffset causes ind to be outside track values, the inds are coerced to the nearest valid value 
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

datamatrix = NaN(length(tp), length(offsetinds));
for j = 1:length(tp)
    fv = tp(j).getDerivedQuantity(field, 'position', centerpos, 'posoffset', offsetinds, varargin{:});
    if (length(fv) < size(datamatrix,2))
        inds = tp(j).getDerivedQuantity('validinds', 'position', centerpos, 'posoffset', offsetinds, varargin{:});
        datamatrix(j,inds) = fv;
    else
        datamatrix(j,:) = fv;
    end
    datamatrix2 = datamatrix;
    datamatrix2(~isfinite(datamatrix)) = 0; %replace unknown values with 0
end
qv =  sum(datamatrix2)./sum(isfinite(datamatrix)); %average is sum of known values / number of known values
    

        
