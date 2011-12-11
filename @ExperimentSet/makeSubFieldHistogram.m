function h = makeSubFieldHistogram(eset, field, subfield, fieldaxis, varargin)
% generates a histogram of all values in eset.expt.track.(field).(subfield)
% function h = makeSubFieldHistogram(eset, field, subfield, fieldaxis, varargin)
%
% generates a histogram of all values in eset.expt.track.(field).(subfield)
% if no return arguments are specified, generates a plot of that histogram
%
% outputs:
%   H: (optional), the histogram;  if no output arguments, plots histogram
% inputs:
%   ESET: a member of the ExperimentSet class
%   FIELD, SUBFIELD: make the histogram over track.FIELD.SUBFIELD 
%   FIELDAXIS: bin centers for the histogram
%   VARARGIN: 
%      optional parameter/value pairs
%      'r2d',true - when plotting, plot fieldaxis in degrees instead of radians
%      'polar',true - adjust data so that the 0 degree bin is the same size as the
%           next bin;  the last bin may be screwed by this, so be careful
%           in other words, passing polar, fieldaxis = deg2rad(0:30:330) is good
%           passing polar,true deg2rad(0:30:360) will nuke your last bin
r2d = false;
polar = false;
varargin = assignApplicable(varargin);

qv = eset.gatherSubField(field, subfield);
if (polar)
    c1 = fieldaxis(1);
    dt = fieldaxis(2) - c1;
    edge = c1 - dt/2;
    qv = mod(qv-edge,2*pi)+edge;
end
h1 = hist(qv, fieldaxis);

if (nargout == 0)
    if (r2d)
        fieldaxis = rad2deg(fieldaxis);
    end
    bar (fieldaxis, h1); title (eset.defaultTitle);
    xlabel([field '.' subfield]); ylabel('counts'); embiggen();
else
    h = h1;
end


