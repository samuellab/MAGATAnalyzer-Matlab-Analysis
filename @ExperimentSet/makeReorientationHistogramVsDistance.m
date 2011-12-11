function h = makeReorientationHistogramVsDistance(eset, fieldname, fieldaxis, varargin)
% generates a plot of mean reorientation rate (in 1/minutes) vs. fieldname
% function h = makeReorientationHistogram(eset, fieldname, fieldaxis, varargin)
% 
% generates a plot of mean reorientation rate (in 1/minutes) vs. fieldname
% if no return arguments are given, plots this 
% as with almost all eset functions, assumes (explicitly here, implicitly
% other places) that interpTime is the same for all experiments
%
% outputs:
%   H: (optional), the reorientation rate vs. bin center;
%   if no output arguments, plots this
% inputs:
%   ESET: a member of the ExperimentSet class
%   FIELDNAME: the name of the field over which to make the reo rate plot 
%   FIELDAXIS: bin centers for the histogram
%   VARARGIN: optional arguments:
%   'r2d','true' causes plot fieldaxis to be displayed in degrees
%       rather than radians
%   'polar',true - adjust data so that the 0 degree bin is the same size as the
%       next bin;  the last bin may be screwed by this, so be careful
%       in other words, passing polar,true fieldaxis = deg2rad(0:30:330) is good
%       passing polar,true fieldaxis = deg2rad(0:30:360) will nuke your
%       last bin
%   'validname', fieldname
%   'validoperation', op
%       reorientation rate is calculated for the subset of the experiment points
%       that satisfy op(gatherField(fieldname)) == true
%   'incllastrun', true/[false]
%       whether to include the last run of a track when calculating
%       statistics
%   'makePlot', true/[false]
%       if true, make the plot even if a value is returned
%   (maggot tracks only)
%   'minHS',n only include reorientations with at least n headsweeps
%   'maxHS',n only include reorientations with at most n headsweeps
%   (worm tracks only)
%   'turnsequence' - only include reorientations with this sequence of
%   turns
minHS = 0;
maxHS = 10000;
validname = [];
validoperation = @(x) logical(x);
r2d = false;
polar = false;
incllastrun = false;
makePlot = false;
turnsequence = [];
varargin = assignApplicable(varargin);
if (isstr(validoperation))
    validoperation = str2func(validoperation);
end
if (isa(eset.expt(1).track(1), 'MaggotTrack'))
    hsexpression = ['[track.reorientation.numHS] >= ' num2str(minHS) ' & [track.reorientation.numHS] <= ' num2str(maxHS)];
    qv = eset.gatherFromSubField('reorientation', fieldname, 'position', 'start', 'indsExpression', hsexpression);
else
    hsexpression = [];
    qv = eset.gatherFromSubField('reorientation', fieldname, 'position', 'start');
end
if ~isempty(turnsequence)
    r = eset.gatherField('reorientation');
    goodreo = r.turnsequenceEquals(turnsequence);
else
    goodreo = true([1 length(qv)]);
end

if (~isempty(validname))
    if (~isempty(hsexpression))
        valid = validoperation(eset.gatherFromSubField('reorientation', validname, 'position', 'start', 'indsExpression', hsexpression));
    else
        valid = validoperation(eset.gatherFromSubField('reorientation', validname, 'position', 'start'));
    end
else
    valid = true([1 length(qv)]);
end
qv = qv(:,valid&goodreo);
if (polar)
    c1 = fieldaxis(1);
    dt = fieldaxis(2) - c1;
    edge = c1 - dt/2;
    qv = mod(qv-edge,2*pi)+edge;
end
h1 = hist(qv, fieldaxis);

if (incllastrun)
    indsExpression = 'true(size(track.run))';
else
    indsExpression = '1:length(track.run) ~= length(track.run)';
end
qv = eset.gatherFromSubField('run', fieldname, 'indsExpression', indsExpression);
if (~isempty(validname))
    valid = validoperation(eset.gatherFromSubField('run', validname, 'indsExpression', indsExpression));
    qv = qv(:,valid);
end
if (polar)
    c1 = fieldaxis(1);
    dt = fieldaxis(2) - c1;
    edge = c1 - dt/2;
    qv = mod(qv-edge,2*pi)+edge;
end

h2 = hist(qv, fieldaxis);

if (nargout > 0)
    h = h1./h2 /eset.expt(1).dr.interpTime * 60;
end
if (nargout == 0 || makePlot)
    htemp = h1./h2/eset.expt(1).dr.interpTime * 60;
    if (r2d)
        fieldaxis = rad2deg(fieldaxis);
    end
    plot (fieldaxis, htemp);
    xlabel (fieldname);
    ylabel ('reorientation rate (min^{-1})');
    title ([eset.defaultTitle ': reorientation rate vs. instantaneous ' fieldname]);
end
