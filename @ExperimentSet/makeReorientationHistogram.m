function [h,eb] = makeReorientationHistogram(eset, fieldname, fieldaxis, varargin)
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
%   'vsDistance' - find reorientation rate in inverse length instead of
%   inverse time
%   'timerange' - only include this time range for analysis             
%   'useprevdir' [true] -- if field is 'theta', use reorientation.prevDir,
%   instead of reorientation.getDerivedQuantity('theta', 'start')
%   'datarow', if field is multi-dimensional (e.g. 'vel'), which row to use
ts1 = tic;
minHS = 0;
maxHS = 10000;
validname = [];
validoperation = @(x) logical(setNonFiniteToZero(x));
r2d = false;
polar = false;
incllastrun = false;
makePlot = false;
turnsequence = [];
vsDistance = false;
timerange = [];
useprevdir = true;
datarow = 1;
varargin = assignApplicable(varargin);

useprevdir = useprevdir && strcmpi(fieldname, 'theta');

if (isstr(validoperation))
    validoperation = str2func(validoperation);
end
if (isa(eset.expt(1).track(1), 'MaggotTrack'))
    hsexpression = ['[track.reorientation.numHS] >= ' num2str(minHS) ' & [track.reorientation.numHS] <= ' num2str(maxHS)];
    if (~useprevdir)
        qv = eset.gatherFromSubField('reorientation', fieldname, 'position', 'start', 'indsExpression', hsexpression);
    else
        nhs = eset.gatherSubField('reorientation', 'numHS');
        qv = eset.gatherSubField('reorientation', 'prevDir');
        qv = qv(:, nhs >= minHS & nhs <= maxHS);
    end
else
    hsexpression = [];
    if (~useprevdir)
        qv = eset.gatherFromSubField('reorientation', fieldname, 'position', 'start');
    else
        qv = eset.gatherSubField('reorientation', 'prevDir');
    end
end

if ~isempty(turnsequence)
    %r = eset.gatherField('reorientation');
    %goodreo = r.turnsequenceEquals(turnsequence);
    goodreo = eset.evaluateTrackExpression(['track.reorientation.turnsequenceEquals(' num2str(turnsequence) ')']);
    while(iscell(goodreo))
        goodreo = [goodreo{:}];
    end
else
    goodreo = true([1 length(qv)]);
end

if ~isempty(timerange)
    if (~isempty(hsexpression))
        reti = eset.gatherFromSubField('reorientation', 'eti', 'position', 'start', 'indsExpression', hsexpression);
    else
        reti = eset.gatherFromSubField('reorientation', 'eti', 'position', 'start');
    end
    validtim = reti >= min(timerange) & reti <= max(timerange);
else
    validtim =  true([1 length(qv)]);
end

valid = true([1 length(qv)]);

if (~isempty(validname))
    if (~iscell(validname))
        validname = {validname};
    end
    if (~iscell(validoperation))
        validoperation = {validoperation};
    end
    for k = 1:length(validname)
        opnum = min(k, length(validoperation));
        if (~isempty(hsexpression))
            valid = valid & validoperation{opnum}(eset.gatherFromSubField('reorientation', validname{k}, 'position', 'start', 'indsExpression', hsexpression));
        else
            valid = valid & validoperation{opnum}(eset.gatherFromSubField('reorientation', validname{k}, 'position', 'start'));
        end
    end
end

qv = qv(:,valid&goodreo&validtim);
if (polar)
    c1 = fieldaxis(1);
    dt = fieldaxis(2) - c1;
    edge = c1 - dt/2;
    qv = mod(qv-edge,2*pi)+edge;
end
h1 = histc(qv(datarow,:), binEdgesFromCenters(fieldaxis));
h1 = h1(1:end-1);
if (incllastrun)
    indsExpression = 'true(size(track.run))';
else
    indsExpression = '1:length(track.run) ~= length(track.run)';
end
qv = eset.gatherFromSubField('run', fieldname, 'indsExpression', indsExpression);
sp = eset.gatherFromSubField('run', 'speed', 'indsExpression', indsExpression);

valid = true([1 length(qv)]);
if (~isempty(validname))
     for k = 1:length(validname)
        opnum = min(k, length(validoperation));
        valid = valid & validoperation{opnum}(eset.gatherFromSubField('run', validname{k}, 'indsExpression', indsExpression));
     end
end

if (~isempty(timerange))
    reti = eset.gatherFromSubField('run', 'eti', 'indsExpression', indsExpression);
    validtim = reti >= min(timerange) & reti <= max(timerange);
else
    validtim =  true([1 length(qv)]);
end


qv = qv(:,valid&validtim);
sp = sp(:,valid&validtim);

if (polar)
    c1 = fieldaxis(1);
    dt = fieldaxis(2) - c1;
    edge = c1 - dt/2;
    qv = mod(qv-edge,2*pi)+edge;
end
%toc
if (~vsDistance)
    h2 = histc(qv(datarow,:), binEdgesFromCenters(fieldaxis)) * eset.expt(1).dr.interpTime / 60;
    h2 = h2(1:end-1);
else
    [~,~,~,~,h2] = meanyvsx (qv(datarow,:), sp,  binEdgesFromCenters(fieldaxis));
    h2 = h2*eset.expt(1).dr.interpTime;
end
%toc

eb = sqrt(h1)./h2;
if (nargout > 0)
    h = h1./h2;
end
if (nargout == 0 || makePlot)
    htemp = h1./h2;
    if (r2d)
        fieldaxis = rad2deg(fieldaxis);
    end
    errorbar (fieldaxis, htemp, eb);
    xlabel (fieldname);
    if (vsDistance)
        ylabel ('reorientation rate (dist^{-1})');
    else
        ylabel ('reorientation rate (min$^{-1}$)');
    end
    title ([eset.defaultTitle ': reorientation rate vs. instantaneous ' fieldname]);
end
%disp (['last line of makeReorientationHistogram - ', num2str(toc(ts1))]);
%{
toc
bob = whos();
nm = {bob.name};
nm = setdiff(nm, {'h'});
clear(nm{:}, 'bob', 'nm');
toc
%}