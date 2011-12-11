function h2 = make2DReorientationHistogram(eset, fieldname1, fieldaxis1, fieldname2, fieldaxis2, varargin)
%function h2 = make2DReorientationHistogram(eset, fieldname1, fieldaxis1, fieldname2, fieldaxis2, varargin)
%
%generates a plot of mean reorientation rate (in 1/minutes) vs. fieldname
%if no return arguments are given, plots this 
%passing 'r2d','true' causes plot fieldaxis to be displayed in degrees
%rather than radians
%
%as will almost all eset functions, assumes (explicitly here, implicitly
%other places) that interpTime is the same for all experiments
%
%optional params:
%'validname', fieldname
%'validoperation', op
%reorientation rate is calculated for the subset of the experiment points
%that satisfy op(gatherField(fieldname)) == true
%
%'minHS',n only include reorientations with at least n headsweeps
%'maxHS',n only include reorientations with at most n headsweeps

minHS = 0;
maxHS = 10000;
validname = [];
validoperation = @(x) logical(x);
r2d = false;
incllastrun = false;
varargin = assignApplicable(varargin);
if (isstr(validoperation))
    validoperation = str2func(validoperation);
end

hsexpression = ['[track.reorientation.numHS] >= ' num2str(minHS) ' & [track.reorientation.numHS] <= ' num2str(maxHS)];

qv1 = eset.gatherFromSubField('reorientation', fieldname1, 'position', 'start', 'indsExpression', hsexpression);
qv2 = eset.gatherFromSubField('reorientation', fieldname2, 'position', 'start', 'indsExpression', hsexpression);

if (~isempty(validname))
    valid = validoperation(eset.gatherFromSubField('reorientation', validname, 'position', 'start', 'indsExpression', hsexpression));
    qv = qv1(:,valid);
    qv = qv2(:,valid);
end

h = makeIm(qv1, qv2, fieldaxis1, fieldaxis2);

if (incllastrun)
    indsExpression = 'true(size(track.run))';
else
    indsExpression = '1:length(track.run) ~= length(track.run)';
end
qv1 = eset.gatherFromSubField('run', fieldname1, 'indsExpression', indsExpression);
qv2 = eset.gatherFromSubField('run', fieldname2, 'indsExpression', indsExpression);

if (~isempty(validname))
    valid = validoperation(eset.gatherFromSubField('run', validname, 'indsExpression', indsExpression));
    qv1 = qv1(:,valid);
    qv2 = qv2(:,valid);
end

hrun = makeIm(qv1, qv2, fieldaxis1, fieldaxis2);


if (nargout > 0)
    h2 = h./hrun /eset.expt(1).dr.interpTime * 60;
else
    htemp = h./hrun/eset.expt(1).dr.interpTime * 60;
    pcolor (fieldaxis1, fieldaxis2, htemp); shading flat
    colorbar vert
    xlabel (fieldname1);
    ylabel (fieldname2);
    title ([eset.defaultTitle ': reorientation rate vs. instantaneous ' fieldname1 '&' fieldname2]);
end
