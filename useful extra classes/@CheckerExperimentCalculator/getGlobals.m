function gq = getGlobals(cec, varargin)
%calculates global quantities to add to experiments based on checkerboard
%layout and cec settings
%
% optional args:
% locfield
% prefix
% example use:  
%   cec.globalQuantities = [cec.getGlobals() cec.getGlobals('locfield', 'shead', 'prefix', 'head_')];
% to get border quanties vs. position and head position
% note that thetaToBorder is only calculated if prefix is empty

locfield = 'sloc';
prefix = ''; %use prefix in combination with locfield to do e.g. head distance to border
varargin = assignApplicable(varargin);

borderpts = cec.borderSizeInCm/cec.imageResolution;
cornerpts = borderpts*1.5;


%do logical fields
gq = GlobalQuantity();
gq.xField = locfield;
gq.derivationMethod = @GlobalQuantity.interpLogicalIm;

aifields = {'bwim', 'border', 'oncorner', 'interior'};
gqnames = {'inlight', 'onborder', 'oncorner', 'interior'};
gq = repmat(gq, [1 length(aifields)]);
for j = 1:length(aifields)
    gq(j).yData = bwpack(cec.ai.(aifields{j}));
    upi = bwunpack(gq(j).yData);
    rx = interp1(cec.rx, 1:size(upi,2), 'linear', 'extrap');
    ry = interp1(cec.ry, 1:size(upi,1), 'linear', 'extrap');
    xdata.x = rx;
    xdata.y = ry;
    gq(j).xData = xdata;
    gq(j).fieldname = [prefix gqnames{j}];
end
gqlog = gq;


maxval = max(cec.ai.distToBorder(:)) - borderpts;
minval = min(cec.ai.distToBorder(:)) + borderpts;

bp = ceil(borderpts);
zmid = [-2*bp:bp/5:-bp -bp:2:bp bp:bp/5:2*bp];
ztop = 2*bp:2*bp:maxval;
zbottom = -2*bp:-2*bp:minval;


z = unique([zbottom ztop zmid]);
c = contourc(double(cec.ai.distToBorder),double(z));

ind = [];
nextind = 1;
while (nextind < length(c))
    ind = [ind, nextind]; %#ok<AGROW>
    nextind = nextind + c(2,nextind) + 1;
end
inds = setdiff(1:length(c), ind);
cpts = round(c(:,inds));
%{
datainds = unique(sub2ind(size(cec.ai.distToDark), cpts(2,:), cpts(1,:)));

[xx,yy] = meshgrid(cec.rx, cec.ry);
edgeinds = (find(xx == min(xx(:)) | xx == max(xx(:)) | yy == min(yy(:)) | yy == max(yy(:))))';

datainds = unique([datainds edgeinds]);
xx = xx(datainds);
yy = yy(datainds);
xData = [xx;yy];
%}
gq = GlobalQuantity();
gq.xField = locfield;
[xx,yy] = meshgrid(cec.rx, cec.ry);
xdata.x = xx;
xdata.y = yy;
gq.xData = xdata;
gq.derivationMethod = @GlobalQuantity.twoDinterpolation;

gq = repmat(gq, 1, 3);

gq(1).fieldname = [prefix 'distToBorder'];
gq(1).yData = cec.ai.distToBorder;

gq(2).fieldname = [prefix 'dirToBorder'];
gq(2).yData = cec.ai.dirToBorder;
gq(2).derivationMethod = @GlobalQuantity.interpAngleIm;
%gq(2).yData = cat(3,cos(cec.ai.dirToBorder), sin(cec.ai.dirToBorder));

gq(3).fieldname = [prefix 'borderNum'];
gq(3).xData = [];
gq(3).yData = [];
onb = [prefix 'onborder'];
dtb = gq(2).fieldname;
gq(3).xField = {onb, dtb};
eval (['gq(3).derivationMethod = @(xin, xdata, ydata) -1*(~xin.' onb ') + xin.' onb '.*(mod(round(mod(xin.' dtb ', 2*pi)*2/pi), 4));']);
n = length(gq);
if (isempty(prefix))
    gq(n+1).xField = {gq(2).fieldname, 'theta'};
    gq(n+1).yData = [];
    gq(n+1).xData = [];
    gq(n+1).derivationMethod = @(xin, xdata, ydata) diff(unwrap([xin.dirToBorder;xin.theta]));
    gq(n+1).fieldname = 'thetaToBorder';
end

gq = [gqlog gq];
