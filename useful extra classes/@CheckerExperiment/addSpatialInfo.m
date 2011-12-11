function addSpatialInfo (expt, varargin)
%function addSpatialInfo (expt, varargin)
%
%override any of these fields with key,value pairs
%
%structvar; %There is no god but structvar and structvar is his prophe
%boundaryImage; 
%lightImage;
%imsize;       
%boundaryThickness;

structvar = expt.structvar;
dirim = [];
distim = [];
lightim = [];
imsize = [];
boundaryThickness = [];

varargin = assignApplicable(varargin);


if (isempty(dirim) || isempty(distim) || isempty(lightim))
    if (isempty('structvar'))
        disp('I need structvar; feed me structvar');
        return
    end
    existsAndDefault('imsize', [1944 2592]);
    existsAndDefault('boundaryThickness', expt.boundaryThickness);
    [dirim,distim,lightim] = CheckerExperiment.generateBoundaryImage(imsize, structvar);
    expt.boundaryThickness = boundaryThickness;
    expt.structvar = structvar;
end

expt.imsize = size(dirim);
    
expt.addBoundary(dirim,distim,lightim);

expt.calculateDerivedQuantity('theta');
expt.evaluateTrackExpression('track.dq.thetaToBound = diff(unwrap([track.dq.boundarytheta;track.dq.theta]));');

[xim,yim,igrid,jgrid] = structvar2grid(structvar);

gq = GlobalQuantity;
gq.derivationMethod = @GlobalQuantity.tri2Dinterpolation;
gq.xData = [xim';yim'];
gq.yData = [igrid';jgrid'];

gq.xField = 'sloc';
gq.fieldname = 'gridLoc';
expt.addGlobalQuantity(gq);

gq.xField = 'shead';
gq.fieldname = 'gridHead';
expt.addGlobalQuantity(gq);

