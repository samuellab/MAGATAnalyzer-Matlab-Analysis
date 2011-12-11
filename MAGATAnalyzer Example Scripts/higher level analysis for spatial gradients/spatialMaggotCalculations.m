function sc = spatialMaggotCalculations (eset, spatial_navigation_options)
%function sc = spatialMaggotCalculations (eset,
%spatial_navigation_options)
%
% spatial_navigation_options -- 1 set of options for all esets


sno.angleBinSize = 30; % in degrees
sno.dtAngleBinSize = 20; % in degrees
sno.hsdtAngleBinSize = 20; % in degrees
sno.preferredDirection = 0;
sno.reoBinSize = 30;
sno.hsBinSize = 90;
sno.hsBinSpacing = 90;
sno.minHS = 1;
sno.minHSTheta = 20;
sno.validname = [];
sno.relativeDirField = [];
sno.dirOffsetField = [];
sno.validoperation = func2str(@(x) logical(setNonFiniteToZero(x)));
sno.confidenceLevel = 0.95;
sno.autocorr_timerange = [];
sno.runTimeBinSize = 10;

if (nargin == 0)
    if (nargout >= 1)
        sc = sno;
    end
    %disp ('spatialNavigationMaggotAnalysis (esets, spatial_navigation_options)');
    return;
end
existsAndDefault('spatial_navigation_options', []);

if (length(eset) > 1)
    for j = 1:length(eset)
        disp (['start ' num2str(j)]); ts1 = tic;
        sc(j) = spatialMaggotCalculations_Janelia(eset(j), spatial_navigation_options);  %#ok<AGROW>
        disp (['end ' num2str(j) ' ' num2str(toc(ts1)) ' s']);
    end
    return;
end

if (~isempty(spatial_navigation_options) && isstruct(spatial_navigation_options))
    fn = fieldnames(spatial_navigation_options);
    for j = 1:length(fn)
        sno.(fn{j}) = spatial_navigation_options.(fn{j});
    end
end

if (xor(isempty(sno.relativeDirField), isempty(sno.dirOffsetField)))
    disp('relativeDirField = theta - dirOffsetField');
    error('you must define both relativeDirField and dirOffsetField or neither');
    return; %#ok<UNRCH>
end

%general statistics of experiment, may be duplicative
sc.eset_stats = calculateStatisticsOfEset(eset, 'validname', sno.validname, 'validoperation', sno.validoperation);
sc.autocorr_tau = eset.getAutocorrTau;
sc.interptime = median(eset.gatherSubField('dr', 'interpTime'));
sc.sno = sno;
if (ischar(sno.validoperation))
    sno.validoperation = str2func(sno.validoperation);
end
% 
% ad.txf = -180:1:180;
% ad.tx = (-180:sno.angleBinSize:(180-sno.angleBinSize)) + sno.preferredDirection;
% ad.txc = (-180:sno.angleBinSize:(180)) + sno.preferredDirection;
% 
% txrad = deg2rad(ad.txf);
% bsa = deg2rad(sno.angleBinSize);
% 

%{

dtx = (-180:sno.dtAngleBinSize:(180-sno.dtAngleBinSize));
dtxc = (-180:sno.dtAngleBinSize:180);

hsdtx = sort([(sno.minHSTheta - sno.hsdtAngleBinSize/2):(-sno.hsdtAngleBinSize):0 (sno.minHSTheta + sno.hsdtAngleBinSize/2):(sno.hsdtAngleBinSize):180]);
hsdtx = sort([-hsdtx hsdtx]);
%hsdtx = hsdtxc(1:(end-1));
%}

%calculate statistics of runs including change within runs and run length
%-------------------------------------------------------------------------%
rs = eset.gatherSubField('run', 'startTheta');
re = eset.gatherSubField('run', 'endTheta');
rt = eset.gatherSubField('run', 'runTime');
rmt = eset.gatherSubField('run', 'meanTheta');

dt = diff(unwrap([rs;re]));

sc.runStartTheta = rs;
sc.runEndTheta = re;
sc.runTime = rt;
sc.runMeanTheta = rmt;
sc.runDeltaTheta = dt;
sc.runLength = eset.gatherSubField('run', 'pathLength');
sc.runEuclidLength = eset.gatherSubField('run', 'euclidLength');
if (~isempty(sno.validname))
    sc.runvalidField = (eset.gatherFromSubField('run', sno.validname, 'position', 'start'));
    sc.runvalid = sno.validoperation(sc.runvalidField);
else
    sc.runvalidField = [];
    sc.runvalid = true(size(sc.runTime));
end


sc.minruntime = median(eset.gatherSubField('so', 'minRunTime'));


%calculate instantaneous change within runs vs. heading angle
%-----------------------------------------------------------------------%
if(~isempty(sno.relativeDirField))
    fieldname = sno.relativeDirField;
else
    fieldname = 'theta';
end

sc.dir = eset.gatherField(fieldname); %dir is either theta or gather field
sc.isrun = eset.gatherField('isrun');
if (~isempty(sno.validname))
    sc.validField = (eset.gatherField(sno.validname));
    sc.valid = sno.validoperation(sc.validField);
else
    sc.validField = [];
    sc.valid = true(size(sc.isrun));
end

for j = 1:length(eset.expt)
    sc.exp(j).fname = eset.expt(j).fname;
    sc.exp(j).vel = eset.expt(j).gatherField('vel'); %dir is either theta or gather field
    sc.exp(j).isrun = eset.expt(j).gatherField('isrun');
    if (~isempty(sno.validname))
        sc.exp(j).validField = (eset.expt(j).gatherField(sno.validname));
        sc.exp(j).valid = sno.validoperation(sc.exp(j).validField);
    else
        sc.exp(j).validField = [];
        sc.exp(j).valid = true(size(sc.exp(j).isrun));
    end
end

sc.lrdtheta = eset.gatherField('lrdtheta');
sc.speed = eset.gatherField('speed');
sc.vel = eset.gatherField('vel');
sc.numhs = eset.gatherSubField('reorientation', 'numHS');

if (~isempty(sno.validname))
    sc.reovalidField = (eset.gatherFromSubField('reorientation', sno.validname, 'position', 'start'));
    sc.reo_prevrunvalidField = (eset.gatherFromSubField('reorientation', sno.validname, 'position', {'prevRun', 'start'}));
    sc.reovalid = sno.validoperation(sc.reovalidField);
    sc.reo_prevrunvalid = sno.validoperation(sc.reo_prevrunvalidField);
else
    sc.reovalidField = [];
    sc.reovalid = true(size(sc.numhs));
    sc.reo_prevrunvalidField = [];
    sc.reo_prevrunvalid = sc.reovalid;
end
if (~isempty(sno.dirOffsetField))
    sc.reo_dirOffset = eset.gatherFromSubField('reorientation', sno.dirOffsetField, 'position', 'mean');
else
    sc.reo_dirOffset = zeros(size(sc.numhs));
end
sc.reo_nextdir = eset.gatherSubField('reorientation', 'nextDir');
sc.reo_prevdir = eset.gatherSubField('reorientation', 'prevDir');
sc.reo_dtheta = diff(unwrap([sc.reo_prevdir;sc.reo_nextdir]));



sc.hs_taildir = eset.gatherSubField('headSwing', 'tailDir');
if (~isempty(sno.dirOffsetField))
    sc.hs_dirOffset = eset.gatherFromSubField('headSwing', sno.dirOffsetField, 'position', 'mean');
else
    sc.hs_dirOffset = zeros(size(sc.hs_taildir));
end

sc.hs_headdir = eset.gatherSubField('headSwing', 'headDir');
sc.hs_acc = eset.gatherSubField('headSwing', 'accepted');
sc.hs_sign = eset.gatherSubField('headSwing', 'sign');
sc.hs_htv = logical(eset.gatherSubField('headSwing', 'valid'));
sc.hs_maxtheta = eset.gatherSubField('headSwing', 'maxTheta');
if (~isempty(sno.validname))
    sc.hsvalidField = (eset.gatherFromSubField('headSwing', sno.validname, 'position', 'start'));
    sc.hsvalid = sno.validoperation(sc.hsvalidField);
else
    sc.hsvalidField = [];
    sc.hsvalid = true(size(sc.hs_htv));
end
sc.hsnum = eset.gatherSubField('headSwing', 'num');

eti = eset.gatherField('eti');

%find the time with the most valid animals
tx = min(eti(sc.valid)):30:max(eti(sc.valid));
h = hist(eti, tx);
[~,I] = max(h);


sl  = eset.gatherField('spineLength', 'mean');
st = eset.evaluateTrackExpression('min(track.getDerivedQuantity(''eti''))');
et = eset.evaluateTrackExpression('max(track.getDerivedQuantity(''eti''))');
sc.spineLengths = sl(st <= tx(I) & et > tx(I) & et-st > 30);



% 
