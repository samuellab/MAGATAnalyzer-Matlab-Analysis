function ad = spatialMaggotAnalysis (sc, spatial_navigation_options)
%function ad = spatialMaggotAnalysis (sc, spatial_navigation_options)
%
% spatial_navigation_options -- 1 set of options for all esets
% plot_options -- each eset gets its own options
% plot_options.
%    lineWidth -- width of lines
%    color -- color for line and marker
%    legendEntry -- what to put in legend
%    marker -- marker 
%    plotOptions -- additional options to pass to plotting function

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
sno.useInvalidHeadSweeps = false;
sno.timeBinSizeForTemporalPlots = 15; %seconds

if (nargin == 0)
    if (nargout >= 1)
        ad = sno;
    end
    %disp ('spatialNavigationMaggotAnalysis (spatialcalcs, spatial_navigation_options)');
    return;
end
existsAndDefault('spatial_navigation_options', sc(1).sno);

if (length(sc) > 1)
    for j = 1:length(sc)
        disp (['start ' num2str(j)]); ts1 = tic;
        ad(j) = spatialMaggotAnalysis(sc(j), spatial_navigation_options);  %#ok<AGROW>
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


ad.sno = sno;
ad.eset_stats = sc.eset_stats;
ad.summary_message = {[num2str(ad.eset_stats.numExpts) ' Experiments, ' num2str(round(ad.eset_stats.numAnimals)), ' animals. Total time = ' num2str(ad.eset_stats.animalTime/3600, '%.1f') ...
    ' hours (' num2str(ad.eset_stats.animalTime/ad.eset_stats.numAnimals/60, '%.1f') ' min/animal)'], ...
    [num2str(sum(ad.eset_stats.numRunsFromDirection)) ' runs, ' num2str(sum(ad.eset_stats.numReosFromDirection)), ' reorientations (' num2str(sum(ad.eset_stats.numReosWithHSFromDirection)) ' with at least 1 headsweep)']};


if (ischar(sno.validoperation))
    sno.validoperation = str2func(sno.validoperation);
end

%simple metrics -- 
% average run length to each quadrant
% average run time to each quadrant
% number of runs to each quadrant
% probability of turning towards better direction from perpendicular directions
% average run speed;
% average speed; 
% percent of time in runs
% mean number of headsweeps/reorientation
% fraction of headsweeps marked as valid; 
% mean and std of spinelength

sm.quadrants = sno.preferredDirection + (0:90:270);
for j = 1:4
    inds = cos(sc.runMeanTheta - deg2rad(sm.quadrants(j))) > 1/sqrt(2) & sc.runvalid;
    sm.numRuns(j) = nnz(inds);
    sm.runLen(j) = mean(sc.runLength(inds));
    sm.runLen_eb(j) = std(sc.runLength(inds))/sqrt(sm.numRuns(j));
    sm.runTime(j) = mean(sc.runTime(inds));
    sm.runTime_eb(j) = std(sc.runTime(inds))/sqrt(sm.numRuns(j));
end

inds = sc.reovalid & sc.numhs >= sno.minHS & abs(cos(sc.reo_prevdir - deg2rad(sno.preferredDirection))) < 1/sqrt(2);
pd = sc.reo_prevdir(inds) - deg2rad(sno.preferredDirection);
nd = sc.reo_nextdir(inds) - deg2rad(sno.preferredDirection);
towards = cos(pd) < cos(nd);
sm.prob_turn_towards = mean(towards);
x = mean(towards);
n = length(towards);
sm.prob_turn_towards_eb = sqrt(x*(1-x)/n);

sm.meanrunspeed = mean(sc.speed(sc.isrun & sc.valid));
sm.meanspeed = mean(sc.speed(sc.valid));
sm.pctTimeInRuns = mean(sc.isrun(sc.valid));
sm.meanHSPerReo = mean(sc.numhs(sc.reovalid & sc.numhs >= sno.minHS));
sm.fractionOfHSMarkedValid = mean(sc.hs_htv(sc.hsvalid));
sm.meanSpineLength = mean(sc.spineLengths);
sm.stdSpineLength = std(sc.spineLengths);

ad.simpleMetrics = sm;

kappa = sqrt((1 + exp(-sc.interptime/sc.autocorr_tau))/(1 - exp(-sc.interptime/sc.autocorr_tau)));

ad.txf = -180:1:180;
ad.tx = (-180:sno.angleBinSize:(180-sno.angleBinSize)) + sno.preferredDirection;
ad.txc = (-180:sno.angleBinSize:(180)) + sno.preferredDirection;

txrad = deg2rad(ad.txf);
bsa = deg2rad(sno.angleBinSize);

rs = sc.runStartTheta(sc.runvalid);
rt = sc.runTime(sc.runvalid);
dt = diff(unwrap([sc.runStartTheta;sc.runEndTheta]));
dt = dt(sc.runvalid);

[~,my, se] = meanyvsx_slidingwindow(rs, dt, txrad, bsa, 'step', true);
ad.meanrunchange_step = my;
ad.meanrunchange_step_eb = se;
[~,my, se] = meanyvsx_slidingwindow(rs, dt, txrad, bsa, 'gaussian', true);
ad.meanrunchange_gauss = my;
ad.meanrunchange_gauss_eb = se;

[~,my, se] = meanyvsx_slidingwindow(rs, dt.^2, txrad, bsa, 'step', true);
ad.magrunchange_step = my;
ad.magrunchange_step_eb = se;
[~,my, se] = meanyvsx_slidingwindow(rs, dt.^2, txrad, bsa, 'gaussian', true);
ad.magrunchange_gauss = my;
ad.magrunchange_gauss_eb = se;


[~,my, se] = meanyvsx_slidingwindow(rs, dt.^2./rt, txrad, bsa, 'step', true);
ad.diffconstrunchange_step = my;
ad.diffconstrunchange_step_eb = se;
[~,my, se] = meanyvsx_slidingwindow(rs, dt.^2./rt, txrad, bsa, 'gaussian', true);
ad.diffconstrunchange_gauss = my;
ad.diffconstrunchange_gauss_eb = se;





rt = sc.runTime(sc.runvalid);
rmt = sc.runMeanTheta(sc.runvalid);
ad.runTimeAxis = sc.minruntime + (sno.runTimeBinSize/2)+ (0:sno.runTimeBinSize:(sno.runTimeBinSize*ceil(600/sno.runTimeBinSize)));
rth = hist(rt(cos(rmt - deg2rad(sno.preferredDirection)) > 1/sqrt(2)), ad.runTimeAxis);
ad.runTimeHistTowards = (rth/sum(rth));
ad.runTimeHistTowards_eb = sqrt((rth/sum(rth)).*(1-(rth/sum(rth)))./sum(rth));
rth = hist(rt(-cos(rmt - deg2rad(sno.preferredDirection)) > 1/sqrt(2)), ad.runTimeAxis);
ad.runTimeHistAway = (rth/sum(rth));
ad.runTimeHistAway_eb = sqrt((rth/sum(rth)).*(1-(rth/sum(rth)))./sum(rth));


th = sc.dir(sc.isrun & sc.valid);
dth = sc.lrdtheta(sc.isrun & sc.valid);
[~,my,se] = meanyvsx_slidingwindow(th, dth, txrad, bsa, 'step', true);
ad.instantaneousdthetavstheta_step = my;
ad.instantaneousdthetavstheta_step_eb = se*kappa;

[~,my,se] = meanyvsx_slidingwindow(th, dth, txrad, bsa, 'gaussian', true);
ad.instantaneousdthetavstheta_gauss = my;
ad.instantaneousdthetavstheta_gauss_eb = se*kappa;

sp = sc.speed(sc.isrun & sc.valid);
[~,my,se] = meanyvsx_slidingwindow(th, sp, txrad, bsa, 'step', true);
ad.speedVsDir_step = my;
ad.speedVsDir_step_eb = se*kappa;

[~,my,se] = meanyvsx_slidingwindow(th, sp, txrad, bsa, 'gaussian', true);
ad.speedVsDir_gauss = my;
ad.speedVsDir_gauss_eb = se*kappa;

%speedVsDir = my;
%speedVsDir_eb = se;



%calculate navigation index for eset as a whole and for individual
%experiments
%------------------------------------------------------------------------%

v = [sc.exp.vel];
v = v(:,[sc.exp.valid]);
M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
v = M*v;
s = sqrt(sum(v.^2));
ad.navind = mean(v,2)/mean(s);
ad.navind_eb = kappa*std(v,0,2)/mean(s)/sqrt(length(s));

for k = 1:length(sc.exp)
    v = sc.exp(k).vel;
    v = v(:,sc.exp(k).valid);
    M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
    v = M*v;
    s = sqrt(sum(v.^2));
    ad.navind_expt(:,k) = mean(v,2)/mean(s);
    ad.navind_expt_eb(:,k) = kappa*std(v,0,2)/mean(s)/sqrt(length(s));
end


v = [sc.exp.vel];
v = v(:,[sc.exp.valid] & [sc.exp.isrun]);
M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
v = M*v;
s = sqrt(sum(v.^2));
ad.navind_run = mean(v,2)/mean(s);
ad.navind_run_eb = kappa*std(v,0,2)/mean(s)/sqrt(length(s));

for k = 1:length(sc.exp)
    v = sc.exp(k).vel;
    v = v(:,sc.exp(k).valid&sc.exp(k).isrun);
    M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
    v = M*v;
    s = sqrt(sum(v.^2));
    ad.navind_run_expt(:,k) = mean(v,2)/mean(s);
    ad.navind_run_expt_eb(:,k) = kappa*std(v,0,2)/mean(s)/sqrt(length(s));
end

%calculate histogram of instantaneous direction
%------------------------------------------------------------------------%

data = adjustForPolarHistogram(sc.dir(sc.isrun & sc.valid), deg2rad(ad.tx));
h = hist(data(isfinite(data)), deg2rad(ad.tx));
n = sum(h);
eb = kappa*sqrt((h/n).*(1-h/n))*sqrt(n);
ad.thetahist_eb = eb/n;
ad.thetahist = h/n;


%calculate reorientation rate vs. direction
%------------------------------------------------------------------------%

h1 = histc(adjustForPolarHistogram(sc.reo_prevdir(sc.reovalid & sc.numhs >= sno.minHS) - sc.reo_dirOffset(sc.reovalid & sc.numhs >= sno.minHS), deg2rad(ad.tx)),binEdgesFromCenters(deg2rad(ad.tx)));
h2 = histc(adjustForPolarHistogram(sc.dir(sc.isrun & sc.valid), deg2rad(ad.tx)),binEdgesFromCenters(deg2rad(ad.tx)));
h1 = h1(1:(end-1));
h2 = h2(1:(end-1));

ad.reohist = h1./h2 * 60/sc.interptime;
ad.reohist_eb = sqrt(h1)./h2 * 60/sc.interptime;


%calculate statistics of reorientations
%-------------------------------------------------------------------------
ad.reobasedirections = [mod(sno.preferredDirection + 180, 360) - 180, mod(sno.preferredDirection, 360) - 180, mod(sno.preferredDirection + 270, 360) - 180, mod(sno.preferredDirection + 90, 360) - 180];

ad.reotx = (-180:sno.reoBinSize:(180-sno.reoBinSize)) + sno.preferredDirection;
ad.reotxc = (-180:sno.reoBinSize:(180)) + sno.preferredDirection;


reo_prevdir = sc.reo_prevdir(sc.reovalid & sc.numhs >= sno.minHS);
reo_nextdir = sc.reo_nextdir(sc.reovalid & sc.numhs >= sno.minHS);
reo_dtheta = sc.reo_dtheta(sc.reovalid & sc.numhs >= sno.minHS);
nhs = sc.numhs(sc.reovalid & sc.numhs >= sno.minHS);

nd = adjustForPolarHistogram(reo_nextdir, deg2rad(ad.tx));  
runStartDirectionHist = hist(nd, deg2rad(ad.tx));
h = runStartDirectionHist;
n = sum(h);

ad.runStartDirectionHist_eb = sqrt(h/n).*sqrt(1-h/n).*sqrt(n) ./ sum(runStartDirectionHist);
ad.runStartDirectionHist = runStartDirectionHist./sum(runStartDirectionHist);
 
%calculate distributions of reorientations for given previous headings
topref = cos(reo_prevdir - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(reo_prevdir - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(reo_prevdir - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(reo_prevdir - deg2rad(sno.preferredDirection)) < -1/sqrt(2);

indset = {topref, frompref, leftofpref, rightofpref};


ad.dtx = (-180:sno.dtAngleBinSize:(180-sno.dtAngleBinSize));
ad.dtxc = (-180:sno.dtAngleBinSize:180);


dt = adjustForPolarHistogram(reo_dtheta, deg2rad(ad.dtx)); 
for j = 1:length(indset)
    ad.reo_dtheta_dist_hs{j} = hist(dt(indset{j} & nhs > 0), deg2rad(ad.dtx));
    ad.reo_dtheta_dist_hs_eb{j} = sqrt(ad.reo_dtheta_dist_hs{j}.*(1-ad.reo_dtheta_dist_hs{j}/nnz(indset{j} & nhs > 0)));
    ad.reo_dtheta_dist_nohs{j} = hist(dt(indset{j} & nhs == 0), deg2rad(ad.dtx));
    ad.reo_dtheta_dist_nohs_eb{j} = sqrt(ad.reo_dtheta_dist_nohs{j}.*(1-ad.reo_dtheta_dist_nohs{j}/nnz(indset{j} & nhs == 0)));
end

%calculate mean reorientation magnitude, direction, and number of head sweeps vs. previous heading

[~,my,se] = meanyvsx_slidingwindow(reo_prevdir, reo_dtheta.^2,  txrad, deg2rad(sno.reoBinSize), 'step', true);
ad.reoMag_step = my;
ad.reoMag_step_eb = se; 
[~,my,se] = meanyvsx_slidingwindow(reo_prevdir, reo_dtheta.^2,  txrad, deg2rad(sno.reoBinSize), 'gaussian', true);
ad.reoMag_gauss = my;
ad.reoMag_gauss_eb = se; 


[~,my,se] = meanyvsx_slidingwindow(reo_prevdir, reo_dtheta,  txrad, deg2rad(sno.reoBinSize), 'step', true);
ad.reoDir_step = my;
ad.reoDir_step_eb = se; 
[~,my,se] = meanyvsx_slidingwindow(reo_prevdir, reo_dtheta,  txrad, deg2rad(sno.reoBinSize), 'gaussian', true);
ad.reoDir_gauss = my;
ad.reoDir_gauss_eb = se; 

[~,my,se] = meanyvsx_slidingwindow(reo_prevdir, nhs, txrad, deg2rad(sno.reoBinSize), 'step', true);
ad.meanNumHS_step = my;
ad.meanNumHS_step_eb = se; 
[~,my,se] = meanyvsx_slidingwindow(reo_prevdir, nhs, txrad, deg2rad(sno.reoBinSize), 'gaussian', true);
ad.meanNumHS_gauss = my;
ad.meanNumHS_gauss_eb = se; 

%create similar distributions to reorientations but within runs instead
runStartTheta = sc.runStartTheta(sc.runvalid);
runEndTheta = sc.runEndTheta(sc.runvalid);
topref = cos(runStartTheta - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(runStartTheta - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(runStartTheta - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(runStartTheta - deg2rad(sno.preferredDirection)) < -1/sqrt(2);

indset = {topref, frompref, leftofpref, rightofpref};

dt = diff(unwrap([runStartTheta;runEndTheta]));
dt = adjustForPolarHistogram(dt, deg2rad(ad.dtx)); 
for j = 1:length(indset)
    ad.run_dtheta_dist{j} = hist(dt(indset{j}), deg2rad(ad.dtx));
    ad.run_dtheta_dist_eb{j} = sqrt(ad.run_dtheta_dist{j}.*(1-ad.run_dtheta_dist{j}/nnz(indset{j})));    
end


%on to head sweeps

ad.hstx = (-180:sno.hsBinSpacing:(180-sno.hsBinSpacing)) + sno.preferredDirection;
ad.hstxc = (-180:sno.hsBinSpacing:(180)) + sno.preferredDirection;
hstx = ad.hstx;
%hstxc = ad.hstxc;
hstxe = ad.hstxc - sno.hsBinSpacing/2;
hbw = sno.hsBinSize/2;

hsinds = sc.hsvalid & (sc.hs_htv | sno.useInvalidHeadSweeps) & abs(sc.hs_maxtheta) > deg2rad(sno.minHSTheta);
hstail = sc.hs_taildir(hsinds);
hshead = sc.hs_headdir(hsinds);
hsoffset = sc.hs_dirOffset(hsinds);
hsacc = sc.hs_acc(hsinds);
hsmt = sc.hs_maxtheta(hsinds);
hssign = sc.hs_sign(hsinds);

hstail = diff(unwrap([hsoffset;hstail]));
hshead = diff(unwrap([hsoffset;hshead]));

toleft = logical(hssign > 0);
toright = logical(hssign < 0);
tdeg = rad2deg(hstail);
tdeg = mod(tdeg, 360);
tdeg (tdeg > max(hstxe)) = tdeg(tdeg > max(hstxe)) - 360;
tdeg (tdeg < min(hstxe)) = tdeg(tdeg < min(hstxe)) + 360;




%[~,bin] = histc(tdeg, hstxe);
bin = -ones(size(tdeg));
for k = 1:length(hstx)
    bin(tdeg >= hstx(k)-hbw & tdeg <= hstx(k)+hbw) = k;
end
towards = cos(deg2rad(sno.preferredDirection) - hshead) > cos(deg2rad(sno.preferredDirection) - hstail);
[~,my,se] = meanyvsx(adjustForPolarHistogram(hstail,deg2rad(hstx)), double(toleft), deg2rad(hstxe));
    
ad.allHSDir = my;
ad.allHSDir_eb = se;
for k = 1:length(hstx)
 
    ad.headSwingAcceptanceRateRight(k) = mean(hsacc(bin == k & toright));
    ad.headSwingAcceptanceRateRight_eb(k) = sqrt(ad.headSwingAcceptanceRateRight(k) * (1-ad.headSwingAcceptanceRateRight(k))) / sqrt(nnz(bin ==k & toright));

    ad.headSwingAcceptanceRateLeft(k) = mean(hsacc(bin == k & toleft));
    ad.headSwingAcceptanceRateLeft_eb(k) = sqrt(ad.headSwingAcceptanceRateLeft(k) * (1-ad.headSwingAcceptanceRateLeft(k))) / sqrt(nnz(bin ==k & toleft));

    ad.headSwingRejectionBiasRight(k) = nnz(~hsacc(bin == k & toright))/nnz(~hsacc(bin == k)) / mean(toright(bin == k));
    ad.headSwingRejectionBiasRight_eb(k) = sqrt(nnz(~hsacc(bin == k & toright)))/nnz(~hsacc(bin == k)) / mean(toright(bin == k)); %this is a kludge for now
    
    ad.headSwingRejectionBiasLeft(k) = nnz(~hsacc(bin == k & toleft))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k));
    ad.headSwingRejectionBiasLeft_eb(k) = sqrt(nnz(~hsacc(bin == k & toleft)))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k)); %this is a kludge for now
    
    ad.meanRejectedHeadSwingDir(k) = mean(hsmt(bin == k & ~hsacc));
    ad.meanRejectedHeadSwingDir_eb(k) = std(hsmt(bin == k & ~hsacc))/sqrt(nnz(bin == k & ~hsacc));

    ad.meanAcceptedHeadSwingDir(k) = mean(hsmt(bin == k & hsacc));
    ad.meanAcceptedHeadSwingDir_eb(k) = std(hsmt(bin == k & hsacc))/sqrt(nnz(bin == k & hsacc));
end

perp = abs(sind(sno.preferredDirection - tdeg)) > cosd(hbw);

ad.headSwingAcceptanceRateTowards = mean(hsacc(perp & towards));
ad.headSwingAcceptanceRateTowards_eb = sqrt(ad.headSwingAcceptanceRateTowards * (1-ad.headSwingAcceptanceRateTowards)) / sqrt(nnz(perp & towards));

ad.headSwingAcceptanceRateAway = mean(hsacc(perp & ~towards));
ad.headSwingAcceptanceRateAway_eb = sqrt(ad.headSwingAcceptanceRateAway * (1-ad.headSwingAcceptanceRateAway)) / sqrt(nnz(perp & ~towards));

ad.headSwingAcceptanceRateIndex = (ad.headSwingAcceptanceRateTowards - ad.headSwingAcceptanceRateAway)/(ad.headSwingAcceptanceRateTowards + ad.headSwingAcceptanceRateAway);
ad.headSwingAcceptanceRateIndex_eb = sqrt(var(hsacc(perp & towards)) * nnz(perp & towards) + var(hsacc(perp & ~towards)) * nnz(perp & ~towards))/nnz(perp)/(ad.headSwingAcceptanceRateTowards + ad.headSwingAcceptanceRateAway);




%repeat head sweep analysis, but only for first head sweeps

hsinds = sc.hsvalid & (sc.hs_htv | sno.useInvalidHeadSweeps) & abs(sc.hs_maxtheta) > deg2rad(sno.minHSTheta) & sc.hsnum == 1;
hstail = sc.hs_taildir(hsinds);

hshead = sc.hs_headdir(hsinds);
hsoffset = sc.hs_dirOffset(hsinds);
hsacc = sc.hs_acc(hsinds);
hsmt = sc.hs_maxtheta(hsinds);
hssign = sc.hs_sign(hsinds);

hstail = diff(unwrap([hsoffset;hstail]));
hshead = diff(unwrap([hsoffset;hshead]));

toleft = logical(hssign > 0);
%toright = logical(hssign < 0);
tdeg = rad2deg(hstail);
tdeg = mod(tdeg, 360);
tdeg (tdeg > max(hstxe)) = tdeg(tdeg > max(hstxe)) - 360;
tdeg (tdeg < min(hstxe)) = tdeg(tdeg < min(hstxe)) + 360;



%[~,bin] = histc(tdeg, hstxe);
bin = -ones(size(tdeg));
for k = 1:length(hstx)
    bin(tdeg >= hstx(k)-hbw & tdeg <= hstx(k)+hbw) = k;
end
    
towards = cos(deg2rad(sno.preferredDirection) - hshead) > cos(deg2rad(sno.preferredDirection) - hstail);

td = adjustForPolarHistogram(hstail, deg2rad(hstx));
[~,my,se] = meanyvsx(td, towards, deg2rad(hstxe));
ad.firstHSBias = my;
ad.firstHSBias_eb =se;
[~,my,se] = meanyvsx(td, toleft, deg2rad(hstxe));
ad.firstHSDir = my;
ad.firstHSDir_eb = se;

[~,my,se] = meanyvsx(td, hsmt, deg2rad(hstxe));
ad.firstHSMeanDir = my;
ad.firstHSMeanDir_eb = se;


perp = abs(sind(sno.preferredDirection - tdeg)) > cosd(hbw);

ad.firstHSProbTowards = mean(towards(perp));
x = mean(towards(perp));
n = nnz(perp);
ad.firstHSProbTowards_eb = sqrt(x*(1-x)/n);


fn = {'meanrunchange','instantaneousdthetavstheta', 'instantaneousdthetavstheta', 'speedVsDir'};
tx = ad.tx;
for j = 1:length(fn)
    yd = ad.([fn{j} '_step']);
    eb = ad.([fn{j} '_step_eb']);
    inds = zeros(size(tx));
    for k = 1:length(tx)
        [~,I] = min(abs(ad.txf - tx(k)));
        inds(k) = I;
    end
    ad.(fn{j}) = yd(inds);
    ad.([fn{j} '_eb']) = eb(inds);
end
fn = {'reoMag','reoDir','meanNumHS'};
tx = ad.reotx;
for j = 1:length(fn)
    yd = ad.([fn{j} '_step']);
    eb = ad.([fn{j} '_step_eb']);
    inds = zeros(size(tx));
    for k = 1:length(tx)
        [~,I] = min(abs(ad.txf - tx(k)));
        inds(k) = I;
    end
    ad.(fn{j}) = yd(inds);
    ad.([fn{j} '_eb']) = eb(inds);
end
if (isfield(sc, 'eti'))
    ad.time_axis = (min(sc.eti(sc.valid))+sno.timeBinSizeForTemporalPlots/2):sno.timeBinSizeForTemporalPlots:(max(sc.eti(sc.valid))-sno.timeBinSizeForTemporalPlots/2);
    txedg = [ad.time_axis-sno.timeBinSizeForTemporalPlots/2 ad.time_axis(end)+sno.timeBinSizeForTemporalPlots/2];
    [~,ad.speed_vs_time,se] = meanyvsx(sc.eti(sc.valid), sc.speed(sc.valid), txedg);
    ad.speed_vs_time_eb = kappa*se;
    
    [~,ad.speed_inruns_vs_time,se] = meanyvsx(sc.eti(sc.valid&sc.isrun), sc.speed(sc.valid&sc.isrun), txedg);
    ad.speed_inruns_vs_time_eb = kappa*se;
    
    h1 = histc(sc.reo_eti(sc.reovalid & sc.numhs >= 1), txedg);
    h2 = histc(sc.eti(sc.valid), txedg);
    
    h2 = h2(1:end-1) * sc.interptime / 60;
    h1 = h1(1:end-1);
    ad.reo_vs_time = h1./h2;
    ad.reo_vs_time_eb = sqrt(h1)./h2;
    if (isfield(sc, 'temperature'))
        [~,ad.temperature_vs_time] = meanyvsx(sc.eti(sc.valid), sc.temperature(sc.valid), txedg);
    end
end


