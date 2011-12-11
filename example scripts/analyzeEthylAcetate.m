fnlist = {'D:\Marc Data\20010205\CS2\CS2_tracks_2.bin', ...
          'D:\Marc Data\20010205\CS3\CS3_tracks.bin',...
          'D:\Marc Data\20010205\CS4\CS4_tracks.bin'};
timfnlist = {'D:\Marc Data\20010205\CS2\CS2_.tim', ...
          'D:\Marc Data\20010205\CS3\CS3_.tim',...
          'D:\Marc Data\20010205\CS4\CS4_.tim'};
      
%load tracks from 3 experiments      
if ~exist ('expt', 'var')
    for j = 1:3
        %load tracks including contours (no images though)
        %don't calibrate to camera (lengths are in pixels)
        %don't load tracks under 400 points (100 seconds)
        expt(j) = Experiment.fromFile(fnlist{j}, timfnlist{j}, true, [], 400);
        
        %fix any gross head tail errors
        expt(j).executeTrackFunction('fixHTOrientation');
    

        validRect = [325 2200 125 1700]; %tracks must start inside this rectangle to be valid

        startPt = zeros([2 length(expt(j).track)]);
        for k = 1:length(expt(j).track)
            startPt(:,k) = expt(j).track(k).pt(1).loc;
        end
        valid = insideRect(validRect, startPt);
        figure(j);
        for k = find(valid)
            expt.track(k).plotPath('sloc', 'b-'); hold on
        end
        for k = find(~valid)
            expt.track(k).plotPath('sloc', 'r-');
        end
        hold off
        expt(j).track = expt(j).track(valid);
    end
end

if (isempty([expt(1).track.run]) || (exist('startover', 'var') && ~isempty(startover) && startover)) 
    so = MaggotSegmentOptions();
    for j = 1:3
        expt(j).so = so;
        expt(j).executeTrackFunction('recalculateDerivedQuantities');
        expt(j).executeTrackFunction('segmentTrack', so);
    end
    startover = false;
end
runstartdir = [];
runenddir = [];
alldir = [];
hsdir = [];
hsaccepted = [];
hsang = [];
hsprevdir = [];
hstaildir = [];
runlen = [];
rundir = [];
for j = 1:3
    runstartdir = [runstartdir rad2deg(expt(j).gatherSubField('run', 'startTheta'))];
    runenddir = [runenddir rad2deg(expt(j).gatherSubField('run', 'endTheta'))];
    alldir = [alldir rad2deg(expt(j).gatherField('theta','run'))];
    hsdir = [hsdir rad2deg(expt(j).gatherSubField('headSwing', 'headDir'))];
    hsaccepted = [hsaccepted (expt(j).gatherSubField('headSwing', 'accepted'))];
    hsang = [hsang rad2deg(expt(j).gatherSubField('headSwing', 'maxTheta'))];
    hsprevdir = [hsprevdir rad2deg(expt(j).gatherSubField('headSwing', 'prevDir'))];
    hstaildir = [hstaildir rad2deg(expt(j).gatherSubField('headSwing', 'tailDir'))];
    runlen = [runlen expt(j).gatherSubField('run', 'pathLength')];
    rundir = [rundir rad2deg(expt(j).gatherSubField('run', 'meanTheta'))];
end
close all;
fignum = 0;
fignum = fignum + 1;
figure(fignum);clf
tx = -180:30:180;
plot (tx,hist(runenddir,tx)./hist(alldir,tx) * 60 / expt(j).dr.interpTime, 'LineWidth', 3);
xlabel ('heading');
ylabel ('reorientation rate (min^{-1})');
title ('probability of ending run vs. instantaneous heading');
embiggen();

fignum = fignum + 1;
figure(fignum);clf
tx = -180:30:180;
h1 = hist(runstartdir,tx);
h1(1) = h1(1) + h1(end); h1(end) = h1(1);
h2 = hist(runenddir,tx);
h2(1) = h2(1) + h2(end); h2(end) = h2(1);
plot (tx, h1, 'g-', tx, h2, 'r-', 'LineWidth', 3); ylim([0 max([h1 h2])]);
xlabel ('heading');
ylabel ('# runs');
title ('distribution of run start and end directions');
legend ('start', 'end');
embiggen();

fignum = fignum + 1;
figure(fignum); clf
rlx = 50:50:1000;
semilogy (rlx, hist(runlen(abs(rundir - 90) < 45), rlx)/sum(abs(rundir - 90) < 45), 'r-', rlx, hist(runlen(abs(rundir + 90) < 45), rlx)/sum(abs(rundir + 90) < 45), 'g-',...
    rlx, hist (runlen(abs(rundir) < 45),rlx)/sum(abs(rundir) < 45), 'm-', rlx, hist(runlen(abs(rundir) > 135),rlx)/sum(abs(rundir) > 135), 'c-',...
    'LineWidth', 3);
legend ('against gradient', 'with gradient', 'to left of gradient', 'to right of gradient');
xlabel ('run length (pixels)');
ylabel ('fraction of runs');
%{
fignum = fignum + 1;
figure(fignum);clf
[x,meany] = meanyvsx (hstaildir, abs(hsang), tx);
plot (x, meany, 'LineWidth', 3);
xlabel ('previous direction');
ylabel ('mean headsweep magnitude (degrees)');
title ('size of head sweep vs. previous direction');
embiggen()

inds = find(hsaccepted);
inds2 = find(~hsaccepted);
fignum = fignum + 1;
figure(fignum);clf
[x1,meany1] = meanyvsx (hstaildir(inds), hsang(inds), tx);
[x2,meany2] = meanyvsx (hstaildir(inds2), hsang(inds2), tx);
plot (x1,meany1, 'g-', x2, meany2, 'r-', 'LineWidth', 3);
title ('mean headsweep angle vs. initial angle');
xlabel ('previous direction');
ylabel ('mean headsweep angle (degrees)');
legend ('accepted head sweeps', 'rejected head sweeps');
%}

rightup = find(abs(hstaildir) < 45 & hsang > 0);
rightdown = find(abs(hstaildir) < 45 & hsang < 0);
leftup = find(abs(hstaildir) > 135 & hsang < 0);
leftdown = find(abs(hstaildir) > 135 & hsang > 0);
