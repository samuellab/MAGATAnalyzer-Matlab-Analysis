function [num2clean,reportstring] = getReport(ecl, eset)
% function getReport(ecl, eset)
%
% displays a report enumerating how many tracks and points will be removed,
% and for which reasons
%
% outputs: 
%   NUM2CLEAN - number of tracks that will be removed
% inputs:
%   ecl < ESetCleaner
%   eset < ExperimentSet
reportstring = {};
sp = eset.gatherField('speed', 'mean');

npts = eset.gatherField('npts');

nt = length(sp);
spx = (0:(2/nt):1) * (max(sp) - min(sp)) + min(sp);
msg = [num2str(sum(sp < ecl.minSpeed)) '/' num2str(length(sp)) ' tracks fail speed test'];

if (ecl.showFigsInReport)
    figure(); hist(sp, spx); hold on;
    yl = get(gca, 'YLim');
    plot (ecl.minSpeed * [1 1], yl, 'r--', 'LineWidth', 3);
    title(msg);
end
disp(msg); reportstring = [reportstring msg];
msg = [num2str(sum(npts(sp < ecl.minSpeed))) '/' num2str(sum(npts)) ' pts removed by speed test (' num2str(sum(npts(sp < ecl.minSpeed))/sum(npts)*100,2), ' %)'];
disp(msg); reportstring = [reportstring msg];

if (eset.expt(1).track(1).validDQName('ihtValid'))
    htv = eset.gatherField('ihtValid', 'mean');
    htvx = (0:(2/nt):1) * (max(htv) - min(htv)) + min(htv);
    msg = [num2str(sum(htv < ecl.minHTValid)) '/' num2str(length(htv)) ' tracks fail head tail valid test'];

    if (ecl.showFigsInReport)
        figure(); hist(htv, htvx); hold on;
        yl = get(gca, 'YLim');
        plot (ecl.minHTValid * [1 1], yl, 'r--', 'LineWidth', 3);
        title(msg);
    end
    disp(msg); reportstring = [reportstring msg];
    msg = [num2str(sum(npts(htv < ecl.minHTValid))) '/' num2str(sum(npts)) ' pts removed by head-tail valid test (' num2str(sum(npts(htv < ecl.minHTValid))/sum(npts)*100,2), ' %)'];
    disp(msg); reportstring = [reportstring msg];
else
    htv = ones(size(sp));
end


dst = eset.evaluateTrackExpression('max(sqrt(sum(track(1).getDerivedQuantity(''displacement'').^2)))');
dstx = (0:(2/nt):1) * (max(dst) - min(dst)) + min(dst);
msg = [num2str(sum(dst < ecl.minDist)) '/' num2str(length(dst)) ' tracks fail displacement test'];

if (ecl.showFigsInReport)
    figure(); hist(dst, dstx); hold on;
    yl = get(gca, 'YLim');
    plot (ecl.minDist * [1 1], yl, 'r--', 'LineWidth', 3);
    title(msg);
end
disp(msg); reportstring = [reportstring msg];
msg = [num2str(sum(npts(dst < ecl.minDist))) '/' num2str(sum(npts)) ' pts removed by displacement test (' num2str(sum(npts(dst < ecl.minDist))/sum(npts)*100,2), ' %)'];
disp(msg); reportstring = [reportstring msg];
nptsx = (0:(2/nt):1) * (max(npts) - min(npts)) + min(npts);
msg = [num2str(sum(npts < ecl.minPts)) '/' num2str(length(npts)) ' tracks fail npts test'];

if (ecl.showFigsInReport)
    figure(); hist(npts, nptsx); hold on;
    yl = get(gca, 'YLim');
    plot (ecl.minPts * [1 1], yl, 'r--', 'LineWidth', 3);
    title(msg);
end
disp(msg); reportstring = [reportstring msg];
msg = [num2str(sum(npts(npts < ecl.minPts))) '/' num2str(sum(npts)) ' pts removed by npts test (' num2str(sum(npts(npts < ecl.minPts))/sum(npts)*100,2), ' %)'];
disp(msg); reportstring = [reportstring msg];

pa = eset.evaluateTrackExpression('max(abs(unwrap(track.getDerivedQuantity(''theta''))))');
et = eset.evaluateTrackExpression('max(track.getDerivedQuantity(''eti'')) - min(track.getDerivedQuantity(''eti''))');

pa = pa / (2*pi);
et = et / 60;
rpm = pa./et;

msg = [num2str(sum(rpm > ecl.rpmCut & pa > ecl.minRevCut)) '/' num2str(length(npts)) ' tracks are too curvy (fail rpm cut)'];
disp(msg);
if (ecl.showFigsInReport)
    rpmx = 0:0.1:max(rpm);
    figure(); hist(rpm, rpmx); hold on;
    yl = get(gca, 'YLim');
    plot (ecl.rpmCut * [1 1], yl, 'r--', 'LineWidth', 3);
    title(msg);
end
msg = [num2str(sum(npts(rpm > ecl.rpmCut & pa > ecl.minRevCut))) '/' num2str(sum(npts)) ' pts removed by rpm test (' num2str(sum(npts(rpm > ecl.rpmCut & pa > ecl.minRevCut))/sum(npts)*100,2), ' %)'];
disp(msg); reportstring = [reportstring msg];


msg = [num2str(sum(npts < ecl.minPts)) '/' num2str(length(npts)) ' tracks fail npts test, ' num2str(sum(sp < ecl.minSpeed)) ' fail speed test, ', num2str(sum(sp < ecl.minSpeed & npts < ecl.minPts)) ' fail both.'];
if (ecl.showFigsInReport)
    I = htv >= ecl.minHTValid;
    figure(); plot (npts(I), sp(I), 'b.', npts(~I), sp(~I), 'r.'); xlabel ('npts'); ylabel('speed'); hold on;
    yl = get(gca, 'YLim');
    plot (ecl.minPts * [1 1], yl, 'r--', 'LineWidth', 3);
    xl = get(gca, 'XLim');
    plot (xl, ecl.minSpeed * [1 1], 'r--', 'LineWidth', 3);
    legend ({'htv passes', 'htv fails'});
    title(msg);
end

msg = [num2str(sum(htv < ecl.minHTValid)) '/' num2str(length(htv)) ' tracks fail head-tail valid test, ' num2str(sum(sp < ecl.minSpeed)) ' fail speed test, ', num2str(sum(sp < ecl.minSpeed & htv < ecl.minHTValid)) ' fail both.'];
if (ecl.showFigsInReport)
    I = npts >= ecl.minPts;
    figure(); plot (htv(I), sp(I), 'b.',htv(~I), sp(~I), 'r.'); xlabel ('htvalid'); ylabel('speed'); hold on;
    yl = get(gca, 'YLim');
    plot (ecl.minHTValid * [1 1], yl, 'r--', 'LineWidth', 3);
    xl = get(gca, 'XLim');
    plot (xl, ecl.minSpeed * [1 1], 'r--', 'LineWidth', 3);
    legend ({'npts passes', 'npts fails'});
    title(msg);
end

msg = [num2str(sum(htv < ecl.minHTValid)) '/' num2str(length(htv)) ' tracks fail head-tail valid test, ' num2str(sum(npts < ecl.minPts)) ' tracks fail npts test, ', num2str(sum(npts < ecl.minPts & htv < ecl.minHTValid)) ' fail both.'];
if (ecl.showFigsInReport)
    figure(); plot (htv(sp > ecl.minSpeed), npts(sp > ecl.minSpeed), 'b.', htv(sp < ecl.minSpeed), npts(sp < ecl.minSpeed),'r.'); xlabel ('htvalid'); ylabel('npts'); hold on;
    yl = get(gca, 'YLim');
    plot (ecl.minHTValid * [1 1], yl, 'r--', 'LineWidth', 3);
    xl = get(gca, 'XLim');
    plot (xl, ecl.minPts * [1 1], 'r--', 'LineWidth', 3);
    title(msg);
    legend ({'speed passes', 'speed fails'});
end

nv = htv < ecl.minHTValid | npts < ecl.minPts | dst < ecl.minDist | sp < ecl.minSpeed | (rpm > ecl.rpmCut & pa > ecl.minRevCut);
num2clean = sum(nv);
msg = [num2str(num2clean) '/' num2str(length(htv)) ' tracks fail at least one test'];
disp(msg); reportstring = [reportstring msg];
msg = [num2str(sum(npts(nv))) '/' num2str(sum(npts)) ' pts removed by eset cleaner (' num2str(sum(npts(nv))/sum(npts)*100,2), ' %)'];
disp(msg); reportstring = [reportstring msg];
