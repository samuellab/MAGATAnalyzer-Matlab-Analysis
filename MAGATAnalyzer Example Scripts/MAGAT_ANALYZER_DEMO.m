% starting with raw image files, we
% (1) run a c++ program to extract maggot positions and postures from the
%     video
% (2) load the maggot positions and postures into matlab, perform
%     associated calibrations, derive behavioral metrics, then save the
%     result to disk, allowing for quick loading later
% (3) segment each track into an alternating series of runs and turns
% (4) calculate metrics of strategy, make example figures, save those
%     figures to disk
% (5) as examples, plots the path of a single track, 
%     plays a movie of a track, 
%     plays a movie of every left turn in a track 
%
% the whole process, start to finish should take < ~1 hr on a desktop computer
% after this script has been run to completion, subsequent executions
% should take <5 minutes
%
% before running this script, or any other MAGAT Analyzer function, you
% must call the function setupDirectories and add MAGATAnalyzer Example Scripts and subdirectories to the path
% (running MAGAT_ANALYZER_START_HERE accomplishes this)
%
% at the conclusion of the script, the Example Data directory should have
% these subdirectories
%       Image Data - raw movies taken of larvae in LADY GAGA, saved in a
%                    special compressed format.  
%       Extracted Tracks\CS EtAc 2ppm per cm spatial gradient - the results
%           of the track extraction software are stored here
%       Extracted Tracks\CS EtAc 2ppm per cm spatial gradient subdirectories:
%           matfiles - extracted tracks processed by matlab for fast
%               loading and computations in the future
%           calculations - spatial analysis results used to make figures
%           figures - pdfs with figures similar to those we use for
%               screening (simpleMetrics) and navigational analysis
%               (strategyMetrics)
datadirectory = fullfile ('..', 'Example Data');
extracteddirectory = fullfile(datadirectory, 'Extracted Tracks');
dataname = 'CS EtAc 2ppm per cm spatial gradient';
basedir = fullfile(extracteddirectory, dataname);

% step 1
d = dir(fullfile(basedir, '*.bin'));
if (isempty(d))
    disp ('extracting tracks from raw image files');
    ExtractTracksFromRawData;
end

esetname = 'CS-EtAc-2ppm-per-cm';

% step 2
d = dir(fullfile(basedir, 'matfiles', [esetname '_experiment*.mat']));
if (isempty(d))
    disp ('loading tracks from bin files');
    [eset,success] = processBinFilesToMatfiles(basedir, esetname);
    if (~success)
        disp('eset processing bins to matfiles failed');
        return
    end
else
    if (~exist('eset', 'var'))
        disp ('loading tracks from mat files');
        eset = ExperimentSet.fromMatFiles(fullfile(basedir, 'matfiles', esetname));
    end
end

% step 3
disp('segmenting tracks into runs and turns');
eset.executeTrackFunction('segmentTrack');

% step 4
disp ('calculating navigational metrics');
spatialNavigationOptions = spatialMaggotAnalysis;
spatialNavigationOptions.validname = 'eti';
spatialNavigationOptions.validoperation = '@(x) x > 120 & x < 1020'; %only consider time between 2 & 17 minutes

spatialCalculations = spatialMaggotCalculations(eset, spatialNavigationOptions);
analyzedData = spatialMaggotAnalysis(spatialCalculations, spatialNavigationOptions);
if (~exist(fullfile(basedir, 'calculations'), 'dir'))
    mkdir(fullfile(basedir, 'calculations'));
end
disp ('saving calculations to disk');
save(fullfile(basedir, 'calculations', 'spatial_analysis.mat'), 'spatialCalculations', 'analyzedData');

simpleMetricFigure(analyzedData, 'CS Navigating 2ppm/cm EtAc', [],[],'fignum', 1);
strategyMetricFigure(analyzedData,[], 'figureTitle', 'CS Navigating 2ppm/cm EtAc', 'fignum', 2);

% save the figures to disk as pdfs
if (~exist(fullfile(basedir, 'figures'), 'dir'))
    mkdir(fullfile(basedir, 'figures'));
end
saveas(1, fullfile(basedir, 'figures', 'simpleMetrics.pdf'), 'pdf');
saveas(2, fullfile(basedir, 'figures', 'strategyMetrics.pdf'), 'pdf');

% step 5

%find a pretty & long track
htv = eset.gatherField('ihtValid', 'mean');
npts = eset.gatherField('npts');
[~,I] = sort(htv, 'descend');
tind = I(find(npts(I) > 4000,1, 'first')); 
track = eset.indToTrack(tind);

% plot the track path
figure(10); clf
track.plotPath('sloc', 'b-', 'highlightinds', 1, 'highlightlinetype', 'g.','MarkerSize',20); axis equal;
%track will be plotted as a blue line with a green dot at the starting
%position;  axes are labeled in cm



%play a short excerpt of the track 
figure(11); clf
track.playMovie('startTime', 485, 'stopTime', 605);

%now play every turn(reorientation) with at least one headsweep where the larva
%turns to the left
track.reorientation([track.reorientation.numHS] > 0 & diff(unwrap([[track.reorientation.prevDir];[track.reorientation.nextDir]])) > 0).playMovie('nopause', true);