existsAndDefault('username', []);
existsAndDefault('segmodelname', 'fithmmmodel2.mat');
if (exist('sm', 'var') && isa (sm, 'SegmentationModel'))
    LabelPointGUI('sm', sm,'username',username,'basedir',fullfile(pwd, 'Segmentation Models', 'N2 Cryo', segmodelname));
    return;
end   

%% Set up cryophilic point labeling
% An example script with annotations
%%


%% LOADING FILES FROM DISK
% loading specific files by name
existsAndDefault('marcmac', false);
if (marcmac)
    basedir =  '~/Documents/lab data/n2cryo/';
else
    basedir = '\\labnas1\Share\David\Extracted\Spatial\N2\18-23GradientC15\OutputFiles\';
end
d = dir([basedir '*.bin']);
nfiles = 1;
for j = 1:nfiles
    fnames{j} = [basedir d(j).name]; %#ok<SAGROW>
end
%fnames = {[basedir '20090226_N2g15_1823_tracks.bin'], [basedir '20090226_w1a_N2g15_1823_tracks.bin']};

% load any track longer than 50 points
minpts = 50;

% this code snippet loads the files if we haven't already loaded them, but
% otherwise skips them; that way we can change the script and rerun it
% without having to reload the files
if (~exist('cryo', 'var'))
    cryo = ExperimentSet.fromFiles(fnames{:}, 'minpts', minpts);
end

%% STITCH TRACKS
% sometimes we miss a frame, so let's stitch together tracks that are close
% by

frameDiff = 3; % stitch together tracks if first ended 3 or fewer frames before second started
maxDist = 7; % stitch together tracks if first ended within 7 pixels of second's start

% For the script, I am executing this function with interactive off, but if
% you set interactive to true, it will show you each potential stitch and
% let you decide whether or not to stitch it
cryo.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', false);

%% CLEAN UP TRACKS
% get rid of any tracks that don't go anywhere

% create an EsetCleaner object

ecl = ESetCleaner();

% now let's look at the autogenerated report
% let's get rid of all tracks less than 750 points and speed less than 0.4
% pixels per second
ecl.minPts = 750;
ecl.minSpeed = 0.9;

%ecl.getReport(cryo);

% the following code just forces the figures to appear in the example documentation
%{
for j = 1:3
    figure(j);
    snapnow; 
end
%}

% we've already shown the report, so we don't need to have it ask us first,
% for the purposes of this script;  generally a good idea to leave this
% enabled
ecl.askFirst = false; 

ecl.clean(cryo);

%% Segmenting Tracks
% separate tracks into runs and reorientations

% the default segmentation options are pretty good
WormSegmentOptions

% but just for the heck of it, let's change the minimum run time

wso = WormSegmentOptions;
wso.minRunTime = 3; % seconds

% segment the tracks
cryo.executeExperimentFunction('segmentTracks', wso);

% pick a track and plot the segmentation

%% add lowpassed theta derivative
%{
if (isempty([cryo.expt.globalQuantity]) || ~any(strcmpi({cryo.expt.globalQuantity.fieldname}, 'lrdtheta')))
    gq = GlobalQuantity();
    gq.fieldname = 'lrdtheta';
    gq.xField = 'theta';
    gq.xData = 5; %smoothing time in seconds
    gq.derivationMethod = @(xin, xData, yData) deriv(unwrap(xin), xData(1));

   % [cryo.expt.globalQuantity] = deal([]);
    cryo.executeExperimentFunction('addGlobalQuantity', gq);
end
%}
%% load segmentation model
sm = SegmentationModel.fromMatFile(fullfile('Segmentation Models', 'N2 Cryo', segmodelname),cryo);
LabelPointGUI('sm', sm,'username',username);