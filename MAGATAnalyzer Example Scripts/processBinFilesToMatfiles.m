function [eset, success] = processBinFilesToMatfiles (basedir, esetname, varargin)
%function eset = loadTrimStitchAndSave(basedir, esetname, ecl, camcalinfo, varargin)
%loads, stitches, trims, cleans, etc. then saves to mat files in a
%'matfiles' subdirectory 
%varargin: parameter = defaultvalue
%          to change, pass 'parameter', value after basedir, esetname
%          checkerboardfname = 'checkerboard.png';
%          checkeroptions = {'xinds', 310:2170, 'yinds', 220:1640, 'flipy', true, 'flipx', false, 'flatten', true};
%          minptsToLoad = 50;
%          frameDiff = 4;
%          maxDist = 0.1;
%          buffer = [-0.1 0.2 0.1 0.1];
%          trimrect = [];
%          fieldsToDerive = {};
%
% explanation of parameters 
%
% checkboardfname - specifies the name of the picture of the 1cm checkerboard that we
%                   will use to calibrate the camera; place the
%                   checkerboard inside basedir
% checkeroptions - options to pass to the checkerboard calibration routine
% 
% minPtsToLoad  - the minimum length a track must be before we will even
%                 load it from disk
%
% frameDiff - stitch together tracks if first ended frameDiff or fewer frames before second started
% maxDist - stitch together tracks if first ended within maxDist (in cm) of second's start
% buffer - distance from leftmost, rightmost, bottommost, topmost point to trim 
%          (to eliminate border region)
% trimrect - you can specify an actual trim rectangle here
%   if both buffer & trimrect are specified, trimrect is used & buffer is
%   ignored; if buffer = [] and trimrect = [], trimming is disabled
% fieldsToDerive - extra data fields to derive before saving to save time
%                  later


%% this section sets the parameters that will be used throughout the script
%  you can change any of these parameters by passing 
%'parameter_name',  parameter value after basedir, esetname in the
% arguments to the function

%this specifies the name of the picture of the 1cm checkerboard that we
%will use to calibrate the camera -- using the checkerboard lets us analyze
%the data - maggot lengths, run speeds, run lengths, etc. in real units
%(cm) that we can compare between experiments and labs
%
%it also takes out any lens distortion 
%
%place the checkerboard inside basedir 
%optionally you can rewrite this script to use a checkerboard you are
%storing somewhere else, but it's probably easiest just to store it with
%the extracted data anyway

checkerboardfname = '1cm_checkerboard.bmp';


%these are options to pass to the checkerboard calibration routine.  xinds,
%yinds specify the region of the image that contains the checkerboard (if
%we are zoomed out enough to see the beyond the edge of the checkerboard)
%'flipy' is true (for cold_arrows) because when we put 0,0 in the lower
%left corner, the image is upside down
%'flipx' is false (for cold_arrows) because when we put 0,0 in the lower
%left corner, warm is to +x
%'flatten', true evens out the illumination over the checkerboard, and is
%generally worth doing
checkeroptions = {'flipy', true, 'flipx', false, 'flatten', true};


% minPtsToLoad is the minimum length a track must be before we will even
% load it from disk
minptsToLoad = 50;


%these are options to use when stitching together tracks -- we do this in
%case the track extraction software dropped a frame or two in the middle of
%an otherwise good track. . . note that if the track is in
%collision with another track in the experiment, we don't stitch 

frameDiff = 4; % stitch together tracks if first ended 4 or fewer frames before second started
maxDist = 0.1; % stitch together tracks if first ended within 1 mm of second's start

buffer = [-0.1 0.2 0.1 0.1];  % [left, right, top, bottom]
                              % don't trim left (no tracks leave the left
                              % edge).  trim right by 2mm, top bottom by 1
                              % mm
trimrect = []; % you can specify an actual trim rectangle here, which will mean buffer is ignored
               % note that the window should be specified in cm, not pixels


%we derive all the fields necessary to segment the tracks and generate navigation figures before saving
%if there are other fields you want derived and saved, place them here
fieldsToDerive = {};
varargin = assignApplicable(varargin);

%we clean the experiment set twice -- once before stitching and once after
%stitching;  the cleaner (eclnukespots) before stitching just looks for spots that are
%clearly not maggot tracks -- head tail is largely invalid; does not move
%more than 2 mm
%
%after stitching we are more selective, insisting on a minimum number of
%points, a fairly high number of valid head tail determinations, a minimum
%distance travelled, and a requirement that the larva not continuously
%circle in one direction

eclnukespots = ESetCleaner();
eclnukespots.minHTValid = 0.6;
eclnukespots.minDist = 0.2;
eclnukespots.askFirst = false;


ecl = ESetCleaner();
ecl.minDist = 1; % track must be more than 1 cm long
mintime = 60; % tracks must be at least 60 seconds -- we will convert this into a minimum number of points later
ecl.minHTValid = 0.9; % discard any track where we couldn't determine a head tail properly at least 90% of the time

ecl.rpmCut = 2; % get rid of any track that circles in one direction with a frequency of at least 2 rev/min
ecl.minRevCut = 3; % provided that track does at least 3 complete circles

varargin = assignApplicable(varargin);
success = false;
ts1 = tic();

%% this section carries out the camera calibration if needed and saves it
% to disk

if (~exist(fullfile(basedir, 'camcalinfo.mat'), 'file'))
    try 
        im = imread(fullfile(basedir, checkerboardfname));
    catch me
        disp ('could not read in checkerboard from: ');
        disp (fullfile(basedir, checkerboardfname));
        disp ('error was: ')
        disp (me.getReport());
        disp ('');
        error ('cannot proceed without camera calibration -- aborting');
        return; %#ok<UNRCH>
    end
    
    disp ('calibrating checkerboard');
    
    camcalinfo = CameraCalibration(im, checkeroptions{:});
    if (isempty(camcalinfo))
        error ('camera calibration failed - check options');
    end
    % make sure that the checkerboard that displays here is correct -- I am
    % not putting a checkoff in the script because you want to run it automatically, but this is a potential place for errors
    
    
    save(fullfile(basedir, 'camcalinfo.mat'), 'camcalinfo');
else
    load(fullfile(basedir, 'camcalinfo.mat'), 'camcalinfo');
end
        
%% this section will look through the directory for .bin files and verify
% that we have the appropriate .mdat files as well
d = dir(fullfile(basedir, '*.bin'));
for j = 1:length(d)
    fname{j} = fullfile(basedir, d(j).name); %#ok<*AGROW>
    [ps,nm] = fileparts(fname{j});
    mdatname{j} = fullfile(ps,[nm '.mdat']);
    if (~exist(mdatname{j}, 'file'))
        disp ('cannot find metadata file');
        disp (['bin file name = ' fname{j}]);
        disp (['metadata file name = ' mdatname{j}]);
        error ('must have metadata file to proceed');
        return %#ok<UNRCH>
    end
end

% now we will load each mdat file and see if there is an ROI offset we need
% to take into account for the camera calibration

rx = camcalinfo.realx;
ry = camcalinfo.realy;
for j = 1:length(mdatname)
    cx = camcalinfo.camx;
    cy = camcalinfo.camy;
    data = importdata(mdatname{j});
    xcol = find(strcmpi(data.colheaders, 'ROIX'));
    if (~isempty(xcol))
        roix = median(data.data(isfinite(data.data(:,xcol)),xcol)); % we should be able just to take the first value, but this should be more robust
    else
        roix = 0;
    end
    ycol = find(strcmpi(data.colheaders, 'ROIY'));
    if (~isempty(ycol))
        roiy = median(data.data(isfinite(data.data(:,ycol)),ycol)); % we should be able just to take the first value, but this should be more robust
    else
        roiy = 0;
    end
    
    cc(j) = CameraCalibration(rx, ry, cx-roix, cy-roiy);
end

%% this command actually loads the experiments from disk -- so short!

eset = ExperimentSet.fromFiles(fname{:}, 'minpts', minptsToLoad, 'camcalinfo', cc, 'parallel', false);

%% we are going to put the rest of the file in try/catch blocks -- that way
% if there's a failure, you don't have to reload the tracks from disk to
% see what the problem is because the function will still return a value

% get rid of obviously bad tracks (dust, crap) before stitching
try
    disp ('cleaning bad spots');
    eclnukespots.clean(eset);
catch me
    disp(me.getReport);
    return;
end

%stitch tracks together
try 
    disp ('stitching tracks');
    eset.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', false);
catch me
    disp(me.getReport);
    return;
end
    
%trim the region near the edge of the experiment
try
    if (isempty(trimrect))
        if (~isempty(buffer))
            il = eset.gatherField('iloc');
            ll = min(il,[],2);
            ur = max(il,[],2);
            trimrect = [ll(1) ll(2) ur(1) ur(2)] - buffer;
        end
    end
    if (~isempty(trimrect))
        disp ('trimming tracks');
        eset.executeExperimentFunction('trimTracks', [], trimrect);
    else
        disp ('skipping trimming');
    end
catch me
    disp(me.getReport);
    return;
end

%clean out marginal tracks that will not contribute to analysis
ecl.minPts = ceil(mintime / eset.expt(1).dr.interpTime);
ecl.askFirst = false; 
ecl.showFigsInReport = false;
[~,rpt] = ecl.getReport(eset);
try
    fid = fopen (fullfile (basedir, 'cleaning report.txt'), 'wt');
    for j = 1:length(rpt)
        fprintf(fid, '%s\n', rpt{j});
    end
catch me
    disp ('trouble writing cleaning report to disk, continuing anyway');
    disp (me.getReport);
end

disp ('cleaning eset');
try
    ecl.clean(eset);
    valid = true(size(eset.expt));
    for j = 1:length(eset.expt)
        if (isempty(eset.expt(j).track))
            disp ([eset.expt(j).fname ' - after cleaning no valid tracks']);
            valid(j) = false;
        end
    end
    eset.expt = eset.expt(valid);
catch me
    disp(me.getReport);
    return;
end


%calculate derived quantities used in segmentation; set segmentation speeds

try
    eset.executeTrackFunction('setSegmentSpeeds');
    dvfields = {'sspineTheta', 'vel_dp', 'spineLength'};
    for j = 1:length(dvfields)
        disp (['calculating ' dvfields{j}]);
        eset.gatherField(dvfields{j});
    end
    eset.gatherField('sspineTheta');
    eset.gatherField('vel_dp');
    eset.gatherField('spineLength');
    for j = 1:length(fieldsToDerive) 
        disp (['calculating ' fieldsToDerive{j}]);
        eset.gatherField(fieldsToDerive{j});
    end
catch me
    disp(me.getReport);
    return;
end
disp('done with loading, stitching and cleaning');
toc(ts1)

%finally, save the eset to disk as matfiles, so it can be loaded quickly later
try
    if (~exist(fullfile(basedir, 'matfiles'), 'dir'))
        mkdir (fullfile(basedir, 'matfiles'));
    end
    eset.toMatFiles(fullfile(basedir, 'matfiles',esetname));
    disp('saved file');
    toc(ts1)
catch me
    disp(me.getReport);
    return;
end

success = true;
