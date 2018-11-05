function eset = loadTrimStitchAndSave(basedir, esetname, ecl, camcalinfo, varargin)
%function eset = loadTrimStitchAndSave(basedir, esetname, ecl, camcalinfo, varargin)
%loads, stitches, trims, cleans, etc. then saves to mat files in a
%'matfiles' subdirectory 
%varargin: minpts - minimum length of a track in frames to load from disk
%          frameDiff - when stitching maximum number of missing frames between tracks
%          maxDist - when stitching maximum physical distance between
%                    tracks
%          trimrect - rectangle outside of which to remove points 
%          buffer - if trimrect is unspecified, distance from edge of analysis rectangle to trim
%          prunebuffer - if prunebuffer is specified, distance from edge of
%                        trimrectangle to prune tracks (remove all that start outside this edge) 
%                        NOTE: this distance is applied to the trimrectangle 
%          default_deltaT - if loading timing fails, what is the default
%                           time interval between frames?
%          checker_calc - if this was a checkerboard experiment, you can add checker data here
%          timerange - remove points outside this time window (usually this is left empty)
%          fieldsToDerive - additional fields to derive before saving
%                       (e.g. 'periAmp')
%          setDRByPeristalsis - whether to automatically set derivation
%               rules by using automatically detected peristalsis frequency
%           isMWT, true/[false] -- if this is an MWT experiment



%Set params
minpts = 50;
frameDiff = 2; % stitch together tracks if first ended 2 or fewer frames before second started
if (isempty(camcalinfo))
    maxDist = 10; % stitch together tracks if first ended within 10 pixels of second's start
else
    maxDist = 0.1; %one mm
end
prunebuffer = [];
buffer = [];
trimrect = [];
default_deltaT = 0.2;

checker_calc = [];
timerange = [];
fieldsToDerive = {};
setDRByPeristalsis = false;
isMWT = false;
varargin = assignApplicable(varargin);
eclexisted = existsAndDefault('ecl', ESetCleaner());
existsAndDefault('camcalinfo', []);

eclnukespots = ESetCleaner();
eclnukespots.minHTValid = 0.5;
if (~isempty(camcalinfo)) 
    eclnukespots.minDist = 0.1;
    eclnukespots.minHTValid = 0.5;
    
end
eclnukespots.askFirst = false;
ts1 = tic;

%LOAD the experiments
if (isMWT)
    eset = ExperimentSet.fromMWTFiles(basedir, camcalinfo, varargin{:});
else
    eset = ExperimentSet.fromFiles(basedir, 'minpts', minpts, 'camcalinfo', camcalinfo, 'parallel', true, varargin{:});
end

%Optionally, generate derivation rules using properties of the track
if (setDRByPeristalsis)
    try
        eset.setDerivationRulesByPeristalsisFrequency;
    catch me
        disp(me.getReport());
        disp ('failed to set derivation rules; aborting now and returning so you still have eset');
        return;
    end     
end

%Add timing info (if it's not already there)
% note: in ExperimentSet.fromFiles, timing info is NOT there by default 
if (~isMWT)
    eset.addTimingByFrameDifference(default_deltaT);
end

%CLEAN the funky-looking tracks out of the eset
eclnukespots.clean(eset);

%Set some ecl field values
if (~eclexisted ) 
    if (~isempty(camcalinfo))
    %real units, assume cm
        ecl.minDist = 0.1; %minimum distance 1 mm
        ecl.minSpeed = 0.01; %minimum average speed 100 microns/sec
    else
        ecl.minDist = 10; %pixels
        ecl.minSpeed = 0.1; %pixels/second
    end
    ecl.minPts = ceil(30 / eset.expt(1).dr.interpTime);
end

%STITCH appropriate tracks together
eset.executeExperimentFunction('stitchTracks', frameDiff, maxDist, 'interactive', false);

%Set some more ecl field values
ecl.askFirst = false; 
ecl.showFigsInReport = false;
ecl.getReport(eset);

%CLEAN the weird-looking tracks out of the experients again (after stitching tracks together)
ecl.clean(eset);

%Clear the experiments which have 0 tracks from the eset
valid = true(size(eset.expt));
for j = 1:length(eset.expt)
    if (isempty(eset.expt(j).track))
        disp ([eset.expt(j).fname ' - after cleaning no valid tracks']);
        valid(j) = false;
    end
end
eset.expt = eset.expt(valid);

%TRIM any points from the track which fall out of the valid time/location
%ranges
if (isempty(trimrect)) 
    il = eset.gatherField('iloc');
    if (isempty(buffer))
        if(isempty(camcalinfo))
            buffer = 25; %pixels
        else
            buffer = 0.15; %cm
        end
    end
    ll = min(il,[],2) + buffer;
    ur = max(il,[],2) - buffer;
    trimrect = [ll(1) ll(2) ur(1) ur(2)];
end
disp ('trimming tracks');
eset.executeExperimentFunction('trimTracks', timerange, trimrect);
if (~isempty(prunebuffer))
    prunerect = trimrect + [prunebuffer prunebuffer -prunebuffer -prunebuffer];
    disp ('pruning tracks');
    eset.executeExperimentFunction('trimTracks', [], prunerect);
end

%CALCULATE derived quantities
disp('calculating derived quantities');
if (isa (eset.expt(1).track(1), 'MaggotTrack'))
    eset.executeTrackFunction('setSegmentSpeeds');
    eset.gatherField('sspineTheta');
    eset.gatherField('vel_dp');
end
for j = 1:length(fieldsToDerive) 
    try
        eset.gatherField(fieldsToDerive{j});
    catch me
        disp(me.getReport());
    end
end


disp('done with loading, stitching and cleaning');
toc(ts1)

%Assign checker track data
if (~isempty(checker_calc))
    try
        disp ('assigning checker track data');
        checker_calc.assignGlobals(eset.expt);
    catch me
        disp(me.getReport());
    end
    toc(ts1);
end

%SAVE the eset as a matfile
ts1 = tic;
if (~exist(fullfile(basedir, 'matfiles'), 'dir'))
    mkdir (fullfile(basedir, 'matfiles'));
end
eset.toMatFiles(fullfile(basedir, 'matfiles',esetname));
disp('saved file');
toc(ts1)

