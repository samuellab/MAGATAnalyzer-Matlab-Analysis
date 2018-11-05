function [eset, success] = processJavaDirectoryToMatfiles (basedir,varargin)
%function eset = loadTrimStitchAndSave(basedir, esetname, ecl, camcalinfo, varargin)
%loads, stitches, trims, cleans, etc. then saves to mat files in a
%'matfiles' subdirectory 
%varargin: parameter = defaultvalue
%          to change, pass 'parameter', value after basedir, esetname
%          checkerLocation = '';
%          checkerboardfname = 'checkerboard.png';
%          checkeroptions = {'xinds', 310:2170, 'yinds', 220:1640, 'flipy', true, 'flipx', false, 'flatten', true};
%          minptsToLoad = 50;
%          frameDiff = 4;
%          maxDist = 0.1;
%          buffer = [-0.1 0.2 0.1 0.1];
%          trimrect = [];
%          fieldsToDerive = {};
%
% explanation of parameters (more detail found inside
%                            processDirectoryToMatfiles_Janelia.m)
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
% maxDist - stitch together tracks if first ended within maxDist (in cm) of
% second's start
% buffer - distance from leftmost, rightmost, bottommost, topmost point to trim 
%          (to eliminate border region)
% trimrect - you can specify an actual trim rectangle here
%   if both buffer & trimrect are specified, trimrect is used & buffer is
%   ignored; if buffer = [] and trimrect = [], trimming is disabled
% fieldsToDerive - extra data fields to derive before saving to save time
%                  later
% ccInSupDataDir - true/[false] -- if true, cameracalibration info is
%   stored in supplemental datadir
%
% individualMatFiles - whether to process the experiment mat files
%                      pseudo-individually using the same name as the .jav files 
%
% requireMaggot = [true]/false -- if true, exits without loading if
% experiment does not contain maggot tracks
%




% this section sets the parameters that will be used throughout the script
%  you can change any of these parameters by passing 
%'parameter_name',  parameter value after basedir, esetname in the
% arguments to the function

%this specifies the name of the picture of the 1cm checkerboard that we
%will use to calibrate the camera -- using the checkerboard lets us analyze
%the data - maggot lengths, run speeds, run lengths, etc. in real units
%(cm) that we can compare between experiments and labs
%
%place the checkerboard inside basedir 
%optionally you can provide a different location for the checkerboard

checkerboardfname = '1cm_checkerboard.png';

checkerLocation = '';


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

% loadParallel is whether to use parallel processes to load mat files.  set
% to false if you are running on a computer with "limited" (say < 8 GB)
% memory
loadParallel = false;

%these are options to use when stitching together tracks -- we do this in
%case the track extraction software dropped a frame or two in the middle of
%an otherwise good track. . . note that if the track is in
%collision with another track in the experiment, we don't stitch 

frameDiff = 7; % stitch together tracks if first ended 4 or fewer frames before second started
maxDist = 0.1; % stitch together tracks if first ended within 1 mm of second's start

buffer = [0.1 0.1 0.1 0.1];  % [left, right, top, bottom]
                             
trimrect = []; % [left bottom right top] you can specify an actual trim rectangle here, which will mean buffer is ignored
               % note that the window should be specified in cm, not pixels

trimrectpixels = []; % [left bottom right top] trim rectangle set in pixels; this is based on full image pixels (ROI is ignored); if both trimrect and trimrectpixels are specified, trimrect wins


               
setRulesByPeri = true; % instead of using default derivation rules, set smoothing, derivative times, and minimum run times by detected peristalsis frequency
               
%we derive all the fields necessary to segment the tracks and generate navigation figures before saving
%if there are other fields you want derived and saved, place them here
fieldsToDerive = {};

individualMatFiles = true;

fixedInterpTime = []; %if fixedInterpTime is a number, then all experiments will have same fixed interpolation time (not determined from frame rate)

camcalinfo = []; %allows you to directly pass in camcalinfo
requireMaggot = true;
%TODO: fix all non-essential displays to respond to verbosity trigger
verbose = false;

varargin = assignApplicable(varargin);


if verbose 
    disp('~~~~~~~~~ processJavaDirectoryToMatfiles ~~~~~~~~~');
end
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
mintime = 30; % tracks must be at least 60 seconds -- we will convert this into a minimum number of points later
ecl.minHTValid = 0.9; % discard any track where we couldn't determine a head tail properly at least 90% of the time
ecl.rpmCut = 2; % get rid of any track that circles in one direction with a frequency of at least 2 rev/min
ecl.minRevCut = 3; % provided that track does at least 3 complete circles


ccInSupDataDir = false;



varargin = assignApplicableAndFixTypes(varargin);
success = false;
ts1 = tic();




% this section will look through the directory for .jav files and verify
% that we have the appropriate .mdat files as well
d = [dir(fullfile(basedir, '*.jav'))];
d = d(~cellfun(@(n) strcmpi(n,'divergedTrackExp.jav'),{d.name}));

if (isempty(d))
    disp (['no .jav files in ' basedir]);
    eset = [];
    return;
end




matfileexists = false(size(d));
for j = 1:length(d)
    fname{j} = fullfile(basedir, d(j).name); %#ok<*AGROW>
    [ps,nm] = fileparts(fname{j});
    mdatname{j} = fullfile(ps,[nm '.mdat']);
    if (~exist(mdatname{j}, 'file'))
        disp ('cannot find metadata file');
        disp (['jav file name = ' fname{j}]);
        disp (['metadata file name = ' mdatname{j}]);
        error ('must have metadata file to proceed');
        return %#ok<UNRCH>
    end
    matfileexists(j) = exist(fullfile(ps, 'matfiles', [nm '.mat']), 'file');
    ptTp = Experiment.getJavaPtType(fname{j});
    ismaggot(j) = isa(ptTp, 'MaggotTrackPoint') || isa(ptTp, 'LarvaTrackPoint');
end
ismaggot = logical(ismaggot);
if (requireMaggot)
    if ~all(ismaggot)
        for j = find(~ismaggot)
            warning ('eset:pd2m', '%s', [fname{j} ' does not contain maggot track points & is being ignored']);
        end
    end
    fname = fname(ismaggot);
    mdatname = mdatname(ismaggot);
    matfileexists = matfileexists(ismaggot);
    
        
end


if (individualMatFiles)
    pfname = fname(matfileexists);
    fname = fname(~matfileexists);
    if (isempty(fname))
        eset = [];
        if (verbose)
            disp (['already processed all jav files in ' basedir]);
        end
        return;
    end
    for j = 1:length(pfname)
        if (verbose)
            disp (['already processed: ' pfname{j}]);
        end
    end
    for j = 1:length(fname)
        if (verbose)
            disp (['need to process: ' fname{j}]);
        end
    end
end

disp('**************************');
disp(['Processing jav->mat (' num2str(length(fname)) ' files) in dir: ' basedir]);
disp('**************************');

%%%%%%%%%%%%%%% CAMCAL STUFF
% this section carries out the camera calibration if needed and saves it
% to disk
disp('***** Carrying out copying camera calibration and/or copying calibration data');
if (ccInSupDataDir)
    for j = 1:length(fname) %#ok<UNRCH>
        [p,f] = fileparts(fname{j});
        if (~exist(fullfile(p, [f ' sup data dir'], 'camcalinfo.mat'), 'file'))
            warning ('es:pd2m', ['you told me to look in ' f ' sup data dir for camcalinfo, but I couldn''t find it']);
            return;
        end
        temp = load(fullfile(p, [f ' sup data dir'], 'camcalinfo.mat'), 'camcalinfo');
        if (j == 1)
            camcalinfo = temp.camcalinfo;
        else
            camcalinfo(j) = temp.camcalinfo;
        end
    end
else
    if (~isempty(camcalinfo) || ~exist(fullfile(basedir, 'camcalinfo.mat'), 'file'))
        if (isempty(camcalinfo))
            try
                if (~isempty(checkerLocation))
                    im = imread(fullfile(checkerLocation, checkerboardfname));
                else
                    im = imread(fullfile(basedir, checkerboardfname));
                end
            catch me
                disp ('could not read in checkerboard from: ');
                disp (fullfile(basedir, checkerboardfname));
                disp ('error was: ')
                disp (me.getReport());
                disp ('');
                warning ('cannot proceed without camera calibration -- aborting');
                eset = [];
                return;
            end
            if (verbose)
                disp ('calibrating checkerboard');
            end        
            camcalinfo = CameraCalibration(im, checkeroptions{:});
            if (isempty(camcalinfo))
                error ('camera calibration failed - check options');
            end
        end
        % make sure that the checkerboard that displays here is correct -- I am
        % not putting a checkoff in the script because you want to run it automatically, but this is a potential place for errors
        
        
        save(fullfile(basedir, 'camcalinfo.mat'), 'camcalinfo');
    else
        load(fullfile(basedir, 'camcalinfo.mat'), 'camcalinfo');
    end    
   close all
end  


% now we will load each mdat file and see if there is an ROI offset we need
% to take into account for the camera calibration

for j = 1:length(mdatname)
     ind = min(j, length(camcalinfo));
    if (isa(camcalinfo(ind), 'SimpleScalingCameraCalibration') || isa(camcalinfo(ind), 'ScalingAndRotationCameraCalibration'))
        cc(j) = camcalinfo(ind);
        continue;
    end
   
    rx = camcalinfo(ind).realx;
    ry = camcalinfo(ind).realy;
    cx = camcalinfo(ind).camx;
    cy = camcalinfo(ind).camy;
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



try
    if (isempty(trimrect) && ~isempty (trimrectpixels))
        [~,trimrect,~] = camcalinfo.realRectFromCamRect(trimrectpixels);
    end
catch me
    disp(me.getReport);
    return;
end


%%%%%%%CAMCALDONE

disp('***** Done carrying out copying camera calibration and/or copying calibration data');

% This section calls the processing functions
args = [varargin, {'frameDiff', frameDiff, 'maxDist', maxDist, 'setRulesByPeri', setRulesByPeri, 'fieldsToDerive', fieldsToDerive}];

if (individualMatFiles)
    if (loadParallel && length(fname) > 1)
        ts1 = tic;
        if (matlabpool('size') == 0)
            if (verbose)
                disp ('opening parallel processes'); 
            end
            matlabpool;
            toc(ts1);
            closepool = true;
        else
            closepool = false;
        end
        
    %%NOTE: when matlabpool is switched to parpool, uncomment lines below 
        
    %AddJavaPathToParPool(poolName) %% UMCOMMENT THIS
         
        parfor n=1:length(fname)
            
            %AddJavaPath(getAttachedFilesFolder) %% UMCOMMENT THIS
            
            if (verbose)
                disp(['Loading file #' num2str(n) ' : ' fname{n}]);
            end        
            try
                %2/7/2015 -- changed camcalinfo to cc(n): camcalinfo is
                %nonvectorized and does not include ROI offsets.
                ss(n) = processSingleJavaExperimentToMatfile (fname{n}, ecl, eclnukespots, cc(n), minptsToLoad, fixedInterpTime, trimrect, buffer, mintime, basedir, args{:});
                disp(['finished loading file #' num2str(n) ' : t = ' num2str(toc(ts1))]);
            catch e
                disp(['Error in processSingleJavaExperimentToMatfile for: ' fname{n}]);
            end
            
        end
        
        if (closepool)
            disp('closing parallel processes'); ts2 = tic;
            matlabpool close;
            toc(ts2);
        end
        disp(['total time to load ' num2str(length(fname)) ' files = ' num2str(toc(ts1))]);
        success = success && all(ss);
    else
        for j = 1:length(fname)
            try
                %2/7/2015 -- changed camcalinfo to cc(n): camcalinfo is
                %nonvectorized and does not include ROI offsets.
                success = success &  processSingleJavaExperimentToMatfile (fname{j}, ecl, eclnukespots, cc(j), minptsToLoad, fixedInterpTime, trimrect, buffer, mintime, basedir, args{:});
            catch me
                disp(['Error in processSingleJavaExperimentToMatfile for: ' fname{j}]);
                disp(me.getReport());
            end
        end
    end
    eset = [];
    
    
    return;

else
    % this command actually loads the experiments from disk -- so short!
    eset = ExperimentSet.fromJavaFiles(fname{:}, 'minpts', minptsToLoad, 'camcalinfo', cc, 'parallel', loadParallel, 'fixedInterpTime', fixedInterpTime);
    if (isempty(eset))
        return;
    end

    disp('***** Post-mat processing');

    % this is NEW NEW NEW -- we will try to determine appropriate smoothing
    % and interpolation times, as well as minimum run times based on
    % peristalsis detected in high frequency movies

    if (setRulesByPeri && median(eset.gatherSubField('dr', 'interpTime')) <= 0.1251) % >= 8 Hz
        [mf, f, ps] = eset.setDerivationRulesByPeristalsisFrequency();
        disp (['calculated peristalis frequency is: ' num2str(mf), ' Hz.']);
        figure();
        plot (f, ps); title(['calculated peristalis frequency is: ' num2str(mf), ' Hz.']);
        xlabel ('freq'); ylabel ('power');
    end

    % we are going to put the rest of the file in try/catch blocks -- that way
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
                trimrect = [ll(1) ll(2) ur(1) ur(2)] + buffer.*[1 1 -1 -1];
            end
        end
    catch me
        disp(me.getReport);
        return;
    end
    try
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

    disp('calculating dq');
    try
        eset.evaluateTrackExpression('track.so.autoset_curv_cut= true;'); %NEW set radius of curvature cutoff to be 1/5 of body length
        eset.executeTrackFunction('setSegmentSpeeds');
        dvfields = {'spineTheta', 'vel_dp', 'spineLength', 'sspineTheta', 'lrdtheta'};
        for j = 1:length(dvfields)
            disp (['calculating ' dvfields{j}]);
            eset.gatherField(dvfields{j});
        end


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
    disp('***** Done with post-mat processing');

    %finally, save the eset to disk as matfiles, so it can be loaded quickly later
    disp('***** Saving mat files');
    try
        if (~exist(fullfile(basedir, 'matfiles'), 'dir'))
            mkdir (fullfile(basedir, 'matfiles'));
        end
        eset.executeExperimentFunction('toMatFile', 'extradir', 'matfiles');
    %{
        if (individualMatFiles)
            eset.executeExperimentFunction('toMatFile', 'extradir', 'matfiles');
        else 
            eset.toMatFiles(fullfile(basedir, 'matfiles',esetname));
        end
        %}
        disp('saved file');
        toc(ts1)
    catch me
        disp(me.getReport);
        return;
    end

    disp('***** Done saving mat files');

    success = true;

    % needed to close all windows and run smoothly in the cluster
    close all;
    if verbose 
        disp('~~~~~~~~~ processJavaDirectoryToMatfiles: DONE ~~~~~~~~~');
    end
end


disp('**************************');
disp(['Done processing jav->mat in dir: ' basedir]);
disp('**************************');

function success = processSingleJavaExperimentToMatfile (fname, ecl, eclnukespots, camcalinfo, minptsToLoad, fixedInterpTime, trimrect, buffer, mintime, basedir, varargin)


frameDiff = 7; % stitch together tracks if first ended 4 or fewer frames before second started
maxDist = 0.1;              
setRulesByPeri = true; % instead of using default derivation rules, set smoothing, derivative times, and minimum run times by detected peristalsis frequency
               
%we derive all the fields necessary to segment the tracks and generate navigation figures before saving
%if there are other fields you want derived and saved, place them here
fieldsToDerive = {};
redomatfile = false;
tracksInSubDir = false;
varargin = assignApplicable(varargin);

success = false;
ts1 = tic();

% this section carries out the camera calibration if needed and saves it
% to disk

[ps,nm] = fileparts(fname);
mdatname = fullfile(ps,[nm '.mdat']);
if (~exist(mdatname, 'file'))
    disp ('cannot find metadata file');
    disp (['jav file name = ' fname]);
    disp (['metadata file name = ' mdatname]);
    error ('must have metadata file to proceed');
    return %#ok<UNRCH>
end
matfileexists = exist(fullfile(ps, 'matfiles', [nm '.mat']), 'file');

% THIS IS ALREADY DONE IN MAIN FUNCTION -- SKIP
% now we will load each mdat file and see if there is an ROI offset we need
% to take into account for the camera calibration
% 
% rx = camcalinfo.realx;
% ry = camcalinfo.realy;
% 
% cx = camcalinfo.camx;
% cy = camcalinfo.camy;
% data = importdata(mdatname);
% xcol = find(strcmpi(data.colheaders, 'ROIX'));
% if (~isempty(xcol))
%     roix = median(data.data(isfinite(data.data(:,xcol)),xcol)); % we should be able just to take the first value, but this should be more robust
% else
%     roix = 0;
% end
% ycol = find(strcmpi(data.colheaders, 'ROIY'));
% if (~isempty(ycol))
%     roiy = median(data.data(isfinite(data.data(:,ycol)),ycol)); % we should be able just to take the first value, but this should be more robust
% else
%     roiy = 0;
% end
% 
% cc = CameraCalibration(rx, ry, cx-roix, cy-roiy);

if (~redomatfile && matfileexists)
    success = true;
    return;
end




% this command actually loads the experiments from disk -- so short!
eset = ExperimentSet.fromJavaFiles(fname, 'minpts', minptsToLoad, 'camcalinfo', camcalinfo, 'parallel', false, 'fixedInterpTime', fixedInterpTime);
if (isempty(eset))
    return;
end

disp('***** Post-mat processing');

% this is NEW NEW NEW -- we will try to determine appropriate smoothing
% and interpolation times, as well as minimum run times based on
% peristalsis detected in high frequency movies
try
    disp ('setting derivation rules by peristalsis frequency');
    if (setRulesByPeri && median(eset.gatherSubField('dr', 'interpTime')) <= 0.1251) % >= 8 Hz
        [mf, f, ps] = eset.setDerivationRulesByPeristalsisFrequency();
     %   disp (['calculated peristalis frequency is: ' num2str(mf), ' Hz.']);
      %  figure();
       % plot (f, ps); title(['calculated peristalis frequency is: ' num2str(mf), ' Hz.']);
        %xlabel ('freq'); ylabel ('power');
    end
catch me
    disp (me.getReport());
    return;
end


% we are going to put the rest of the file in try/catch blocks -- that way
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
            trimrect = [ll(1) ll(2) ur(1) ur(2)] + buffer.*[1 1 -1 -1];
        end
    end
catch me
    disp(me.getReport);
    return;
end

try
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
try
    ecl.minPts = ceil(mintime / eset.expt(1).dr.interpTime);
    ecl.askFirst = false; 
    ecl.showFigsInReport = false;
    [~,rpt] = ecl.getReport(eset);
catch me
    disp (me.getReport())
    return;
end

[basedir,nm] = fileparts(fname);

try
    fid = fopen (fullfile (basedir, [nm ' cleaning report.txt']), 'wt');
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

if (isempty(eset.expt))
    return;
end

%calculate derived quantities used in segmentation; set segmentation speeds

try
    eset.evaluateTrackExpression('track.so.autoset_curv_cut= true;'); %NEW set radius of curvature cutoff to be 1/5 of body length
    eset.executeTrackFunction('setSegmentSpeeds');
    dvfields = {'spineTheta', 'vel_dp', 'spineLength', 'sspineTheta', 'lrdtheta'};
    for j = 1:length(dvfields)
        disp (['calculating ' dvfields{j}]);
        eset.gatherField(dvfields{j});
    end
    
        
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


disp('***** Done with post-mat processing');

disp('***** Saving mat files');

%finally, save the eset to disk as matfiles, so it can be loaded quickly later
try
    if (~exist(fullfile(basedir, 'matfiles'), 'dir'))
        mkdir (fullfile(basedir, 'matfiles'));
    end
    eset.executeExperimentFunction('toMatFile', 'extradir', 'matfiles', 'tracksInSubDir', tracksInSubDir);
    disp('saved file');
    imdir = fullfile(basedir, 'diagnostics');
    if (~exist(imdir, 'dir'))
        mkdir(imdir);
    end
    for j = 1:length(eset.expt)
        im = eset.expt(j).diagnosticImage();
        if (~isempty(im))
            [~,fn] = fileparts(eset.expt(j).fname);
            imwrite(im/255, fullfile(imdir, [fn ' diagnostic foreground.bmp']), 'bmp');
        end
    end
    
    toc(ts1)
catch me
    disp(me.getReport);
    return;
end

disp('***** Done saving mat files');
success = true;





