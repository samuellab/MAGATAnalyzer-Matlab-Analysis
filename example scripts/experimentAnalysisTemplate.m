%this is the directory where your .bin and .tim files are stored
basedir = 'd:\marc Processed\maggots\ethyl acetate 4 pct 20 2000\';

minLoadPts = 50; %don't load any tracks shorter than this

%load all the files, this will take a while, so we only do it
%if we haven't already loaded them
if (~exist('eset', 'var'))
    eset = ExperimentSet.fromFiles(basedir,'minpts',minLoadPts); 
    disp('files loaded; fixing head tail orientation; this may take a while');
    tic
    eset.executeTrackFunction('fixHTOrientation');
    toc
end

existsAndDefault('restitchAndTrim',true);
stitchDist = 10; %distance that tracks can be apart to be stitched
stitchFrameDiff = 5; %number of frames that can be between end and next start when stitching
minTrackDist = 75; %minimum distance a track must travel in order to survive trimming
minTrackLength = 400; %after stitching, minimum number of frames a track must have 
minTrackSpeed = 0; %minimum average speed to avoid trimming

if (restitchAndTrim)
    disp('stitching tracks and trimming slow and short ones'); ts = tic(); 
    ntstart = length([eset.expt.track]);
    eset.executeExperimentFunction('stitchTracks', stitchFrameDiff, stitchDist);
    npoststitch = length([eset.expt.track]);
    for j = 1:length(eset.expt)
       meanspeed = eset.expt(j).evaluateTrackExpression('mean(track.getDerivedQuantity(''speed''))');
       dt = sqrt(eset.expt(j).evaluateTrackExpression('max(sum(track.getDerivedQuantity(''displacement'').^2))'));
       valid = [eset.expt(j).track.npts] > minTrackLength & meanspeed > minTrackSpeed & dt > minTrackDist;
       eset.expt(j).track = eset.expt(j).track(valid);
    end
    nposttrim = length([eset.expt.track]);
    disp(['num tracks start: ' num2str(ntstart) ', after stitch: ' num2str(npoststitch) ' after trim: ' num2str(nposttrim)]);
    restitchAndTrim = false;
    toc(ts)
end

%set segmentation options here
so = MaggotSegmentOptions();
so.curv_cut = 1;

existsAndDefault('resegment', true);
if (resegment)
    %synch up segment options
    [eset.expt.so] = deal(so);
    for j = 1:length(eset.expt)
        [eset.expt(j).track.so] = deal(eset.expt(j).so);
    end
%    eset.executeExperimentFunction('[expt.track.so] = deal(expt.so)');
    
    %recommend running setSegmentSpeeds to let track set its own segment
    %speed based on high curvature regions
    eset.executeTrackFunction('setSegmentSpeeds');
    
    disp('segmenting tracks, this can take a while'); tic
    eset.executeTrackFunction('segmentTrack');
    toc
    resegment = false;
end

varsToClean = {'minLoadPts', 'stitchDist','stitchFrameDiff','minTrackDist = 75','minTrackLength = 400','minTrackSpeed',...
    'ntstart', 'npoststitch', 'meanspeed', 'dt', 'valid', 'nposttrim', 'so', 'varsToClean', 'ts'};
clear(varsToClean{:});

existsAndDefault('saveMe', false);
if (saveMe)
    disp('saving to disk, this can take an ass long time');
    tic
    save ([basedir 'importedToMatlab.mat']);
    toc
    saveMe = false;
end

