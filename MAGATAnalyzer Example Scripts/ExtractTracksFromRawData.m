datadirectory = fullfile ('..', 'Example Data');
imagedirectory = fullfile (datadirectory, 'Image Data');
extracteddirectory = fullfile(datadirectory, 'Extracted Tracks');
dataname = 'CS EtAc 2ppm per cm spatial gradient';
zippeddata = fullfile (datadirectory, [dataname '.zip']);


%This section unzips the raw data so it can be extracted by the track
%extraction software
if (~exist (imagedirectory, 'dir'))
    disp ('extracting example mmf data files from compressed zip archive');
    disp ('you should only see this message once');
    
    disp ('on completion, there should be 3 mmf files in ' );
    try
        file = java.io.File(fullfile(cd, imagedirectory));
        disp (char(file.getCanonicalPath));
    catch me
        disp (imagedirectory);
    end
    try
        unzip (zippeddata, imagedirectory);
    catch me
        disp ('unzipping files failed, with error: ');
        disp(me.getReport());
        disp ('');
        disp (['please delete the image directory ' imagedirectory ' and try again']);
        return;
    end
    disp ('unzipping finished!');
end

%This section creates the extraction settings prior to calling the track
%extractor

be = defaultExtractionProcessingParams;
f2p = be.files_to_process;
pp = f2p.processing_params;

% these are fields you may need to adjust when using the track extractor on
% your own data

pp.verbosity_level = 1; %only show errors in log file
pp.analysis_rectangle = [0 0 2592 1944]; %region of the image to analyze -- set to exclude edges etc.

% these parameters determine whether a spot that is brighter than the
% background might represent a larva and therefore should be added to a
% track
pp.minArea = 20; %minimum area (in pixels) a spot must be in order to be added to a track
pp.maxArea = 1600; %maximum area (in pixels) a spot may be in order to be added to a track
pp.overallThreshold = 30; %threshold above background for a spot to be added to a track

% as we select bright moving spots from the background and add them to
% tracks, we save a square window of the image around that spot.  winSize
% must be >= the longest maggot in the data set.  with that constraint,
% winSize should be as small as possible 
pp.winSize = 40;

% these parameters determine how the background is determined
% the background is the minimum of at least nBackgroundFrames 
% (when extracting from mmfs, the background computation is different,
% resulting in the inclusion of more than just nBackgroundFrames)
% the nBackgroundFrames evenly span the background_resampleInterval
% the frame being examined will be somewhere within the
% background_resampleInterval, but not necessarily in its center
pp.nBackgroundFrames = 5;
pp.background_resampleInterval = 500;

% maxExtractDist - the maximum distance between the centers of mass (in
% pixels) of spots in two successive frames that will allow those spots to
% be considered part of the same track.  
pp.maxExtractDist = 15;

% showExtraction - whether to show a window displaying the track extraction
% in progress.  if true (1), a window shows the background subtracted
% image; blue lines indicate active tracks; red lines indicate tracks that
% have ended (due to collisions, maggot leaving the analysis rectangle, dropping below
% threshold intensity etc.)
%
% after all the tracks are extracted, a window will pop up showing selected
% maggots with the head tail midline and contour overlaid.  

pp.showExtraction = true;
pp.extension = 'mmf'; %verify that we are set to process mmf files


% this code section finds the raw data files, creates an output
% directory to store the extracted tracks, sets the batch extraction
% settings, and saves the batch extraction settings to disk

srcdir = fullfile(imagedirectory, dataname);
d = dir (fullfile(srcdir, '*.mmf'));
if (isempty(d))
    disp ('expected to find mmf files in:')
    file = java.io.File(fullfile(cd, srcdir));
    disp (char(file.getCanonicalPath));
    disp ('but they were not found');
    return;
end

dstdir = fullfile(extracteddirectory, dataname);
if (~exist(dstdir, 'dir'))
    mkdir(dstdir);
end


f2p.processing_params = pp;

f2p = repmat(f2p, size(d));
for j = 1:length(d)
    [~,fn] = fileparts(d(j).name);
    
    file = java.io.File(fullfile(cd, srcdir, fn));
    f2p(j).file_stub = char(file.getCanonicalPath);
    fn = regexprep(fn, '_stack$', '');
    file = java.io.File(fullfile(cd, dstdir, [fn '.bin']));
    f2p(j).output_file = char(file.getCanonicalPath);
    
end

be.files_to_process = f2p;
str = batchExtractorToString(be);

file = java.io.File(fullfile(cd,'..', 'Track-Extraction-Software', 'Command Line Program', 'Command Line Extractor', 'extract-stack.exe'));
progname = char(file.getCanonicalPath);

file = java.io.File(fullfile(cd, dstdir, 'extraction_settings.bxx'));
bxxname = char(file.getCanonicalPath);
fid = fopen(bxxname, 'wt');
fprintf(fid, '%s', str);
fclose(fid);

% copy the checkerboard image into the extracted directory
copyfile(fullfile(datadirectory, '1cm_checkerboard.bmp'), fullfile(dstdir, '1cm_checkerboard.bmp'));


disp ('running the track extraction.  this will take about 45 minutes');
disp ('when the extractor has finished, there will be 5 files for each experiment in this directory:')
file = java.io.File(fullfile(cd, dstdir));
disp(char(file.getCanonicalPath));
disp ('xxx_log.txt - messages from the extraction software');
disp ('xxx_header.txt - a text header describing how data is stored in the bin files');
disp ('xxx_foreground.bmp - a simulated long exposure image of the experiment');
disp ('xxx.bin - the extracted track data');
disp ('xxx.mdat - associated meta data including camera timing information');
disp ('when the track extraction software has finished running, this script will continue');

[s,w] = dos (['"' progname '" "' bxxname '"']);

disp ('extraction completed!');