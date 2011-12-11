function af = fromFile (fname, timfname, loadContour, camcalinfo, minTrackLength)
%loads an experiment from a bin file
%expt = fromFile (fname, timfname, loadContour, camcalinfo, minTrackLength)
%this is a static method of the Experiment class (Experiment.fromFile)
%
%outputs: 
%EXPT, a member of the experiment class
%inputs:
%FNAME: name of .bin file to load
%TIMFNAME: timining information file (.tim);
%   default: change extension of fname to .tim
%LOADCONTOUR: whether to load the contour if this is a maggot track:
%   default TRUE
%CAMCALINFO: camera calibration struct (ask Marc); pass empty ([]) to ignore
%   default: []
%MINTRACKLENGTH: minimum length of a track (in points) to load from disk
%   default: 1
tic
if (~exist ('loadContour', 'var') || isempty ('loadContour'))
    loadContour = true;
end

if (~exist ('camcalinfo', 'var'))   
    camcalinfo = [];
end

if (~exist ('minTrackLength', 'var') || isempty(minTrackLength))
    minTrackLength = 1;
end

af = Experiment();

af.fname = fname;
af.camcalinfo = camcalinfo;
if (iscell(fname))
    [~,~,ext] = fileparts(fname{1});
else
    [~,~,ext] = fileparts(fname);
end

if (strcmpi (ext, '.blob') || strcmpi (ext, '.blobs'))
    af.track = MWTTrack.fromFile(fname, camcalinfo);
    return;
end
if (iscell(fname))
    error ('multiple file names only supported for MWT blob(s) files');
end


d = dir(fname);
totalsize = d.bytes;

fid = fopen(fname, 'r');
code = fread(fid, 1, 'int32');
disp (['code = ' num2str(code, '%x')]);

switch (bitshift(code, -16))
    case 1
        ptType = TrackPoint();
    case 2
        ptType = ImTrackPoint();
    case 3
        ptType = OldMaggotTrackPoint();
        af.so = MaggotSegmentOptions();
        %af.dr.smoothTime = 1;
    case 4
        ptType = MaggotTrackPoint();
        af.so = MaggotSegmentOptions();
    otherwise
      disp('invalid code: I don''t know what kind of point I''m loading');
      return
end  


ntracks = fread(fid, 1, 'int32');
if (isa(ptType, 'MaggotTrackPoint'))
    bob = repmat(MaggotTrack(), [1 ntracks]);
else
    bob = repmat(WormTrack(), [1 ntracks]);
end
%af.track = repmat(Track(), [1 ntracks]);
ts = tic;
lastelapsed = 0;
reportEvery = 60;
for j = 1:ntracks
    elapsed = toc(ts);
    if (elapsed - lastelapsed > reportEvery)
        lastelapsed = elapsed;
        disp ([num2str(elapsed) 's: ' num2str(ftell(fid)) '/' num2str(totalsize) ' bytes (' num2str(100*ftell(fid)/totalsize, 2) '%) loaded' ...
            num2str(elapsed*(totalsize - ftell(fid))/ftell(fid)) ' s remain']);
    end
    bob(j) = Track.fromFile(fid, ptType, [], loadContour, camcalinfo, minTrackLength);
    if (mod(j,100) == 0)
     %   toc(ts);
    end
end
%disp('now you have to wait, because matlab is slow and sucks nuts');
%ts = tic;

%modified by MHG 03/17 to show report
showReport = false;
if (showReport)
    semilogy (0:10:max([bob.npts]), hist([bob.npts], 0:10:max([bob.npts])), 'b-'); hold on;
    title ('distribution of track length in frames');
end
af.track = bob([bob.npts] > minTrackLength);

if (showReport)
    semilogy (0:10:max([af.track.npts]), hist([af.track.npts], 0:10:max([af.track.npts])), 'r-'); hold off;
    legend ('pre trim', 'post trim');
end
%toc(ts);

[af.track.expt] = deal(af);

fclose(fid);

disp('adding timing information'); 

try
    af.addtime(timfname);
catch me
    disp(me.getReport());
end
toc(ts)
