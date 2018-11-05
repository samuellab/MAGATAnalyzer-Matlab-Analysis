function track2 = reloadTrack(expt, track, varargin)
%loads a copy of track from disk, or loads images into an existing track
%function track2 = reloadTrack(expt, track, varargin)
%
%outputs:
%TRACK2: the track, as loaded from disk
%inputs:
%EXPT: a member of the Experiment class
%TRACK: either a handle to a track, or the index of a track in EXPT.track
%   if track is a track handle, then 
%   creates a new track handle and reloads track into new handle (with all
%   images).  NB: you must delete track2 handle when done
%
%   if track is an index to a track, then we reload the track into the old
%   place in memory
%
% optional parameter: 
%   'usecamcalinfo', [true]/false
%    whether to use camera calibration when reloading file
%   use camcalinfo is always true if track is a number (to prevent mixed
%   representation in experiment)

usecamcalinfo = true;
varargin = assignApplicable(varargin);

if (~isa (track, 'Track')) 
    trackid = track;
    track = expt.track(trackid);
    usecamcalinfo = true;
else
    trackid = 0;
end

try 
    invalid = (expt.fid == 0 || ftell(expt.fid) < 0);
catch
    expt.openDataFile();
    invalid = (expt.fid == 0 || ftell(expt.fid) < 0);
end
if invalid
    expt.openDataFile();
end

if (usecamcalinfo)
    cc = expt.camcalinfo;
else
    cc = [];
end

if (isa(track, 'MaggotTrack'))   
    track2 = repmat(MaggotTrack(), [1 track.nt]);
else
    track2 = repmat(Track(), [1 track.nt]);
end
for j = 1:track.nt
    fseek(expt.fid, track.locInFile(j), -1);
    track2(j) = Track.fromFile(expt.fid, track.pt(1), true, true, cc);
end
for j = (track.nt-1):-1:1
    track2(j).merge(track2(j+1));
end
delete (track2(2:end));
track2 = track2(1);
track2.dr = track.dr;
if (~isempty(expt.elapsedTime))    
    indx = (1:length(expt.elapsedTime)) - 1;
    track2.addTime (indx, expt.elapsedTime);
end
track2.expt = track.expt;

if (trackid > 0)
    %take care of case where track may have been trimmed previously
    startind = find([track2.pt.ind] == expt.track(trackid).startFrame);
    endind = find([track2.pt.ind] == expt.track(trackid).endFrame);
    [expt.track(trackid).pt.imData] = track2.pt(startind:endind).imData;
end