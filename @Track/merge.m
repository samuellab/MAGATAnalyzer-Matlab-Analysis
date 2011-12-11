function merge(track, track2)
% merges two tracks into a single track; does not change track2
% function merge(track, track2)
%
% this utility function is called by stitchTracks and not by end user
%
% outputs: none
% inputs:
%   TRACK: first track; modified to contain track and track2
%   TRACK2: second track; unchanged, but should be deleted later
%
% note after calling this function, derived quantities, segmentation, etc. are
% invalid

track.pt = [track.pt track2.pt];
track.npts = track.npts + track2.npts;
track.nt = track.nt + track2.nt;
track.locInFile = [track.locInFile track2.locInFile];
track.endFrame = track2.endFrame;

track.recalculateDerivedQuantities;
    
