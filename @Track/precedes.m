function c = precedes(track1, track2, maxFrameDiff, maxDist)
% true if track1 precedes track2 in a stitch tracks sense
% function c = precedes(track1, track2, maxFrameDiff, maxDist)
% utility function:
% returns true iff track2 starts within maxFrameDiff after the end of track1
% AND the distance between the start of track2 and the end of track 1 <
% maxDist
%
% outputs: 
%   C: TRUE OR FALSE 
% inputs:
%   TRACK1, TRACK2 < Track
%   MAXFRAMEDIFF: maximum number of frames between end of track1 & start of
%       track2 for C to be true
%   MAXDIST: maximum distance between the end of track1 & start of track 2
%       for C to be true

c = (track1.pt(end).ind < track2.pt(1).ind && track1.pt(end).ind + maxFrameDiff >= track2.pt(1).ind ...
        && track1.pt(end).distance(track2.pt(1)) < maxDist);
