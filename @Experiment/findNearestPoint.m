function [pt, track, trackind, ptind] = findNearestPoint (expt, loc)
%searches through all points in tracks in expt to find point(s) nearest loc
%function [pt, track, trackind, ptind] = findNearestPoint (expt, loc)
%
%outputs:
%PT: an array of point objects that duplicates the point(s) in expt nearest
%   loc(s)
%TRACK: an array of track (handle) object(s) containing the nearest
%   point(s)
%TRACKIND: the indices of the track(s) containing the nearest point(s)
%PTIND: the indices of the nearest point(s) within track(s)
%   i.e. EXPT.track(TRACKIND) = TRACK
%   EXPT.track(TRACKIND(j)).pt(PTIND(j)) is the same as PT(j)
%inputs:
%EXPT: a member of the experiment class
%LOC: a 2xN list of points

if (~exist ('x', 'var') || isempty (x))
    [x,y] = getpts();   
    %use last point selected
    x = x(end);
    y = y(end);
else
    x = loc(1,end);
    y = loc(2,end);
end

pts = expt.gatherField('loc');
dist = (pts(1,:)-x).^2 + (pts(2,:)-y).^2;
[blah,I] = min(dist);
n = cumsum([expt.track.npts]);
trackind = find(n >= I, 1, 'first');
if (trackind > 1)
    ptind = I - n(trackind-1);
else
    ptind = I;
end
track = expt.track(trackind);
pt = track.pt(ptind);