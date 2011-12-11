function [pt, ind, dist] = nearestPoint(track, loc)
% finds the nearest point(s) in a track to given location(s)
% function [pt, ind, dist] = nearestPoint(track, loc)
%
% outputs:
%   PT: TrackPoint(s) that duplicate the point(s) closest to LOC(s)
%   IND: the index(ices) of the closest point(s) in TRACK
%   DIST: the distance(s) between LOC(s) and PT(s)
% inputs:
%   TRACK: a member of the track class
%   LOC: (optional) a 2xN location vector
%       if loc is empty, user is prompted for graphical input using getpts
updatept = false;
if (~exist('loc', 'var') || isempty(loc))
    [x,y] = getpts();
    loc = ([x y])';
    size(loc)
    updatept = true;
end

if (size(loc,2) > 1)
    for j = 1:size(loc,2)
        [p,i,d] = nearestPoint (track, loc(:,j));
        pt(j) = p;
        ind(j) = i;
        dist(j) = d;
    end
else 
    pt = [track.pt];
    [dist, ind] = min (sum(([pt.loc] - repmat(loc, [1 length(track.pt)])).^2, 1));
    pt = track.pt(ind);
    dist = sqrt(dist);
end

if (updatept)
    hh = ishold;
    hold on;
    lnew = [pt.loc];
    plot (loc(1,:), loc(2,:), 'gx', lnew(1,:), lnew(2,:), 'ro');
    if (~hh)
        hold off
    end
end