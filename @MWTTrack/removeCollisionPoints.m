function valid = removeCollisionPoints(track, maxArea, varargin)
% maintains only the portion of the track where the contour area is below
% maxArea
% function valid = removeCollisionPoints(track, maxArea, varargin)

valid = true;
pt = [track.pt];
area = [pt.area];
if (mean(area) > maxArea)
    valid = false;
    return;
end
good = (area/mean(area) < 1.5);
if (all(good))   
    return;
end

start = find (diff(good) > 0) + 1;
if (good(1))
    start = [1 start];
end
stop = find(diff(good) < 0);
if (good(end))
    stop = [stop length(good)];
end

len = stop - start;
[~,I] = max(len);
track.pt = pt(start(I):stop(I));
track.npts = length(track.pt);
track.dq = [];


