function addTime(track, indx, elapsedTime)
%adds timing information to a track
%function addTime(track, indx, elapsedTime)
%track.addTime(indx, elapsedTime)
%
%inputs:
%TRACK: a member of the track class
%INDX: a list of frame indices
%ELAPSEDTIME: a list of times in seconds
%frame indx(j) occurred at elapsedTime(j)
if (isempty(elapsedTime))
    return;
end

pt = [track.pt];
et = interp1(indx, elapsedTime, double([pt.ind]),'linear', NaN);
pt = pt(isfinite(et));
et = num2cell(et(isfinite(et)));
[pt.et] = et{:};
track.pt = pt;
track.npts = length(track.pt);
track.startFrame = pt(1).ind;
track.endFrame = pt(end).ind;