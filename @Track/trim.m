function deleteMe = trim (track, timerange, validrect)
% removes excess points and determines if any are left
% function deleteMe = trim (track, timerange, validrect)
%
% removes all points from the track outside timerange and 
% all points from the first time track location leaves validrect
% until the end of the track 
%
% if this would result in no points remaining in the track, we return
% true; trimming of points may or may not have already occurred
% otherwise, we return false
%
% if timerange or validrect is empty, no comparison is performed
% validrect has format [x0 y0 x1 y1]
%
% outputs:
%   DELETEME: true if trimming would remove all points
% inputs:
%   TRACK < Track
%   TIMERANGE: [t0 t1]; removes all points outside this range of times
%       if empty, this step is skipped
%   VALIDRECT: [x0 y0 x1 y1]; removes all points from the first time track
%       leaves VALIDRECT to end of track; if empty, this step is skipped

if (exist ('timerange', 'var') && ~isempty(timerange))
    pt = [track.pt];
    validTime = [pt.et] >= min(timerange) & [pt.et] <= max(timerange);
    if (~any(validTime))
        deleteMe = true;
        return
    end
    track.pt = track.pt(validTime);
    track.npts = length(track.pt);
end

if (exist ('validrect', 'var') && ~isempty(validrect))
    pt = [track.pt];
    x = [pt.loc];
    validLoc = (x(1,:) >= validrect(1) & x(1,:) <= validrect(3) & x(2,:) >= validrect(2) & x(2,:) <= validrect(4));
    ind = find(~validLoc, 1, 'first');
    if (~isempty(ind))
        if (ind == 1)
            deleteMe = true;
            return;
        end
        track.pt = track.pt(1:(ind-1));
    end
end
track.npts = length(track.pt);
track.startFrame = track.pt(1).ind;
track.endFrame = track.pt(end).ind;
deleteMe = false;
