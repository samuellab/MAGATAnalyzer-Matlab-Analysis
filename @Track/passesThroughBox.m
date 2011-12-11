function tf = passesThroughBox (track, rect, timeInterval, varargin)
% whether a track passed through a region of space in a given time window
% function tf = passesThroughBox (track, rect, timeInterval, varargin)
%
% outputs: TF(j) - whether track(j) passes through rect in timeInterval
% inputs: TRACK < Track
%         RECT  - [xmin xmax ymin ymax]
%         TIMEINTERVAL - [tmin tmax]
% optional:
%         'locField' - 'iloc' default

if (length(track) > 1)
    tf = false(size(track));
    for j = 1:length(tf)
        tf(j) = track(j).passesThroughBox(rect, timeInterval, varargin{:});
    end
    return;
end

locField = 'iloc';
varargin = assignApplicable(varargin);

l = track.getDerivedQuantity(locField);
t = track.getDerivedQuantity('eti');

tf = any(l(1,:) >= rect(1) & l(1,:) <= rect(2) & l(2,:) >= rect(3) & l(2,:) <= rect(4)...
    & t >= timeInterval(1) & t <= timeInterval(2));
