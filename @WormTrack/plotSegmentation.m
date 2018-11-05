function plotSegmentation(track, varargin)
% plots the segments of the track (runs & reorientations) graphically
% function plotSegmentation(track, varargin)
%
% outputs: none
% inputs:
%   TRACK < Track;
%   VARARGIN: optional parameter/value pairs
%       'Axes', ax : which axes to plot to, default gca
%       'fieldName', fieldName : what field to plot, default 'sloc'
%       'stMarkerSize', pts : size in points of the markers for sharp turns,
%           default 10
%       'multicolor', tf : default true, whether to plot each track in a
%           different color
%       anything that can be passed to PLOT, affects plotting of runs only

Axes = [];
fieldName = 'sloc';
stMarkerSize = 10;
multicolor = true;
varargin = assignApplicable(varargin);
if (isempty(Axes))
    Axes = gca;
end

ih = ishold(Axes);

for j = 1:length(track.run)
    if (multicolor)
        col = nthcolor(j);
    else
        col = 'b-';
    end
    track.plotPath(fieldName, col, 'inds', track.run(j).inds, 'Axes', Axes, varargin{:});
    hold(Axes, 'on');
end

axes(Axes);
for j = 1:length(track.sharpTurn)
    switch(track.sharpTurn(j).typeCode)
        case -1
            letter = '\Omega';
            col = 'm';
        case 0
            letter = 'b';
            col = 'g';
        case 1
            letter = 'R';
            col = 'r';
        case 2
            letter = 'R2';
            col = 'g';
    end
    loc = double(track.dq.sloc(:,track.sharpTurn(j).centralInd));
    text(loc(1), loc(2), letter, 'Color', col, 'HorizontalAlignment', 'Center');
%    track.plotPath(fieldName, mc, 'inds', track.sharpTurn(j).centralInd, 'Axes', Axes, 'MarkerSize', stMarkerSize);
end

for j = 1:length(track.reorientation)
    track.plotPath(fieldName, 'k-', 'inds', track.reorientation(j).inds,'Axes',Axes,'Color', [0.7 0.7 0.7]);
end
if (~ih)
    hold(Axes,'off');
end
%for j = 1:length(
