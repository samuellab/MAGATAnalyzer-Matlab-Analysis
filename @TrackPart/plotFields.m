function h = plotFields (tp, xfield, yfields, varargin)
% plots relevant fields from track part 
% see Track.plotFields
% tp < TrackPart
% xfield < field that forms x axis
% yfields < field(s) to plot vs. xfield
% varargin:
%   'Axes', axes to plot in (default gca)
%   'makeLegend', whether to display the legend (default true)
%   anything that can be passed to plot
%   'ptbuffer' = 30: number of points on either side to plot

Axes = [];
ptbuffer = 30;
varargin = assignApplicable(varargin);
if (isempty(Axes))
    Axes = gca;
end



si = max(tp.startInd - ptbuffer, 1);
ei = min(tp.endInd + ptbuffer, length(tp.track.getDerivedQuantity('eti')));

ih = ishold (Axes);
h{1} = tp.track.plotFields(xfield, yfields, 'inds', si:ei, 'LineWidth', 1, 'Axes', Axes, varargin{:});
hold (Axes, 'on');
h{2} = tp.track.plotFields(xfield, yfields, 'inds', tp.startInd:tp.endInd, 'LineWidth', 3, 'Axes', Axes,  varargin{:}, 'makeLegend', false);
if (~ih)
    hold(Axes, 'off');
end