function plotPath(track, pathType, linetype, varargin)
% plots the track path
% function plotPath(track, pathType, linetype, varargin)
%
% outputs: none
% inputs:
%   TRACK < Track; a single track, or a list of tracks.  if a list, same
%       as calling track(j).plotPath(...) on each track with hold on
%   PATHTYPE: the field that gives the location: default 'sloc'
%   LINETYPE: a line type specifier, see PLOT.  default 'b-'
%   VARARGIN: optional paramter/value pairs:
%       'inds', inds
%       'highlightinds', inds
%           'highlightinds' can also be a field of track, e.g.
%           'iscollision'
%       'highlightlinetype', lt (default 'r.')
%           plots with highlightlinetype at highlight inds on top of line
%       'Axes', Axes -- which axes to plot in
%       any parameter/value pair that can be passed to PLOT
%       'indsExpression', an expression evaluated to provide inds
%                         e.g. 'track.getDerivedQuantity(''eti'') < 10'
%       'highlightindsExpression' 
Axes = [];
plotEnds = false;
varargin = assignApplicable (varargin);
if (isempty(Axes))
    Axes = gca;
end 

if (length(track) > 1)
    ih = ishold(Axes);
    existsAndDefault('pathType', 'sloc');
    existsAndDefault('linetype', 'b-');
    for j = 1:length(track)
        track(j).plotPath(pathType, linetype, 'Axes', Axes, varargin{:});
        hold (Axes,'on');
    end
    if (~ih)
        hold (Axes,'off');
    end
    return;
end

if (~exist ('pathType', 'var') || isempty(pathType))
    pathType = 'sloc';
end

x = track.getDerivedQuantity(pathType);

inds = 1:length(x);

highlightinds = [];
highlightlinetype = 'r.';

if (~exist ('linetype', 'var') || isempty(linetype))
    linetype = 'b-';
end
indsExpression = [];
highlightindsExpression = [];
varargin = assignApplicable (varargin);
if (~isempty(indsExpression))
    inds = eval(indsExpression);
end
if (~isempty(highlightindsExpression))
    highlightinds = eval(highlightindsExpression);
end

if (ischar(highlightinds) && any(strcmp(properties(track), highlightinds)))
    highlightinds = track.(highlightinds);
    if (length(highlightinds) == length(x))
        highlightinds = logical(highlightinds);
    end
    
end

%{
if (~isfield (track.dq, 'sloc'))
    x = [track.pt.loc];
else
    x = track.dq.sloc;
end
%}
inds = inds(inds > 0 & inds <= length(track.dq.eti));

if (isempty(highlightinds))
    plot (Axes, x(1,inds), x(2,inds),linetype,varargin{:});
else
    plot (Axes, x(1,inds), x(2,inds),linetype, x(1,highlightinds), x(2,highlightinds), highlightlinetype, varargin{:});
end

if (plotEnds)
    try
        color = linetype(1);    
        hold on;
        scatter(Axes, x(1,1), x(2,1), [color 'o']);                  
        scatter(Axes, x(1,end), x(2,end), [color 'o']);
        scatter(Axes, x(1,end), x(2,end), [color 'x']);
        hold off;
    catch e
        hold on;
        scatter(Axes, x(1,1), x(2,1), 'bo');                  
        scatter(Axes, x(1,end), x(2,end), 'bo');
        scatter(Axes, x(1,end), x(2,end), 'bx');
        hold off;
    end
end
