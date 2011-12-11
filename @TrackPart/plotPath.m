function plotPath(tp, pathType, linetype, varargin)
%function plotPath(tp, pathType, linetype, varargin)
%
%optional args, in addition to those passed to Track.plotPath
%'ptbuffer', n (default 10): number of points before and after inds to plot
%

Axes = [];
varargin = assignApplicable (varargin);
if (isempty(Axes))
    Axes = gca;
end 

existsAndDefault('linetype', 'k-');
existsAndDefault('pathType', 'sloc');
    
if (length(tp) > 1)
    ih = ishold(Axes);
    for j = 1:length(tp)
        tp(j).plotPath(pathType, linetype, 'Axes', Axes, varargin{:});
        hold (Axes,'on');
    end
    if (~ih)
        hold (Axes,'off');
    end
    return;
end
ptbuffer = 10;
highlightlinetype = 'r-';
varargin = assignApplicable(varargin);

si = max(tp.startInd - ptbuffer, 1);
ei = min(tp.endInd + ptbuffer, length(tp.track.getDerivedQuantity('eti')));

tp.track.plotPath(pathType, linetype, 'inds', si:ei, 'highlightinds', tp.startInd:tp.endInd, 'highlightlinetype', highlightlinetype, ...
    'Axes', Axes, varargin{:});
