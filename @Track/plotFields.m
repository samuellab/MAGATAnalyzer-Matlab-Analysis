function h = plotFields (track, xfield, yfields, varargin)
%plots multiple yfields against one x field
%function plotFields (track, xfield, yfields, varargin)
%
%inputs:
%   track < Track
%   xfield: the field defining the x variable
%   yfield: the fields defining the y variables
%outputs:
%   h < graphics handle
%optional arguments:
%   'inds', inds to plot (default all)
%   'Axes', axes to plot in (default gca)
%   'makeLegend', whether to display the legend (default true)
%   anything that can be passed to plot
%   'labeled', if true, label according to segmentation
%   
%   if length(xfield) = length(yfields) = length(Axes), plots yfields{j}
%   vs. xfield{j} in Axes(j)
Axes = [];
varargin = assignApplicable(varargin);
if (iscell(xfield) && iscell(yfield) && length(xfield) == length(yfields) && length(xfield) == length(Axes) && length(Axes) > 1)
    for j = 1:length(Axes)
        plotFields(track, xfield{j}, yfields{j}, 'Axes', Axes(j), varargin{:});
    end
    return
end
inds = [];
makeLegend = true;
labeled = false;
varargin = assignApplicable(varargin);

if (isempty(Axes))
    Axes = gca;
end


if (isempty(inds))
    x = track.getDerivedQuantity(xfield);
    inds = 1:length(x);
else
    x = track.getDerivedQuantity(xfield,false,inds);
end
y = [];
if (~iscell(yfields))
    yfields = {yfields};
end
leg = {};
for j = 1:length(yfields)
    yf = track.getDerivedQuantity(yfields{j},false,inds);
    if size(yf,1) == 1
        leg = [leg, yfields{j}];
    else
        for k = 1:size(yf,1)
            leg = [leg, [yfields{j} '(' num2str(k) ')']];
        end
    end
    y = [y;yf];
end

hh = plot (Axes, x, y, varargin{:});

hl = [];
lleg = {};
if (labeled && ~isempty(track.run) && isa(track.run, 'Run') && any(track.isrun(inds)))
    hold(Axes,'on');
    hlt = plot(Axes, x(track.isrun(inds)), y(:, track.isrun(inds)),'m.','MarkerSize',5);
    hl = [hl hlt(1)];
    lleg = [lleg, 'run'];
    hold(Axes,'off');
end
if (labeled && ~isa(track, 'MaggotTrack') && ~isempty(track.sharpTurn) && isa(track.sharpTurn, 'WormSharpTurn'))
    hold(Axes,'on');
    st = track.sharpTurn;
    mcolors = {'c*', 'yo', 'rh', 'gh'};
    for j = -1:1:2
        if (~any([st.typeCode] == j))
            continue;
        end
        [~,indsom] = intersect(inds, [st([st.typeCode] == j).inds]);
        if (~isempty(indsom))
            hlt = plot(Axes, x(indsom), y(:, indsom),mcolors{j+2},'MarkerSize',5);
            hl = [hl hlt(1)];
            lleg = [lleg, st(find([st.typeCode]==j,1)).type];
        end
    end
    hold(Axes,'off');
end
    
if (makeLegend)
    xlabel(Axes, xfield);
    if (~isempty(hl));
        legend(Axes, [hh; hl'], [leg, lleg]);
    else
        legend(Axes, hh,leg);
    end
end
if (nargout > 0)
    h = hh;
end